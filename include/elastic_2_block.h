#ifndef ELASTIC_2_BLOCK_H
#define ELASTIC_2_BLOCK_H

#include <deal.II/base/convergence_table.h>
#include <deal.II/base/function.h>
#include <deal.II/base/logstream.h>
#include <deal.II/base/quadrature_lib.h>
#include <deal.II/base/timer.h>
#include <deal.II/base/utilities.h>
#include <deal.II/dofs/block_info.h>
#include <deal.II/dofs/dof_accessor.h>
#include <deal.II/dofs/dof_handler.h>
#include <deal.II/dofs/dof_renumbering.h>
#include <deal.II/dofs/dof_tools.h>
#include <deal.II/fe/fe_q.h>
#include <deal.II/fe/fe_system.h>
#include <deal.II/fe/fe_values.h>
#include <deal.II/grid/grid_generator.h>
#include <deal.II/grid/tria.h>
#include <deal.II/grid/tria_accessor.h>
#include <deal.II/grid/tria_boundary_lib.h>
#include <deal.II/grid/tria_iterator.h>
#include <deal.II/lac/block_sparse_matrix.h>
#include <deal.II/lac/block_vector.h>
#include <deal.II/lac/constraint_matrix.h>
#include <deal.II/lac/full_matrix.h>
#include <deal.II/lac/precondition.h>
#include <deal.II/lac/solver_cg.h>
#include <deal.II/lac/solver_gmres.h>
#include <deal.II/lac/sparse_direct.h>
#include <deal.II/lac/sparse_ilu.h>
#include <deal.II/lac/trilinos_block_sparse_matrix.h>
#include <deal.II/lac/trilinos_block_vector.h>
#include <deal.II/lac/trilinos_precondition.h>
#include <deal.II/lac/trilinos_sparse_matrix.h>
#include <deal.II/lac/trilinos_vector.h>
#include <deal.II/numerics/data_out.h>
#include <deal.II/numerics/data_out_faces.h>
#include <deal.II/numerics/error_estimator.h>
#include <deal.II/numerics/matrix_tools.h>
#include <deal.II/numerics/vector_tools.h>

#include <fstream>
#include <sstream>
#include <string>
#include <typeinfo>

#include "elastic_base.h"
#include "parameters.h"
#include "preconditioner_2block.h"

using namespace std;
using namespace dealii;
namespace Elastic
{
template <int dim>
class Elastic2Blocks : public ElasticBase<dim> {
public:
    Elastic2Blocks (const unsigned int degree, const int _info);

private:
    // Setup Algebraic multigrid(AMG)
    virtual void setup_AMG ();
    // Solve the system
    virtual void solve ();

    std_cxx1x::shared_ptr<typename Preconditioner::inner> A_preconditioner;
    std_cxx1x::shared_ptr<typename Preconditioner::schur> S_preconditioner;
};
}

/*
 ------------- IMPLEMENTATION --------------
 */
template <int dim>
Elastic::Elastic2Blocks<dim>::Elastic2Blocks (const unsigned int degree, const int _info)
    : Elastic::ElasticBase<dim>(degree, _info, 2){}

/**
 * Setup preconditioners
 */
template <int dim>
void
Elastic::Elastic2Blocks<dim>::setup_AMG ()
{
    A_preconditioner
            = std_cxx1x::shared_ptr<typename Preconditioner::inner>(new typename Preconditioner::inner());

    S_preconditioner
            = std_cxx1x::shared_ptr<typename Preconditioner::schur>(new typename Preconditioner::schur());

    std::vector<std::vector<bool> > constant_modes;
    std::vector<bool>  displacement_components (ElasticBase<dim>::n_components,true);

    // A
    //    displacement_components[0] = true;
    displacement_components[ElasticBase<dim>::n_components-1] = false;
    DoFTools::extract_constant_modes (ElasticBase<dim>::dof_handler,
                                      displacement_components,
                                      constant_modes);

    TrilinosWrappers::PreconditionAMG::AdditionalData amg_A;

    amg_A.constant_modes = constant_modes;
    amg_A.elliptic = true;
    amg_A.higher_order_elements = false;
    amg_A.smoother_sweeps = 2;
    amg_A.aggregation_threshold = ElasticBase<dim>::par->threshold;
//    amg_A.output_details = true;


    A_preconditioner->initialize( ElasticBase<dim>::system_preconditioner.block(0,0),amg_A);
    //
    TrilinosWrappers::PreconditionAMG::AdditionalData amg_S;

    amg_S.elliptic = true;
    amg_S.higher_order_elements = false;
    amg_S.smoother_sweeps = 2;
    amg_S.aggregation_threshold = ElasticBase<dim>::par->threshold;

    // elem-by-elem Schur
    S_preconditioner->initialize( ElasticBase<dim>::system_preconditioner.block(1,1), amg_S);
}

template <int dim>
void
Elastic::Elastic2Blocks<dim>::solve ()
{
    const Preconditioner2Blocks< typename Preconditioner::inner, // A, schur
            typename Preconditioner::schur>
            preconditioner( ElasticBase<dim>::system_preconditioner, *A_preconditioner, *S_preconditioner); // system_matrix

    SolverControl solver_control (ElasticBase<dim>::system_matrix.m(),
                                  ElasticBase<dim>::par->TOL*ElasticBase<dim>::system_rhs.l2_norm());

#ifdef LOGRUN
    solver_control.enable_history_data();
    solver_control.log_history(true);
    solver_control.log_result(true);
#endif

    SolverFGMRES<TrilinosWrappers::BlockVector>
            solver (solver_control,
                    SolverFGMRES<TrilinosWrappers::BlockVector >::AdditionalData(100)); // With restart of 100

#ifdef LOGRUN
    deallog.push("Outer");
#endif
    solver.solve(ElasticBase<dim>::system_matrix, ElasticBase<dim>::solution, ElasticBase<dim>::system_rhs, preconditioner);
#ifdef LOGRUN
    deallog.pop();
#endif

    ElasticBase<dim>::par->system_iter = solver_control.last_step();
}


#endif // ELASTIC_2_BLOCK_H

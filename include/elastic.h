/* TODO

 */

#ifndef ELASTIC_H
#define ELASTIC_H

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
#include "boundary.h"
#include "coefficient.h"
#include "exact.h"
#include "parameters.h"
#include "preconditioner.h"
#include "rhs.h"
#include "SurfaceDataOut.h"

using namespace std;
using namespace dealii;
namespace Elastic
{
template <int dim>
class ElasticProblem : public ElasticBase<dim> {
public:
    ElasticProblem (const unsigned int degree, const int _info);

private:
    // Setup Algebraic multigrid(AMG)
    virtual void setup_AMG ();
    // Solve the system
    virtual void solve ();

    std_cxx1x::shared_ptr<typename Preconditioner::inner> A0_preconditioner;
    std_cxx1x::shared_ptr<typename Preconditioner::inner> A1_preconditioner;
    std_cxx1x::shared_ptr<typename Preconditioner::schur> S_preconditioner;
};
}

/*
 ------------- IMPLEMENTATION --------------
 */
template <int dim>
Elastic::ElasticProblem<dim>::ElasticProblem (const unsigned int degree, const int _info)
    : Elastic::ElasticBase<dim>(degree, _info, dim+1){}

/**
 * Setup preconditioners
 */
template <int dim>
void
Elastic::ElasticProblem<dim>::setup_AMG ()
{
    // Reset preconditioners and matrices
    A0_preconditioner.reset ();
    A1_preconditioner.reset ();
    S_preconditioner.reset ();

    A0_preconditioner
            = std_cxx1x::shared_ptr<typename Preconditioner::inner>(new typename Preconditioner::inner());

    A1_preconditioner
            = std_cxx1x::shared_ptr<typename Preconditioner::inner>(new typename Preconditioner::inner());
    
    S_preconditioner
            = std_cxx1x::shared_ptr<typename Preconditioner::schur>(new typename Preconditioner::schur());

    std::vector<std::vector<bool> > constant_modes;
    std::vector<bool>  displacement_components (dim+1,false);
    
    // A00
    displacement_components[0] = true;
    DoFTools::extract_constant_modes (ElasticBase<dim>::dof_handler,
                                      displacement_components,
                                      constant_modes);
    
    TrilinosWrappers::PreconditionAMG::AdditionalData amg_A0;
    amg_A0.constant_modes = constant_modes;
    
    amg_A0.elliptic = true;
    amg_A0.higher_order_elements = false;
    amg_A0.smoother_sweeps = 2;
    amg_A0.aggregation_threshold = ElasticBase<dim>::par->threshold;
    
    
    A0_preconditioner->initialize( ElasticBase<dim>::system_preconditioner.block(0,0),amg_A0);
    //
    
    //A11
    displacement_components[0] = false;
    displacement_components[1] = true;
    DoFTools::extract_constant_modes (ElasticBase<dim>::dof_handler,
                                      displacement_components,
                                      constant_modes);
    
    TrilinosWrappers::PreconditionAMG::AdditionalData amg_A1;
    amg_A1.constant_modes = constant_modes;
    
    amg_A1.elliptic = true;
    amg_A1.higher_order_elements = false;
    amg_A1.smoother_sweeps = 2;
    amg_A1.aggregation_threshold = ElasticBase<dim>::par->threshold;
    
    A1_preconditioner->initialize( ElasticBase<dim>::system_preconditioner.block(1,1),amg_A1);
    //
    
    TrilinosWrappers::PreconditionAMG::AdditionalData amg_S;
    
    amg_S.elliptic = true;
    amg_S.higher_order_elements = false;
    amg_S.smoother_sweeps = 2;
    amg_S.aggregation_threshold = ElasticBase<dim>::par->threshold;
    
    // elem-by-elem Schur
    S_preconditioner->initialize( ElasticBase<dim>::system_preconditioner.block(2,2), amg_S);
}

template <int dim>
void
Elastic::ElasticProblem<dim>::solve ()
{
    
    const BlockSchurPreconditioner<typename Preconditioner::inner, // A0, schur
            typename Preconditioner::schur>
            preconditioner( ElasticBase<dim>::system_preconditioner,
                            *A0_preconditioner,
                            *A1_preconditioner,
                            *S_preconditioner); // system_matrix
    
    SolverControl solver_control (ElasticBase<dim>::system_matrix.m(),
                                  ElasticBase<dim>::par->TOL*ElasticBase<dim>::system_rhs.l2_norm());
    
    SolverGMRES<TrilinosWrappers::BlockVector>
            solver (solver_control,
                    SolverGMRES<TrilinosWrappers::BlockVector >::AdditionalData(100));
    
    solver.solve(ElasticBase<dim>::system_matrix, ElasticBase<dim>::solution, ElasticBase<dim>::system_rhs, preconditioner);
    
    ElasticBase<dim>::par->system_iter = solver_control.last_step();
    deallog << "\t\tSchur " << ElasticBase<dim>::par->system_iter << std::endl;
}

#endif

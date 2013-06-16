/* TODO
 - Revise print_matlab method code
 - Revise Matlab_print_matrix code
*/

#include <deal.II/base/timer.h>
#include <deal.II/base/quadrature_lib.h>
#include <deal.II/base/logstream.h>
#include <deal.II/base/function.h>
#include <deal.II/base/utilities.h>
#include <deal.II/base/convergence_table.h>

#include <deal.II/lac/block_vector.h>
#include <deal.II/lac/full_matrix.h>
#include <deal.II/lac/block_sparse_matrix.h>

#include <deal.II/lac/solver_cg.h>
#include <deal.II/lac/solver_gmres.h>

#include <deal.II/lac/precondition.h>
#include <deal.II/lac/constraint_matrix.h>

#include <deal.II/grid/tria.h>
#include <deal.II/grid/grid_generator.h>
#include <deal.II/grid/tria_accessor.h>
#include <deal.II/grid/tria_iterator.h>
#include <deal.II/grid/tria_boundary_lib.h>
#include <deal.II/grid/grid_tools.h>

#include <deal.II/dofs/dof_handler.h>
#include <deal.II/dofs/dof_renumbering.h>
#include <deal.II/dofs/dof_accessor.h>
#include <deal.II/dofs/dof_tools.h>

#include <deal.II/fe/fe_q.h>
#include <deal.II/fe/fe_system.h>
#include <deal.II/fe/fe_values.h>

#include <deal.II/numerics/vector_tools.h>
#include <deal.II/numerics/matrix_tools.h>
#include <deal.II/numerics/data_out.h>
#include <deal.II/numerics/error_estimator.h>
//#include <deal.II/numerics/vector_tools.h> // replaced by: <deal.II/numerics/vectors.h>

#include <deal.II/lac/trilinos_sparse_matrix.h>
#include <deal.II/lac/trilinos_block_sparse_matrix.h>
#include <deal.II/lac/trilinos_vector.h>
#include <deal.II/lac/trilinos_block_vector.h>
#include <deal.II/lac/trilinos_precondition.h>

#include <deal.II/lac/sparse_direct.h>
#include <deal.II/lac/sparse_ilu.h>

#include <fstream>
#include <sstream>
#include <string>
#include <typeinfo>

#include "parameters.h"

#ifndef PRECONDITIONER_H
#define PRECONDITIONER_H

#define ZERO 1.0e-8

using namespace std;
using namespace dealii;
namespace Elastic
{
    // Preconditioner structure
    struct Preconditioner
    {
        //typedef SparseDirectUMFPACK		inner;
        //typedef SparseILU<double>		inner; // long time computing preconditioner
        typedef TrilinosWrappers::PreconditionAMG schur;
        typedef TrilinosWrappers::PreconditionAMG inner; // set AMG true
    };
    /*
     ************** PRECONDITIONER **************
     */
	// new code step-31
	template <class PreconditionerA, class PreconditionerS>
	class BlockSchurPreconditioner : public Subscriptor
    {
    public:
        BlockSchurPreconditioner (const TrilinosWrappers::BlockSparseMatrix     &S,
                                  const PreconditionerA           &Apreconditioner,
                                  const PreconditionerS           &Spreconditioner,
                                  parameters						*_par);
        
        void vmult (TrilinosWrappers::BlockVector       &dst,
                    const TrilinosWrappers::BlockVector &src) const;
        
    private:
    	// pointer to parameter object
    	parameters *par;
        const SmartPointer<const TrilinosWrappers::BlockSparseMatrix> s_matrix;
        const PreconditionerA &a_preconditioner;
        const PreconditionerS &s_preconditioner;
        
        mutable TrilinosWrappers::Vector tmp;
    };
}

/*
     ------------- IMPLEMENTATION --------------
*/
template <class PreconditionerA, class PreconditionerS>
Elastic::BlockSchurPreconditioner<PreconditionerA, PreconditionerS>::
BlockSchurPreconditioner(const TrilinosWrappers::BlockSparseMatrix  &S,
                         const PreconditionerA                      &Apreconditioner,
                         const PreconditionerS                      &Spreconditioner,
                         parameters									*_par)
    :
      par 					(_par),
      s_matrix				(&S),
      a_preconditioner        (Apreconditioner),
      s_preconditioner        (Spreconditioner),
      tmp                     (s_matrix->block(1,1).m())
{}

template <class PreconditionerA, class PreconditionerS>
void
Elastic::BlockSchurPreconditioner<PreconditionerA, PreconditionerS>::
vmult (TrilinosWrappers::BlockVector       &dst,
       const TrilinosWrappers::BlockVector &src) const
{

    // Solver control for solving the block with A^{-1}
    SolverControl control_inv (s_matrix->block(0,0).m(),
                               par->InvMatPreTOL*src.block(0).l2_norm());
    control_inv.enable_history_data ();
    control_inv.log_history (true);
    control_inv.log_result (true);

    SolverGMRES<TrilinosWrappers::Vector> // SchurTOL
            solver (control_inv, SolverGMRES<TrilinosWrappers::Vector >::AdditionalData(100));

    // Solve the block system for A^{-1}
    solver.solve(s_matrix->block(0,0), dst.block(0), src.block(0), a_preconditioner);

    // Push number of inner iterations to solve first block.
    par->inv_iterations.push_back(control_inv.last_step());

    // Write number of inner iterations to log file.
    deallog << "\t\tInner " << control_inv.last_step() << ", with TOL = "<< par->InvMatPreTOL*src.block(0).l2_norm() << std::endl;

    // a_preconditioner.vmult (dst.block(0), src.block(0));

    s_matrix->block(1,0).residual(tmp, dst.block(0), src.block(1));
    tmp *= -1;

    if(par->one_schur_it){// Use one iteration to find Schure complement
        s_preconditioner.vmult (dst.block(1), tmp);
    }
    else{
        SolverControl control_s (s_matrix->block(1,1).m(),
                                 par->SchurTOL*tmp.l2_norm());
        control_s.enable_history_data ();
        control_s.log_history (true);
        control_s.log_result (true);

        SolverGMRES<TrilinosWrappers::Vector> // SchurTOL
                solver (control_s, SolverGMRES<TrilinosWrappers::Vector >::AdditionalData(100));

        solver.solve(s_matrix->block(1,1), dst.block(1), tmp, s_preconditioner);

        // Push number of inner iterations for computing Schure complement.
        par->schur_iterations.push_back(control_s.last_step());

        // Write number of inner iterations for computing Schure complement to file.
        deallog << "\t\tSchur " << control_s.last_step() << ", with TOL = "<< par->SchurTOL*tmp.l2_norm() << std::endl;
    }
}
/*
     **************************************
*/

#endif

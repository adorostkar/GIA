/*! TODO
 */
#ifndef PRECONDITIONER_H
#define PRECONDITIONER_H

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

using namespace std;
using namespace dealii;
namespace Elastic
{
	// new code step-31
    template <class PreconditionerA, class PreconditionerS>
	class BlockSchurPreconditioner : public Subscriptor
    {
    public:
        BlockSchurPreconditioner (const TrilinosWrappers::BlockSparseMatrix     &S,
                                  const PreconditionerA            &A0preconditioner,
                                  const PreconditionerA            &A1preconditioner,
                                  const PreconditionerS            &Spreconditioner);
        
        void vmult (TrilinosWrappers::BlockVector       &dst,
                    const TrilinosWrappers::BlockVector &src) const;
        
    private:
    	// pointer to parameter object
        parameters *par;
        const SmartPointer<const TrilinosWrappers::BlockSparseMatrix> s_matrix;
        const PreconditionerA  &a0_preconditioner;
        const PreconditionerA  &a1_preconditioner;
        const PreconditionerS  &s_preconditioner;
        
        mutable TrilinosWrappers::Vector tmp;
    };
}

/*
     ------------- IMPLEMENTATION --------------
*/
template <class PreconditionerA, class PreconditionerS>
Elastic::BlockSchurPreconditioner<PreconditionerA, PreconditionerS>::
BlockSchurPreconditioner(const TrilinosWrappers::BlockSparseMatrix  &S,
                         const PreconditionerA                       &A0preconditioner,
                         const PreconditionerA                       &A1preconditioner,
                         const PreconditionerS                       &Spreconditioner)
    :
      s_matrix				(&S),
      a0_preconditioner        (A0preconditioner),
      a1_preconditioner        (A1preconditioner),
      s_preconditioner        (Spreconditioner),
      tmp                     (s_matrix->block(2,2).m())
{
    par = parameters::getInstance();
}

template <class PreconditionerA, class PreconditionerS>
void
Elastic::BlockSchurPreconditioner<PreconditionerA, PreconditionerS>::
vmult (TrilinosWrappers::BlockVector       &dst,
       const TrilinosWrappers::BlockVector &src) const
{

    // Solver for solving the block with A0^{-1}
    SolverControl control_inv0 (s_matrix->block(0,0).m(),
                               par->InvMatPreTOL*src.block(0).l2_norm());

#ifdef LOGRUN
    control_inv0.enable_history_data ();
    control_inv0.log_history (true);
    control_inv0.log_result (true);
#endif

    SolverFGMRES<TrilinosWrappers::Vector> // SchurTOL
            solver0 (control_inv0, SolverFGMRES<TrilinosWrappers::Vector >::AdditionalData(100));

#ifdef LOGRUN
    deallog.push("A1");
#endif
    solver0.solve(s_matrix->block(0,0), dst.block(0), src.block(0), a0_preconditioner);
#ifdef LOGRUN
    deallog.pop();
#endif

    // Solve the block system for A1^{-1}
    SolverControl control_inv1 (s_matrix->block(1,1).m(),
                               par->InvMatPreTOL*src.block(1).l2_norm());
#ifdef LOGRUN
    control_inv1.enable_history_data ();
    control_inv1.log_history (true);
    control_inv1.log_result (true);
#endif

    SolverFGMRES<TrilinosWrappers::Vector> // SchurTOL
            solver1 (control_inv1, SolverFGMRES<TrilinosWrappers::Vector >::AdditionalData(100));

#ifdef LOGRUN
    deallog.push("A2");
#endif
    solver1.solve(s_matrix->block(1,1), dst.block(1), src.block(1), a1_preconditioner);
#ifdef LOGRUN
    deallog.pop();
#endif

    // Push number of inner iterations to solve first block.
    par->inv_iterations.push_back((control_inv0.last_step() + control_inv1.last_step())/2);

    s_matrix->block(2,0).residual(tmp, dst.block(0),src.block(2));
    tmp *= -1;
    s_matrix->block(2,1).vmult_add(tmp, dst.block(1));


    SolverControl control_s (s_matrix->block(2,2).m(),
                             par->SchurTOL*tmp.l2_norm());
#ifdef LOGRUN
    control_s.enable_history_data ();
    control_s.log_history (true);
    control_s.log_result (true);
#endif

    SolverFGMRES<TrilinosWrappers::Vector> // SchurTOL
            solver (control_s, SolverFGMRES<TrilinosWrappers::Vector >::AdditionalData(100));

#ifdef LOGRUN
    deallog.push("Schur");
#endif
    solver.solve(s_matrix->block(2,2), dst.block(2), tmp, s_preconditioner);
#ifdef LOGRUN
    deallog.pop();
#endif

    // Push number of inner iterations for computing Schure complement.
    par->schur_iterations.push_back(control_s.last_step());
}


#endif

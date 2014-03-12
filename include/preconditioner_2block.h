/*! TODO
 */

#ifndef PRECONDITIONER_2BLOCK_H
#define PRECONDITIONER_2BLOCK_H

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
class Preconditioner2Blocks : public Subscriptor
{
public:
    Preconditioner2Blocks (const TrilinosWrappers::BlockSparseMatrix     &S,
                           const PreconditionerA           &Apreconditioner,
                           const PreconditionerS            &Spreconditioner);

    void vmult (TrilinosWrappers::BlockVector       &dst,
                const TrilinosWrappers::BlockVector &src) const;

private:
    // pointer to parameter object
    parameters *par;
    const SmartPointer<const TrilinosWrappers::BlockSparseMatrix> s_matrix;
    const PreconditionerA &a_preconditioner;
    const PreconditionerS  &s_preconditioner;

    mutable TrilinosWrappers::Vector tmp;
};
}

/*
     ------------- IMPLEMENTATION --------------
*/
template <class PreconditionerA, class PreconditionerS>
Elastic::Preconditioner2Blocks<PreconditionerA, PreconditionerS>::
Preconditioner2Blocks(const TrilinosWrappers::BlockSparseMatrix  &S,
                      const PreconditionerA                      &Apreconditioner,
                      const PreconditionerS                       &Spreconditioner)
    :
      s_matrix				(&S),
      a_preconditioner        (Apreconditioner),
      s_preconditioner        (Spreconditioner),
      tmp                     (s_matrix->block(1,1).m())
{
    par = parameters::getInstance();
}

template <class PreconditionerA, class PreconditionerS>
void
Elastic::Preconditioner2Blocks<PreconditionerA, PreconditionerS>::
vmult (TrilinosWrappers::BlockVector       &dst,
       const TrilinosWrappers::BlockVector &src) const
{
    double inv_relative_tol = par->InvMatPreTOL*src.block(0).l2_norm(),
            schur_relative_tol;
    // Solver for solving the block with A0^{-1}
    SolverControl control_inv (s_matrix->block(0,0).m(),
                               inv_relative_tol);
    control_inv.enable_history_data ();
    control_inv.log_history (true);
    control_inv.log_result (true);

    SolverGMRES<TrilinosWrappers::Vector> // SchurTOL
            solver_A (control_inv, SolverGMRES<TrilinosWrappers::Vector >::AdditionalData(100));

    deallog.push("A1");
    solver_A.solve(s_matrix->block(0,0), dst.block(0), src.block(0), a_preconditioner);
    deallog.pop();

    // Push number of inner iterations to solve first block.
    par->inv_iterations.push_back(control_inv.last_step());// + control_inv1.last_step());

    s_matrix->block(1,0).residual(tmp, dst.block(0),src.block(1));
    tmp *= -1;

    schur_relative_tol = par->SchurTOL*tmp.l2_norm();

    if(par->one_schur_it){// Use one iteration to find Schure complement
        s_preconditioner.vmult (dst.block(1), tmp);
        par->schur_iterations.push_back(1);
    }
    else{
        SolverControl control_s (s_matrix->block(1,1).m(),
                                 schur_relative_tol);
        control_s.enable_history_data ();
        control_s.log_history (true);
        control_s.log_result (true);

        SolverGMRES<TrilinosWrappers::Vector> // SchurTOL
                solver (control_s, SolverGMRES<TrilinosWrappers::Vector >::AdditionalData(100));

        deallog.push("Schur");
        solver.solve(s_matrix->block(1,1), dst.block(1), tmp, s_preconditioner);
        deallog.pop();

        // Push number of inner iterations for computing Schure complement.
        par->schur_iterations.push_back(control_s.last_step());
    }
}

#endif // PRECONDITIONER_2BLOCK_H

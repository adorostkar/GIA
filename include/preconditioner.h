/* TODO

*/

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
    template <class PreconditionerA0, class PreconditionerA1, class PreconditionerS>
	class BlockSchurPreconditioner : public Subscriptor
    {
    public:
        BlockSchurPreconditioner (const TrilinosWrappers::BlockSparseMatrix     &S,
                                  const PreconditionerA0           &A0preconditioner,
                                  const PreconditionerA1           &A1preconditioner,
                                  const PreconditionerS            &Spreconditioner,
                                  parameters                       *_par);
        
        void vmult (TrilinosWrappers::BlockVector       &dst,
                    const TrilinosWrappers::BlockVector &src) const;
        
    private:
    	// pointer to parameter object
    	parameters *par;
        const SmartPointer<const TrilinosWrappers::BlockSparseMatrix> s_matrix;
        const PreconditionerA0 &a0_preconditioner;
        const PreconditionerA1 &a1_preconditioner;
        const PreconditionerS  &s_preconditioner;
        
        mutable TrilinosWrappers::Vector tmp;
    };
}

/*
     ------------- IMPLEMENTATION --------------
*/
template <class PreconditionerA0, class PreconditionerA1, class PreconditionerS>
Elastic::BlockSchurPreconditioner<PreconditionerA0, PreconditionerA1, PreconditionerS>::
BlockSchurPreconditioner(const TrilinosWrappers::BlockSparseMatrix  &S,
                         const PreconditionerA0                      &A0preconditioner,
                         const PreconditionerA1                      &A1preconditioner,
                         const PreconditionerS                       &Spreconditioner,
                         parameters									 *_par)
    :
      par 					(_par),
      s_matrix				(&S),
      a0_preconditioner        (A0preconditioner),
      a1_preconditioner        (A1preconditioner),
      s_preconditioner        (Spreconditioner),
      tmp                     (s_matrix->block(2,2).m())
{}

template <class PreconditionerA0, class PreconditionerA1, class PreconditionerS>
void
Elastic::BlockSchurPreconditioner<PreconditionerA0, PreconditionerA1, PreconditionerS>::
vmult (TrilinosWrappers::BlockVector       &dst,
       const TrilinosWrappers::BlockVector &src) const
{

    // Solver for solving the block with A0^{-1}
    SolverControl control_inv0 (s_matrix->block(0,0).m(),
                               par->InvMatPreTOL*src.block(0).l2_norm());
    control_inv0.enable_history_data ();
    control_inv0.log_history (true);
    control_inv0.log_result (true);

    SolverGMRES<TrilinosWrappers::Vector> // SchurTOL
            solver0 (control_inv0, SolverGMRES<TrilinosWrappers::Vector >::AdditionalData(100));

    solver0.solve(s_matrix->block(0,0), dst.block(0), src.block(0), a0_preconditioner);

    // Solve the block system for A1^{-1}
    SolverControl control_inv1 (s_matrix->block(1,1).m(),
                               par->InvMatPreTOL*src.block(1).l2_norm());
    control_inv1.enable_history_data ();
    control_inv1.log_history (true);
    control_inv1.log_result (true);

    SolverGMRES<TrilinosWrappers::Vector> // SchurTOL
            solver1 (control_inv1, SolverGMRES<TrilinosWrappers::Vector >::AdditionalData(100));
    solver1.solve(s_matrix->block(1,1), dst.block(1), src.block(1), a1_preconditioner);

    // Push number of inner iterations to solve first block.
    par->inv_iterations.push_back(control_inv0.last_step());// + control_inv1.last_step());

    // Write number of inner iterations to log file.
    deallog << "\t\tInner First block" << control_inv0.last_step() << ", with TOL = "<< par->InvMatPreTOL*src.block(0).l2_norm() << std::endl;
    deallog << "\t\tInner second block" << control_inv1.last_step() << ", with TOL = "<< par->InvMatPreTOL*src.block(1).l2_norm() << std::endl;

    // a_preconditioner.vmult (dst.block(0), src.block(0));

//    s_matrix->block(1,0).residual(tmp, dst.block(0), src.block(1));
//    tmp *= -1;
    s_matrix->block(2,0).residual(tmp, dst.block(0),src.block(2));
    tmp *= -1;
    s_matrix->block(2,1).vmult_add(tmp, dst.block(1));


    if(par->one_schur_it){// Use one iteration to find Schure complement
        s_preconditioner.vmult (dst.block(2), tmp);
        par->schur_iterations.push_back(1);
    }
    else{
        SolverControl control_s (s_matrix->block(2,2).m(),
                                 par->SchurTOL*tmp.l2_norm());
        control_s.enable_history_data ();
        control_s.log_history (true);
        control_s.log_result (true);

        SolverGMRES<TrilinosWrappers::Vector> // SchurTOL
                solver (control_s, SolverGMRES<TrilinosWrappers::Vector >::AdditionalData(100));

        solver.solve(s_matrix->block(2,2), dst.block(2), tmp, s_preconditioner);

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

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

#ifndef EXACT_H
#define EXACT_H

#define ZERO 1.0e-8

using namespace std;
using namespace dealii;
namespace Elastic
{
    // Exact solution, measured in par->cases = 2 (Uniform load)
	template <int dim>
	class ExactSolution : public Function<dim>
	{
    public:
        ExactSolution (parameters *_par) : Function<dim>(dim+1), par(_par) {}
		
		virtual void vector_value (const Point<dim> &p,
								   Vector<double>   &value) const;
	private:
		// pointer to parameter object
		parameters *par;
	};
}
	
/*
     ------------- IMPLEMENTATION --------------
*/
template <int dim>
void
Elastic::ExactSolution<dim>::vector_value (const Point<dim> &p, Vector<double>   &values) const{
    Assert (values.size() == dim+1, ExcDimensionMismatch (values.size(), dim+1));

    // Variables fetched from "par" object for readability.
    double L = par->L, alpha = par->alpha, delta = par->delta, gamma = par->gamma;
    double rho_r = par->rho_r, g0 = par->g0;


    const double yb = par->y1*L; // scaled bottom
    const double x  = p[0]*L;
    const double y  = p[1]*L;
    const double A  = par->delta * g0 * par->rho_i * par->h;
    const double p0 = gamma * par->rho_i * g0 * par->h;

    if(par->cases() == 2){ // Uniform load
        values(0) = 0;
        if(par->adv_enabled){
            if(par->div_enabled){ // adv = 1, div = 1, complete
                values(1) = A*(yb-y);
                values(2) = -p0 * L;
            }
            else{				// adv = 1, div = 0, pre-stress of adv
                values(1) =  alpha * (exp(rho_r * g0 * delta * yb) - exp(rho_r * g0 * delta *y));
                values(2) = -p0 * exp(rho_r * g0 * delta * y) * L;
            }
        }else{					// adv = 0, div = 0, simple case
            values(1) = A*(yb-y);
            values(2) = -p0 * L;
        }
    }else{ // no exact solution available
        values(0) = 1e10;
        values(1) = 1e10;
        values(2) = 1e10;
    }
}
// end Exact solution
/*
     **************************************
*/

#endif

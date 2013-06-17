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

#ifndef COEFFICIENT_H
#define COEFFICIENT_H

#define ZERO 1.0e-8

using namespace std;
using namespace dealii;
namespace Elastic
{
	// ------------- coefficients mu, mu^2/lambda
	template <int dim>
	class Coefficients : public Function<dim>
	{
	public:
		double mu, beta;
		
        Coefficients (double E, double v, parameters *_par)  : Function<dim>(), mu( get_mu(E,v) ), beta( get_beta(E,v) ), par(_par) {}
		
		double get_mu(double E, double v);
		
		double get_beta(double E, double v);
		
		virtual double mu_value (const Point<dim>   &p, const unsigned int  component = 0) const;
		
		virtual void mu_value_list (const std::vector<Point<dim> > &points,
									std::vector<double>            &values,
									const unsigned int              component = 0) const;
		
		// beta = mu^2/alpha
		virtual double beta_value (const Point<dim>   &p, const unsigned int  component = 0) const;
		
		virtual void beta_value_list (const std::vector<Point<dim> > &points,
									  std::vector<double>            &values,
									  const unsigned int              component = 0) const;
	private:
		// pointer to parameter object
		parameters *par;
	};
}
/*
     ------------- IMPLEMENTATION --------------
*/
// coefficient transformation: get mu, mu^2/alpha with respect to E, v
template <int dim>
double
Elastic::Coefficients<dim>::get_mu(double E, double v)
{
    mu = E/( 2*(1+v) );
    return mu;
}

template <int dim>
double
Elastic::Coefficients<dim>::get_beta(double E, double v)
{
    beta = E*(1-2*v)/( 4*v*(1+v) );
    return beta;
}

template <int dim>
double
Elastic::Coefficients<dim>::mu_value (const Point<dim> &p,
                                      const unsigned int /*component*/) const
{
    /*// reserved for variable mu(x,y)
         if (p.square() < 0.5*0.5)
         return 20;
         else
         return 1;
         */
    return mu;
}
template <int dim>
void
Elastic::Coefficients<dim>::mu_value_list (const std::vector<Point<dim> > &points,
                                           std::vector<double>            &values,
                                           const unsigned int              component) const
{
    Assert (values.size() == points.size(), ExcDimensionMismatch (values.size(), points.size()));
    Assert (component == 0, ExcIndexRange (component, 0, 1));
    const unsigned int n_points = points.size();

    for (unsigned int i=0; i<n_points; ++i)
    {
        // // reserved for variable mu(x,y)
        // if (points[i].square() < 0.5*0.5)
        // values[i] = 20;
        // else
        // values[i] = 1;

        values[i] = mu;
    }
}

template <int dim>
double
Elastic::Coefficients<dim>::beta_value (const Point<dim> &p,
                                        const unsigned int /*component*/) const
{
    // // reserved for variable beta(x,y)
    // if (p.square() < 0.5*0.5)
    // 	return 20;
    // else
    // 	return 1;

    return beta;
}

template <int dim>
void
Elastic::Coefficients<dim>::beta_value_list (const 	std::vector<Point<dim> >		&points,
                                             std::vector<double>				&values,
                                             const unsigned int				component) const
{
    Assert (values.size() == points.size(), ExcDimensionMismatch (values.size(), points.size()));
    Assert (component == 0, ExcIndexRange (component, 0, 1));
    const unsigned int n_points = points.size();

    for (unsigned int i=0; i<n_points; ++i)
    {
        /* // reserved for variable beta(x,y)
             if (points[i].square() < 0.5*0.5)
             values[i] = 20;
             else
             values[i] = 1;
             */
        values[i] = beta;
    }
}
// end coefficients
/*
     **************************************
*/

#endif

/* TODO

*/

#include <deal.II/base/function.h>

#include <fstream>
#include <sstream>
#include <string>
#include <typeinfo>

#ifndef COEFFICIENT_H
#define COEFFICIENT_H

#include "parameters.h"

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
		
        Coefficients (double E, double v)  : Function<dim>(), mu( get_mu(E,v) ), beta( get_beta(E,v) ) {par = parameters::getInstance();}
		
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

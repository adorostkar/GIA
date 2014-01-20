/* TODO

*/

#include <deal.II/base/function.h>

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
        ExactSolution () : Function<dim>(dim+1) {par = parameters::getInstance();}
		
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

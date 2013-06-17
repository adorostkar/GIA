/* TODO

*/

#include <deal.II/base/function.h>

#include <fstream>
#include <sstream>
#include <string>
#include <typeinfo>

#include "parameters.h"

#ifndef BOUNDARY_H
#define BOUNDARY_H

#define ZERO 1.0e-8

using namespace std;
using namespace dealii;
namespace Elastic
{
	// Neumann Boundary conditions
	template <int dim>
    class BoundaryValues: public Function<dim>
	{
    public:
        BoundaryValues (parameters *_par) : Function<dim>(dim+1), par(_par) {}
		
		virtual double value (const Point<dim>   &p, const unsigned int  component = 0) const;
		
		virtual void vector_value (const Point<dim> &p, Vector<double>   &value) const;
	private:
		// pointer to parameter object
		parameters *par;
    };
}

/*
     ------------- IMPLEMENTATION --------------
*/

template <int dim>
double
Elastic::BoundaryValues<dim>::value (const Point<dim>  &p, const unsigned int component) const
{
    Assert (component < this->n_components, ExcIndexRange (component, 0, this->n_components));

    if ( (par->load_enabled) && (component == 1) ){
        if( (std::fabs(p[1] - par->y2) < ZERO) && ( p[0] <= par->Ix ) ){
            return par->load;
        }
    }
    return 0;
}

template <int dim>
void
Elastic::BoundaryValues<dim>::vector_value (const Point<dim> &p, Vector<double>   &values) const
{
    for (unsigned int c=0; c < this->n_components; ++c)
        values(c) = BoundaryValues<dim>::value (p, c);
}
// End Dirichlet Boundary conditions

/*
     **************************************
*/

#endif
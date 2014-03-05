/** TODO
*/

#ifndef RHS_H
#define RHS_H

#include <deal.II/base/function.h>

#include "parameters.h"

using namespace std;
using namespace dealii;
namespace Elastic
{
	// Right hand side // not working!
	template <int dim>
	class RightHandSide : public Function<dim>
	{
    public:
        RightHandSide () : Function<dim>(dim+1) {par = parameters::getInstance();}
		
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
Elastic::RightHandSide<dim>::value (const Point<dim>  &p, const unsigned int component) const
{
    if (component == 1) // y-component
        return (par->weight);

    return 0;
}

template <int dim>
void
Elastic::RightHandSide<dim>::vector_value (const Point<dim> &p, Vector<double>   &values) const
{
    for (unsigned int c=0; c<this->n_components; ++c)
        values(c) = RightHandSide<dim>::value (p, c);
}
// End right hand side
/*
     **************************************
*/

#endif

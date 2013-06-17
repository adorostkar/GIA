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

#ifndef RHS_H
#define RHS_H

#define ZERO 1.0e-8

using namespace std;
using namespace dealii;
namespace Elastic
{
	// Right hand side // not working!
	template <int dim>
	class RightHandSide : public Function<dim>
	{
    public:
        RightHandSide (parameters *_par) : Function<dim>(dim+1), par(_par) {}
		
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

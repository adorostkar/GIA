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

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

#ifndef ELASTIC_H
#define ELASTIC_H

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

	// Solver structure
	struct Solver
	{
		// typedef SolverCG<>		inner;
		// typedef SolverCG<>		schur;
		typedef SolverGMRES<TrilinosWrappers::Vector>	inner;
		typedef SolverGMRES<>	schur;
	};

	// Neumann Boundary conditions
	template <int dim>
    class BoundaryValues: public Function<dim>
	{
    public:
		BoundaryValues (parameters *_par) : par(_par), Function<dim>(dim+1) {}
		
		virtual double value (const Point<dim>   &p, const unsigned int  component = 0) const;
		
		virtual void vector_value (const Point<dim> &p, Vector<double>   &value) const;
	private:
		// pointer to parameter object
		parameters *par;
	};

	// Right hand side // not working!
	template <int dim>
	class RightHandSide : public Function<dim>
	{
    public:
		RightHandSide (parameters *_par) : par(_par), Function<dim>(dim+1) {}
		
		virtual double value (const Point<dim>   &p, const unsigned int  component = 0) const;
		
		virtual void vector_value (const Point<dim> &p, Vector<double>   &value) const;
	private:
		// pointer to parameter object
		parameters *par;
	};

	// ------------- coefficients mu, mu^2/lambda
	template <int dim>
	class Coefficients : public Function<dim>
	{
	public:
		double mu, beta;
		
		Coefficients (double E, double v, parameters *_par)  : par(_par), mu( get_mu(E,v) ), beta( get_beta(E,v) ), Function<dim>() {}
		
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

	// new code step-31
	template <class PreconditionerA, class PreconditionerS>
	class BlockSchurPreconditioner : public Subscriptor
    {
    public:
        BlockSchurPreconditioner (const TrilinosWrappers::BlockSparseMatrix     &S,
                                  const PreconditionerA           &Apreconditioner,
                                  const PreconditionerS           &Spreconditioner,
                                  parameters						*_par);
        
        void vmult (TrilinosWrappers::BlockVector       &dst,
                    const TrilinosWrappers::BlockVector &src) const;
        
    private:
    	// pointer to parameter object
    	parameters *par;
        const SmartPointer<const TrilinosWrappers::BlockSparseMatrix> s_matrix;
        const PreconditionerA &a_preconditioner;
        const PreconditionerS &s_preconditioner;
        
        mutable TrilinosWrappers::Vector tmp;
    };

    // Exact solution, measured in par->cases = 2 (Uniform load)
	template <int dim>
	class ExactSolution : public Function<dim>
	{
    public:
		ExactSolution (parameters *_par) : par(_par), Function<dim>(dim+1) {}
		
		virtual void vector_value (const Point<dim> &p,
								   Vector<double>   &value) const;
	private:
		// pointer to parameter object
		parameters *par;
	};
	
	template <int dim>
	class ElasticProblem
	{
    public:
        ElasticProblem (const unsigned int degree, parameters *_par);
        void run ();
        
    private:
        // conditional outputs
        ConditionalOStream info_0, info_1, info_2;
        // timer
        TimerOutput computing_timer;
    	// pointer to parameter object
    	parameters *par;
    	// Change string to uppercase
    	std::string to_upper(const std::string str);
    	// Create matlab file with the matrix
    	void matlab_print_matrix(const TrilinosWrappers::SparseMatrix &M, string filename );
    	// Create matlab file with the matrix
    	void matlab_print_matrix(const FullMatrix<double> &M, string filename );
    	// print matlab code for matrices
    	void print_matlab( string filename);
    	// Setup degree of freedom (DOF) of the system.
        void setup_dofs ();
        // Assemble the system
        void assemble_system ();
        // Setup Algebraic multigrid(AMG)
        void setup_AMG ();
        // Solve the system
        void solve ();
        // Computer the error
        void compute_errors () const;
        void output_results ();
        void surface_values ();
        
        const unsigned int						degree;
        std::vector<unsigned int>				boundary;
        
        Triangulation<dim>						triangulation;
        FESystem<dim>							fe;
        DoFHandler<dim>							dof_handler;
        
        BlockSparsityPattern					sparsity_pattern;
        TrilinosWrappers::BlockSparseMatrix		system_matrix;
        TrilinosWrappers::BlockSparseMatrix		system_preconditioner; 		// preconditioner [A 0;Bt S]
        
        TrilinosWrappers::BlockVector			solution;
        TrilinosWrappers::BlockVector			system_rhs, load, body_force, precond_rhs;
        
        std_cxx1x::shared_ptr<typename Preconditioner::inner> A_preconditioner;
        std_cxx1x::shared_ptr<typename Preconditioner::schur> S_preconditioner;
        //std_cxx1x::shared_ptr<TrilinosWrappers::PreconditionAMG> A_preconditioner;
	};
}

/* ==========================================================
   ======================= DECLARATION ======================
   ========================================================== */

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


template <class PreconditionerA, class PreconditionerS>
Elastic::BlockSchurPreconditioner<PreconditionerA, PreconditionerS>::
BlockSchurPreconditioner(const TrilinosWrappers::BlockSparseMatrix  &S,
						 const PreconditionerA                      &Apreconditioner,
						 const PreconditionerS                      &Spreconditioner,
						 parameters									*_par)
:
par 					(_par),
s_matrix				(&S),
a_preconditioner        (Apreconditioner),
s_preconditioner        (Spreconditioner),
tmp                     (s_matrix->block(1,1).m())
{}

template <class PreconditionerA, class PreconditionerS>
void
Elastic::BlockSchurPreconditioner<PreconditionerA, PreconditionerS>::
vmult (TrilinosWrappers::BlockVector       &dst,
	   const TrilinosWrappers::BlockVector &src) const
{
	
    // Solver control for solving the block with A^{-1}
	SolverControl control_inv (s_matrix->block(0,0).m(),
                               par->InvMatPreTOL*src.block(0).l2_norm());
	control_inv.enable_history_data ();
	control_inv.log_history (true);
	control_inv.log_result (true);
	
	SolverGMRES<TrilinosWrappers::Vector> // SchurTOL
            solver (control_inv, SolverGMRES<TrilinosWrappers::Vector >::AdditionalData(100));
	
    // Solve the block system for A^{-1}
	solver.solve(s_matrix->block(0,0), dst.block(0), src.block(0), a_preconditioner);
	
    // Push number of inner iterations to solve first block.
	par->inv_iterations.push_back(control_inv.last_step());

    // Write number of inner iterations to log file.
	deallog << "\t\tInner " << control_inv.last_step() << ", with TOL = "<< par->InvMatPreTOL*src.block(0).l2_norm() << std::endl;
	
	// a_preconditioner.vmult (dst.block(0), src.block(0));
	
	s_matrix->block(1,0).residual(tmp, dst.block(0), src.block(1));
	tmp *= -1;
	
    if(par->one_schur_it){// Use one iteration to find Schure complement
		s_preconditioner.vmult (dst.block(1), tmp);
	}
	else{
		SolverControl control_s (s_matrix->block(1,1).m(),
                                 par->SchurTOL*tmp.l2_norm());
		control_s.enable_history_data ();
		control_s.log_history (true);
		control_s.log_result (true);
		
		SolverGMRES<TrilinosWrappers::Vector> // SchurTOL
                solver (control_s, SolverGMRES<TrilinosWrappers::Vector >::AdditionalData(100));
		
		solver.solve(s_matrix->block(1,1), dst.block(1), tmp, s_preconditioner);
		
        // Push number of inner iterations for computing Schure complement.
		par->schur_iterations.push_back(control_s.last_step());

        // Write number of inner iterations for computing Schure complement to file.
		deallog << "\t\tSchur " << control_s.last_step() << ", with TOL = "<< par->SchurTOL*tmp.l2_norm() << std::endl;
	}
}

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


template <int dim>
Elastic::ElasticProblem<dim>::ElasticProblem (const unsigned int degree, parameters *_par)
:
  par(_par),
  info_0(std::cout, _par->info == 0),
  info_1(std::cout, _par->info == 1),
  info_2(std::cout, _par->info == 2),
  degree (degree),
  triangulation (Triangulation<dim>::maximum_smoothing),
  fe (FE_Q<dim>(degree+1), dim,
      FE_Q<dim>(degree), 1),
  dof_handler (triangulation),
  boundary (4, 0),
  computing_timer (info_0,
                   TimerOutput::summary,
                   TimerOutput::wall_times)
{
//    info_0.set_condition(par->info == 0);
//    info_1.set_condition(par->info == 1);
//    info_2.set_condition(par->info == 2);
}

template <int dim>
void
Elastic::ElasticProblem<dim>::setup_dofs ()
{
	
	//GridGenerator::hyper_cube (triangulation, par->left, par->right);
	
	// JC: labeling the faces of the gometry with boundary conditions
    par->load_enabled = true; // PROBLEM?? change in parameter is made.
	
	par->print_variables();
	
	std::vector<unsigned int> subdivisions (dim, 1);
	subdivisions[0] = par->xdivisions;
	subdivisions[1] = par->ydivisions;
    
	const Point<dim> bottom_left = (dim == 2 ?
                                        Point<dim>(par->x1,par->y1) :
                                        Point<dim>(par->x1,0,par->y1));

	const Point<dim> top_right   = (dim == 2 ?
                                        Point<dim>(par->x2,par->y2) :
                                        Point<dim>(par->x2,1,par->y2));
    
	GridGenerator::subdivided_hyper_rectangle (triangulation,
											   subdivisions,
											   bottom_left,
											   top_right);
	
	for (typename Triangulation<dim>::active_cell_iterator
         cell = triangulation.begin_active();
         cell != triangulation.end(); ++cell)
		for (unsigned int f=0; f < GeometryInfo<dim>::faces_per_cell; ++f){
			if (cell->face(f)->at_boundary()){
				
				const Point<dim> face_center = cell->face(f)->center();
				
                // If x component of the face's center is on the left boundary,
                // the face is one the left boundary
                if (face_center[dim-2] == par->x1)
					cell->face(f)->set_boundary_indicator(par->b_left);
                // If x component of the face's center is on the right boundary,
                // the face is on the right boundary.
                else if (face_center[dim-2] == par->x2)
					cell->face(f)->set_boundary_indicator(par->b_right);
                // If y component of the face's center is on the bottom boundary,
                // the face is on the bottom boundary.
                else if (face_center[dim-1] == par->y1)
					cell->face(f)->set_boundary_indicator(par->b_bottom);
                // If y component of the face's center is on the top boundary,
                // the face is on the top boundary.
                else if (face_center[dim-1] == par->y2)
					cell->face(f)->set_boundary_indicator(par->b_up);
			}// at boundary
		}// for faces
	
    // Refine the mesh with the number of refinement in the parameters.
	triangulation.refine_global (par->refinements);
	
	dof_handler.distribute_dofs (fe);
//    DoFRenumbering::Cuthill_McKee (dof_handler); // PROBLEM?? can we use it?
	
    std::vector<unsigned int> block_component (dim+1,0);
    block_component[dim] = 1;
    DoFRenumbering::component_wise (dof_handler, block_component);
    
	A_preconditioner.reset ();
	S_preconditioner.reset ();
	
	system_matrix.clear ();
	system_preconditioner.clear ();
	
	std::vector<unsigned int> dofs_per_block (2);
	DoFTools::count_dofs_per_block (dof_handler, dofs_per_block, block_component);
	const unsigned int n_u = dofs_per_block[0],
                       n_p = dofs_per_block[1];
	

    info_0 << "   Number of active cells: "
           << triangulation.n_active_cells()
           << ", Number of degrees of freedom: "
           << dof_handler.n_dofs()
           << " (" << n_u << '+' << n_p << ')'
           << std::endl;

	par->dofs << triangulation.n_active_cells() << "\t(" << n_u << '+' << n_p << ')';
	
	const unsigned int
            n_couplings = dof_handler.max_couplings_between_dofs();
	
	sparsity_pattern.reinit (2,2);
	sparsity_pattern.block(0,0).reinit (n_u, n_u, n_couplings);
	sparsity_pattern.block(1,0).reinit (n_p, n_u, n_couplings);
	sparsity_pattern.block(0,1).reinit (n_u, n_p, n_couplings);
	sparsity_pattern.block(1,1).reinit (n_p, n_p, n_couplings);
	sparsity_pattern.collect_sizes();
	
	DoFTools::make_sparsity_pattern (dof_handler, sparsity_pattern);
	sparsity_pattern.compress();
	
    // Print sparsity pattern of the matrix
    if(par->print_matrices){
        std::ofstream out ("sparsity_pattern.1");
        sparsity_pattern.print_gnuplot (out);
    }


	system_matrix.reinit (sparsity_pattern);
	system_preconditioner.reinit (sparsity_pattern);
	
	solution.reinit (2);
	solution.block(0).reinit (n_u);
	solution.block(1).reinit (n_p);
	solution.collect_sizes ();
	
	system_rhs.reinit (2);
	system_rhs.block(0).reinit (n_u);
	system_rhs.block(1).reinit (n_p);
	system_rhs.collect_sizes ();
	
	precond_rhs.reinit (2);
	precond_rhs.block(0).reinit (n_u);
	precond_rhs.block(1).reinit (n_p);
	precond_rhs.collect_sizes ();
	
	load.reinit (2);
	load.block(0).reinit (n_u);
	load.block(1).reinit (n_p);
	load.collect_sizes ();
	
	body_force.reinit (2);
	body_force.block(0).reinit (n_u);
	body_force.block(1).reinit (n_p);
	body_force.collect_sizes ();
}

template <int dim>
void
Elastic::ElasticProblem<dim>::assemble_system ()
{
	system_matrix=0;
	system_rhs=0;
	
	QGauss<dim>   quadrature_formula(degree+2);
	QGauss<dim-1> face_quadrature_formula(degree+2);
	
	FEValues<dim> fe_values (fe, quadrature_formula,
							 update_values    |
							 update_quadrature_points  |
							 update_JxW_values |
							 update_gradients);
	
	FEFaceValues<dim> fe_face_values (fe, face_quadrature_formula,
									  update_values    | update_normal_vectors |
									  update_quadrature_points  | update_JxW_values);
	
	const unsigned int   dofs_per_cell   = fe.dofs_per_cell;
	
	const unsigned int   n_q_points      = quadrature_formula.size();
	const unsigned int   n_face_q_points = face_quadrature_formula.size();
	
	// deal.ii local structure
	const unsigned int	 	dim_disp = 9,				// Q2 nodes
    dim_u = dim*dim_disp,		// dim * Q2 nodes
    dim_p = 4;					// 1 * Q1 nodes
	
	unsigned int 			l_u[] = {0,3,6,9,12,14,16,18,20, 1,4,7,10,13,15,17,19,21};
	unsigned int 			l_p[] = {2,5,8,11};
	
	unsigned int 			order[] = {0,3,6, 9,12,14,16,18,20, 1,4,7,10,13,15,17,19,21, 2,5,8,11};// use guido recommendation
	
	FullMatrix<double>	cell_matrix  (dofs_per_cell, dofs_per_cell),
    					cell_ordered (dofs_per_cell, dofs_per_cell),
    					cell_precond (dofs_per_cell, dofs_per_cell),
                        l_A          (dim_u,dim_u),
                        l_Bt         (dim_u,dim_p),
                        l_B          (dim_p,dim_u),
                        l_C          (dim_p,dim_p),
                        l_S          (dim_p,dim_p),
                        l_Ainv       (dim_u,dim_u),
                        l_Adiag      (dim_u,dim_u), // laplacian
                        aux          (dim_u,dim_u);
	
	Vector<double>       cell_rhs (dofs_per_cell);
	
	std::vector<unsigned int> local_dof_indices (dofs_per_cell);
	
	const RightHandSide<dim>			right_hand_side(par);
	const BoundaryValues<dim>			boundaries(par);
	std::vector<Vector<double> >		rhs_values (n_q_points, Vector<double>(dim+1));
	std::vector<Vector<double> >		boundary_values (n_face_q_points, Vector<double>(dim+1));
	
	Coefficients<dim> 				  	 coeff(par->YOUNG,par->POISSON, par);
	std::vector<double>     		  	 mu_values (n_q_points);
	std::vector<double>     		  	 beta_values (n_q_points);	// mu^2/alpha
	
	const FEValuesExtractors::Vector displacements (0);
	const FEValuesExtractors::Scalar pressure (dim);
	
	Tensor<1,dim> e(true);
    e[0] = 0;
    e[1] = 1.0; // 3rd dim??
	
	std::vector<SymmetricTensor<2,dim> > symgrad_phi_u	(dofs_per_cell);
	std::vector<Tensor<2,dim> >          grad_phi		(dofs_per_cell);
	std::vector<Tensor<1,dim> >			 phi_u			(dofs_per_cell);
	std::vector<double>                  div_phi_u		(dofs_per_cell);
	std::vector<double>                  phi_p			(dofs_per_cell);
	
	bool first = true;
	unsigned int counter = 0;
	double h, aux_p;
	
	typename DoFHandler<dim>::active_cell_iterator
	cell = dof_handler.begin_active(),
	endc = dof_handler.end();
	for (; cell!=endc; ++cell)
	{
		fe_values.reinit (cell);
		cell_matrix		= 0;
		cell_rhs		= 0;
		cell_precond	= 0;
		l_Adiag			= 0;
		
		right_hand_side.vector_value_list(fe_values.get_quadrature_points(),
										  rhs_values);
		
		boundaries.vector_value_list(fe_face_values.get_quadrature_points(),
                                     boundary_values);
		
		coeff.mu_value_list     (fe_values.get_quadrature_points(), mu_values);
		coeff.beta_value_list   (fe_values.get_quadrature_points(), beta_values);
		
		for (unsigned int q=0; q<n_q_points; ++q)
		{
			for (unsigned int k=0; k < dofs_per_cell; ++k)
			{
				symgrad_phi_u[k] = fe_values[displacements].symmetric_gradient (k, q);
				grad_phi[k]		 = fe_values[displacements].gradient (k, q);
				phi_u[k]		 = fe_values[displacements].value (k, q);
				div_phi_u[k]     = fe_values[displacements].divergence (k, q);
				phi_p[k]         = fe_values[pressure].value (k, q);
			}
			
			for (unsigned int i=0; i<dofs_per_cell; ++i)
			{
				const unsigned int component_i =
				fe.system_to_component_index(i).first;
				
				for (unsigned int j=0; j < dofs_per_cell; ++j)
				{
					const unsigned int component_j =
					fe.system_to_component_index(j).first;
					
					cell_matrix(i,j) += (symgrad_phi_u[i] * symgrad_phi_u[j] * 2 * mu_values[q] // A
                                         - grad_phi[j]  * e * phi_u[i] * par->scale3 * par->adv	// A-adv
                                         + div_phi_u[j] * e * phi_u[i] * par->scale3 * par->div	// A-div
                                         + div_phi_u[i] * phi_p[j] * mu_values[q]				// Bt
                                         + phi_p[i] * div_phi_u[j] * mu_values[q]				// B
                                         - phi_p[i] * phi_p[j] * beta_values[q] )				// C
                                         * fe_values.JxW(q);
					
					cell_precond(i,j) += (
                                          phi_p[i] * div_phi_u[j] * mu_values[q]				// B
                                          )*
                                          fe_values.JxW(q);
					
				}// end j
				
                cell_rhs(i) +=  phi_u[i] * e * par->weight * fe_values.JxW(q); // load vector, body force
			}// end i
		} // end q
		
		// Neumann Boundary conditions (Ice-Load and free surface)
        for (unsigned int face_num=0; face_num<GeometryInfo<dim>::faces_per_cell; ++face_num){
            if (cell->face(face_num)->at_boundary() /* && (cell->face(face_no)->boundary_indicator() == par->b_up )*/ ){
                //if (cell->at_boundary(face_no))
                Point<dim> face_center = cell->face(face_num)->center();
				if(face_center[dim-1] == par->y2){ // y2 = surface or TOP boundary
					
                    fe_face_values.reinit (cell, face_num);
                    
					for (unsigned int q=0; q<n_face_q_points; ++q)
						for (unsigned int i=0; i<dofs_per_cell; ++i){
                            const unsigned int
                                    component_i = fe.system_to_component_index(i).first;
                            
							cell_rhs(i) +=  fe_face_values.shape_value(i, q) *
                                            boundary_values[q](component_i) *
                                            fe_face_values.JxW(q);
						}
				}// end if
			}// end if at boundary
		}// end face
		
		// Local assemble and Schur generation
		// extract here using velocities, pressure
        for (unsigned int i=0; i<dofs_per_cell; ++i){// shape index i
            for (unsigned int j=0; j < dofs_per_cell; ++j){// shape index j
				
				cell_ordered(i,j) = cell_matrix(order[i],order[j]); // local matrix ordered by u0,u1...un,v0,v1...vn,p0,p1...pm
				
				if( i < dim_u ){
                    if(j < dim_u){// i < dim_u, j <  dim_u
                        l_A(i,j)  = cell_ordered(i,j);
						aux(i,j)  = l_A(i,j);
						
						// bgn Diagonal 11 block
						if(par->precond == 0){
                            if(i < dim_disp && j < dim_disp){
                                l_Adiag(i,j) = cell_ordered(i,j); // first block
							}else if(j >= dim_disp){
								l_Adiag(i,j) = cell_ordered(i,j); // second block
							}
						}
						// end Diagonal 11 block
						
                    }else{// i < dim_u, j >= dim_u
                        l_Bt(i,j-dim_u) = cell_ordered(i,j);
					}
				}else{
                    if(j < dim_u){// i >= dim_u, j <  dim_u
                        l_B(i-dim_u,j)  = cell_ordered(i,j);
                    }else{// i >= dim_u, j >= dim_u
                        l_C(i-dim_u,j-dim_u)  = cell_ordered(i,j);
						l_S(i-dim_u,j-dim_u)  = l_C(i-dim_u,j-dim_u); // -C ... look at the sign
					}
				}
			}
		}
		// boundary conditions
		
		// Local Schur calculation	l_A(k,k) += h*h;// A(i,j) + h²I
		h = cell->diameter();
		aux.diagadd(h*h);
		l_Ainv.invert(aux);
		
        l_S.triple_product 	(	l_Ainv,l_B,l_Bt,false,false, -1.0  );
		// End Schur calculation
		
        // begin Schur assembly preconditioner
        for (unsigned int i=0; i< dim_p; ++i){// shape index i
			for (unsigned int j=0; j < dim_p; ++j){
				cell_precond(l_p[i],l_p[j]) = l_S(i,j);
			}
		}
		// end Schur assembly preconditioner
		
        // begin assembly A preconditioner
        if(par->precond == 0){
            for (unsigned int i=0; i< dim_u; ++i){// shape index i
                for (unsigned int j=0; j < dim_u; ++j){
                    cell_precond(l_u[i],l_u[j]) = l_Adiag(i,j);
                }
            }
        }else if(par->precond == 1){
            for (unsigned int i=0; i< dim_u; ++i){// shape index i
                for (unsigned int j=0; j < dim_u; ++j){
                    cell_precond(l_u[i],l_u[j]) = l_A(i,j);
                }
            }
        }
//      for (unsigned int i=0; i< dim_u; ++i){// shape index i
//			for (unsigned int j=0; j < dim_u; ++j){
//				if(par->precond == 0){
//					cell_precond(l_u[i],l_u[j]) = l_Adiag(i,j);
//				}else if(par->precond == 1){
//					cell_precond(l_u[i],l_u[j]) = l_A(i,j);
//				}
//			}
//		}
		// end assembly A preconditioner
		
		// printing local matrices
		if(par->print_local && first && (par->POISSON == 0.2) ){ // print_local, first
			string name = "l_m";
			//name += Utilities::int_to_string (counter, 2);
			matlab_print_matrix(cell_matrix,name);
			
			name = "l_p";
			//name += Utilities::int_to_string (counter, 2);
			matlab_print_matrix(cell_precond,name);
		}
		
		// local-to-global
		cell->get_dof_indices (local_dof_indices);
		for (unsigned int i=0; i<dofs_per_cell; ++i)
		{
			for (unsigned int j=0; j<dofs_per_cell; ++j){
				
				system_matrix.add (local_dof_indices[i],
                                   local_dof_indices[j],
                                   cell_matrix(i,j));
				
				system_preconditioner.add (local_dof_indices[i],
                                           local_dof_indices[j],
                                           cell_precond(i,j));
			}
			system_rhs(local_dof_indices[i]) += cell_rhs(i);
		}
		
		// end local-to-global
		
		first = false;
		counter++;
	} // end cell
}

// PROBLEM?? check the paper.
template <int dim>
void
Elastic::ElasticProblem<dim>::setup_AMG ()
{
	A_preconditioner
    = std_cxx1x::shared_ptr<typename Preconditioner::inner>(new typename Preconditioner::inner());
	
	S_preconditioner
    = std_cxx1x::shared_ptr<typename Preconditioner::schur>(new typename Preconditioner::schur());
	
	std::vector<std::vector<bool> > constant_modes;
	std::vector<bool>  displacement_components (dim,true);
	//displacement_components[dim] = false;
    DoFTools::extract_constant_modes (dof_handler,
                                      displacement_components,
									  constant_modes);
	
	TrilinosWrappers::PreconditionAMG::AdditionalData amg_A;
	amg_A.constant_modes = constant_modes;
	
	amg_A.elliptic = true;
	amg_A.higher_order_elements = false;
	amg_A.smoother_sweeps = 2;
	amg_A.aggregation_threshold = par->threshold;
	
    // A
    A_preconditioner->initialize( system_preconditioner.block(0,0),amg_A);
	
	TrilinosWrappers::PreconditionAMG::AdditionalData amg_S;
	
	amg_S.elliptic = true;
	amg_S.higher_order_elements = false;
	amg_S.smoother_sweeps = 2;
	amg_S.aggregation_threshold = par->threshold;
	
    // elem-by-elem Schur
    S_preconditioner->initialize( system_preconditioner.block(1,1), amg_S);
}

template <int dim>
void
Elastic::ElasticProblem<dim>::solve ()
{
	
	const BlockSchurPreconditioner<typename Preconditioner::inner, // A, schur
									typename Preconditioner::schur>
    preconditioner( system_preconditioner, *A_preconditioner, *S_preconditioner, par); // system_matrix
    
	SolverControl solver_control (system_matrix.m(),
								  par->TOL*system_rhs.l2_norm());
    
	SolverGMRES<TrilinosWrappers::BlockVector>
	solver (solver_control,
            SolverGMRES<TrilinosWrappers::BlockVector >::AdditionalData(100));
	
	solver.solve(system_matrix, solution, system_rhs, preconditioner);
	
	par->system_iter = solver_control.last_step();
	deallog << "\t\tSchur " << par->system_iter << std::endl;
}


template <int dim>
void
Elastic::ElasticProblem<dim>::compute_errors () const
{
	const ComponentSelectFunction<dim>
    pressure_mask (dim, dim+1);
	const ComponentSelectFunction<dim>
    velocity_mask(std::make_pair(0, dim), dim+1);
	
	ExactSolution<dim> exact_solution(par);
	Vector<double> cellwise_errors (triangulation.n_active_cells());
	
    
	QIterated<dim> quadrature ( QTrapez<1>(), degree+1);
	
	// L2-norm
	VectorTools::integrate_difference (dof_handler, solution, exact_solution,
									   cellwise_errors, quadrature,
									   VectorTools::L2_norm,
									   &pressure_mask);
	
	const double p_l2_error = cellwise_errors.l2_norm();
	
	VectorTools::integrate_difference (dof_handler, solution, exact_solution,
									   cellwise_errors, quadrature,
									   VectorTools::L2_norm,
									   &velocity_mask);
	
	const double u_l2_error = cellwise_errors.l2_norm();
	// end L2-norm
	

    info_0 << "Errors: ||e_u||_L2, ||e_p||_L2 = " << u_l2_error
           << "," << p_l2_error
           << std::endl;
    info_1 << u_l2_error << "\t" << p_l2_error <<"\t";
    info_2 << u_l2_error << "\t" << p_l2_error <<"\t";
}

template <int dim>
void
Elastic::ElasticProblem<dim>::output_results ()
{
	std::vector<std::string> solution_names (dim, "displacements");
	solution_names.push_back ("pressure");
	
	std::vector<DataComponentInterpretation::DataComponentInterpretation>
	data_component_interpretation
	(dim, DataComponentInterpretation::component_is_part_of_vector);
	data_component_interpretation
	.push_back (DataComponentInterpretation::component_is_scalar);
	
	DataOut<dim> data_out;
	data_out.attach_dof_handler (dof_handler);
	data_out.add_data_vector (solution, solution_names,
							  DataOut<dim>::type_dof_data,
							  data_component_interpretation);
	data_out.build_patches ();
	
	std::ostringstream filename;
	filename << "solution_" << par->POISSON	<< ".vtk";
	
	std::ofstream output (filename.str().c_str());
	data_out.write_vtk (output);
}

template <int dim>
void
Elastic::ElasticProblem<dim>::run ()
{
	double t_ass=0, t_solve=0, t_tot=0;

	setup_dofs ();

    info_0 << "   Assembling ... " << std::endl;

    computing_timer.enter_section("Assembling");
    assemble_system ();
    computing_timer.exit_section("Assembling");
	
    info_0 << "   AMG preconditioners ... " << std::endl;

    computing_timer.enter_section("AMG preconditioners");
	setup_AMG ();
    computing_timer.exit_section("AMG preconditioners");
	
	// applying the Dirichlet BC
	
	std::vector<bool> ns_mask (dim+1, true); // NO_SLIP
	
	ns_mask[0] = true;
	ns_mask[1] = true;
	ns_mask[2] = false;
	
	std::map<unsigned int,double> boundary_values;
	VectorTools::interpolate_boundary_values (dof_handler,
											  NO_SLIP,
											  ZeroFunction<dim>(),
											  boundary_values,
											  ns_mask);
	
	std::vector<bool> vs_mask (dim+1, true); // V_SLIP
	
	vs_mask[0] = true;
	vs_mask[1] = false;
	vs_mask[2] = false;
	
	VectorTools::interpolate_boundary_values (dof_handler,
											  V_SLIP,
											  ZeroFunction<dim>(),
											  boundary_values,
											  vs_mask);
	
	MatrixTools::apply_boundary_values (boundary_values,
										system_preconditioner,
										solution,
										precond_rhs,
										false);
	
	MatrixTools::apply_boundary_values (boundary_values,
										system_matrix,
										solution,
										system_rhs,
										false);
	
	if(par->print_matrices ){ // define print_data
		// printing matrices in matlab form, after Dirichlet B.C
		std::cout << "   ** Printing matrices in Matlab form **" << std::endl << std::flush;
		
		matlab_print_matrix(system_matrix.block(0,0),"a00");
		matlab_print_matrix(system_matrix.block(0,1),"a01");
		matlab_print_matrix(system_matrix.block(1,0),"a10");
		matlab_print_matrix(system_matrix.block(1,1),"a11");
		
		matlab_print_matrix(system_preconditioner.block(0,0),"p00");
		matlab_print_matrix(system_preconditioner.block(0,1),"p01");
		matlab_print_matrix(system_preconditioner.block(1,0),"p10");
		matlab_print_matrix(system_preconditioner.block(1,1),"p11");
	}
	
	int inv_iter = 0, schur_iter = 0;
	if(par->solve){
        info_0 << "   system solver ... " << std::endl;

        computing_timer.enter_section("System solver");
        solve ();
        computing_timer.exit_section("System solver");
		
		// printing solver info
		for(int i = 0; i < par->inv_iterations.size(); i++){
			inv_iter   += par->inv_iterations[i];
			schur_iter += par->schur_iterations[i];
		}
		
		// average inner iterations
		inv_iter = (int)(1.0*inv_iter/par->inv_iterations.size());
		schur_iter = (int)(1.0*schur_iter/par->schur_iterations.size());
		
		output_results ();
		surface_values ();
		
		if(par->cases() == 2 )
			compute_errors ();
	}
	
	// printing: matrices, par->info
	if(par->POISSON == 0.2 ){
		print_matlab("gia");
	}
	if(par->solve){
		int tempSpace = 15;
        info_0   << "GMRES iterations: system(<inv>,<schur>) = "
                 << par->system_iter
				 << "(" << inv_iter << ", " << schur_iter << ")"
				 << endl;

        info_1   << left << setw(tempSpace) << setfill(' ') << par->system_iter
				 << left << setw(tempSpace) << setfill(' ') << schur_iter
				 << left << setw(tempSpace) << setfill(' ') << t_ass
				 << left << setw(tempSpace) << setfill(' ') << t_solve
				 << left << setw(tempSpace) << setfill(' ') << par->dofs.str()
				 << endl;

        info_2   << left << setw(tempSpace) << setfill(' ') << par->system_iter
				 << left << setw(tempSpace) << setfill(' ') << inv_iter
				 << left << setw(tempSpace) << setfill(' ') << schur_iter
				 << left << setw(tempSpace) << setfill(' ') << t_ass
				 << left << setw(tempSpace) << setfill(' ') << t_solve
				 << left << setw(tempSpace) << setfill(' ') << par->dofs.str()
				 << endl;
	}
}

template <int dim>
void
Elastic::ElasticProblem<dim>::surface_values () {
	
	const unsigned int n = par->surf_samples;
	double dx;
	
	ostringstream filename;
	filename << "surface_values" << par->str_poisson << ".m";
	std::ofstream detector_data(filename.str().c_str());
	std::vector<Point<dim> > detector_locations;
	
	Vector<double> value;
	
	switch (par->cases()){
		case 2: // uniform load
			dx = (par->y2 - par->y1)/n;
			for (unsigned int i = 0; i < n; ++i){
				detector_locations.push_back (Point<dim> ( (par->x1 + par->x2 )*0.5,
                                                          par->y2 - dx*i ));
			}
			break;
		default:
			dx = (par->x2 - par->x1)/n;
			for (unsigned int i = 0; i < n; ++i){
				detector_locations.push_back (Point<dim> ( par->x1 + dx*i,
														  0.0 ));
			}
			break;
	}
	
	detector_data << "sol"<< par->str_poisson <<" = [";
	for (unsigned int i=0 ; i<detector_locations.size(); ++i){
		
		VectorTools::point_value (dof_handler,
								  solution,
								  detector_locations[i],
								  value);
		detector_data << (detector_locations[i](0)*par->L/1e3) << ", " << (detector_locations[i](1)*par->L/1e3) << ", " << value(0) << ", " << value(1) << ", " << (value(2)/par->L) <<";\n";
	}
	detector_data << "];" << std::endl;
	
    info_0 << "... extracting surface values ..."
           << " dx = " << dx << ", n = " << n
           << std::endl;
}

template <int dim>
void
Elastic::ElasticProblem<dim>::print_matlab( string filename){
	//Printing in matlab form
	string extension = filename + ".m";
	ofstream myfile(extension.c_str());
	
	int figure = 1;
	
	if(!myfile.is_open()){
		cout << "Print_matlab: Unable to open file...";
		return;
	}

	myfile << "close all;clear;" << endl;
	
	if(par->print_matrices){
		myfile << "data_a00" << endl;
		myfile << "data_a01" << endl;
		myfile << "data_a10" << endl;
		myfile << "data_a11" << endl;
		myfile << "" << endl;
		
		myfile << "A = [A00 A01; A10 A11];" << endl;
		myfile << "" << endl;
		
		myfile << "data_p00" << endl;
		myfile << "data_p01" << endl;
		myfile << "data_p10" << endl;
		myfile << "data_p11" << endl;
		myfile << "" << endl;
		
		myfile << "P = [P00 P01; P10 P11];" << endl;
		myfile << "" << endl;
		
		myfile << "Ainv = inv(A00);" << endl;
		myfile << "S = A11 - A10*Ainv*A01;" << endl;
		
		myfile << "figure(" << figure << ");mesh(S);title('Full Schur');" << endl;
		figure++;
		
		myfile << "figure(" << figure << ");mesh(P11);title('Element-by-element Schur');" << endl;
		figure++;
		
		myfile << "data_l_m;" << endl;
		myfile << "data_l_p;" << endl;
		//myfile << "data_l_o;" << endl;
		//myfile << "data_l_po;" << endl;
		myfile << "" << endl;
		
		myfile	<< "eigS=eig(S);\n"
				<< "eigSe=eig(P11);\n"
				<< "eigM=eig(A11);\n\n"
				<< "r_part=[real(eigS) real(eigSe) real(eigM)];\n"
				<< "i_part=[imag(eigS) imag(eigSe) imag(eigM)];\n"
				<< "figure("<< figure <<");plot(r_part)\n";
		figure++;
		
		myfile << "figure(" << figure << ");plot(r_part);title('Quality comparison of $S,\hat S$ and $C$','Interpreter','latex');"<<endl;
		figure++;

		myfile	<< "legend('$S$','$\hat S$','$C$','Interpreter','latex');\n"
				<< "h1 = legend;\n"
				<< "set(h1, 'interpreter', 'latex')\n"
				<< "figure("<< figure <<");plot(i_part);title('Quality comparison of $S,\hat S$ and $C$','Interpreter','latex');"<<endl;
		figure++;
        
	}
	
	std::stringstream xplot;
	xplot << "set(gca,'XTick',";
	
	if(par->cases() == 2)
		xplot << "-4000:500:0)\nxlabel('y(xt) ";
	else if(par->cases() == 1)
		xplot << "0:5000:"<<par->x2*1e4<<")\nxlabel('x(y=0) ";
	else if(par->cases() == 0)
		xplot << "0:2000:"<<par->x2*1e4<<")\nxlabel('x(y=0) ";
	
	xplot << "in (Km)');\n";
	
	myfile	<< "surface_values0_2;\n"
			<< "surface_values0_3;\n"
			<< "surface_values0_4;\n"
			<< "surface_values0_5;\n\n"
			<< "x = sol0_2(:,1);\n"
			<< "y = sol0_2(:,2);\n"
			<< "domain = " 
			<< ((par->cases() == 2)?"y(:,1)":"x(:,1)") << ";\n\n";

	myfile	<< "horizontal = [sol0_2(:,3) sol0_3(:,3) sol0_4(:,3) sol0_5(:,3)];\n"
			<< "vertical   = [sol0_2(:,4) sol0_3(:,4) sol0_4(:,4) sol0_5(:,4)];\n"
			<< "pressure   = [sol0_2(:,5) sol0_3(:,5) sol0_4(:,5) sol0_5(:,5)];\n\n";

	myfile	<< "% Change default axes fonts.\n"
			<< "set(0,'DefaultAxesFontName', 'Times New Roman');\n"
			<< "set(0,'DefaultAxesFontSize', 16);\n\n";

	myfile	<< "% Change default text fonts.\n"
			<< "set(0,'DefaultTextFontname', 'Times New Roman');\n"
			<< "set(0,'DefaultTextFontSize', 16);\n\n";

	myfile	<< "figure("<<figure<<"); plot(domain,horizontal(:,1:4));%title('horizontal');\n"
			<< xplot.str() 
			<< "ylabel('Horizontal displacement in (m)');\n"
			<< "hleg1 = legend('\\nu=0.2','\\nu=0.3','\\nu=0.4','\\nu=0.5', 'Location','SouthEast');\n"
			<< "grid on\n"
			<< "%set(gca, 'GridLineStyle', '-');\n" << endl;
	figure++;
	
	myfile	<< "figure("<< figure <<"); plot(domain,  vertical(:,1:4));%title('vertical');\n"
			<< xplot.str()
			<< "ylabel('Vertical displacement in (m)');\n"
			<< "hleg1 = legend('\\nu=0.2','\\nu=0.3','\\nu=0.4','\\nu=0.5', 'Location','SouthEast');\n"
			<< "grid on\n"
			<< "%set(gca, 'GridLineStyle', '-');\n\n" << endl;
	figure++;
	
	myfile	<< "% analytical constants, alpha = rho_i*h/rho_r, delta = (1+v)(1-2v)/(E(1-v)), 1/(2mu + lambda) = (1+v)(1-2v)/(E(1-v))\n"
			<< "v = 0.2;E = " << (par->S*par->YOUNG) << "; g = " << par->g0 << "; h = "<< par->h <<";\n"
			<< "rho_i = "<< par->rho_i << "; rho_r = " << par->rho_r << ";\n"
			<< "alpha = "<< par->alpha <<";\n"
			<< "beta  = "<< par->beta <<";\n"
			<< "delta = (1+v)*(1-2*v)/(E*(1-v));\n"
			<< "gamma = "<< par->gamma <<";\n";

	myfile	<< "yb    = " << (par->L*par->y1) << ";\n\n"
			<< "A1 = (1+v)*(1-2*v)/(E*(1-v))*g*rho_i*h; % rho_i\n\n"
			<< "mu     = E/(2*(1+v));\n"
			<< "lambda = E*v/((1+v)*(1-2*v));\n";
	
	myfile	<< "y = 1000*domain;%yb:abs(yb)/1000:0.0;\n"
			<< "v1 = A1.*(yb-y);\n"
			<< "v2 = alpha.*(exp(rho_r*g*delta.*yb)-exp(rho_r*g*delta.*y));\n"
			<< "p = -rho_i*g*h*gamma;\n\n";

	myfile	<< "%figure("<<figure<<");plot(y,v1,y,v2);\n" << endl;
	figure++;

	myfile	<< "p1 = sol0_2(:,4)\n;"
			<< "%figure("<<figure<<");plot(domain,p1,domain,p);\n"
			<< endl;
	
	myfile.close();
}

// Create a matlab file with matrices
template <int dim>
void
Elastic::ElasticProblem<dim>::matlab_print_matrix(const TrilinosWrappers::SparseMatrix &M, string filename ){
	//Printing in matlab form
	string extension = "data_" + filename + ".m";
	string up_filename = to_upper(filename);
	
	ofstream myfile(extension.c_str());
	
	if (!myfile.is_open()){
		cout << "Unable to open file";
		return;
	}
	myfile << up_filename.c_str() << " = [" << endl;
	
	for(unsigned int i = 0;i < M.m();i++){
		for(unsigned int j = 0;j < M.n();j++){
			
			myfile << M.el(i,j) << "\t";
		}
		myfile << ";" << endl;
	}
	myfile << "];" << endl;
	myfile << up_filename.c_str() << " = sparse(" << up_filename.c_str() << ");" << endl;
	myfile.close();
}

template <int dim>
void
Elastic::ElasticProblem<dim>::matlab_print_matrix(const FullMatrix<double> &M, string filename ){
	//Printing in matlab form
	string extension = "data_" + filename + ".m";
	to_upper(filename);
	
	ofstream myfile(extension.c_str());
	
	if (!myfile.is_open()){
		cout << "Unable to open file";
		return;
	}
	myfile << filename.c_str() << " = [" << endl;
	
	for(unsigned int i = 0;i < M.m();i++){
		for(unsigned int j = 0;j < M.n();j++){
			
			myfile << M(i,j) << "\t";
		}
		myfile << ";" << endl;
	}
	myfile << "];" << endl;
	myfile.close();
}

// Change string to uppercase
template <int dim>
std::string
Elastic::ElasticProblem<dim>::to_upper(const std::string str){
	std::string out_str(str);
	for (int i = 0; i < str.size(); ++i) 
		out_str[i] = toupper(str[i]);
	return out_str;
}

#endif

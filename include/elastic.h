/* TODO

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

#include <deal.II/lac/trilinos_sparse_matrix.h>
#include <deal.II/lac/trilinos_block_sparse_matrix.h>
#include <deal.II/lac/trilinos_vector.h>
#include <deal.II/lac/trilinos_block_vector.h>
#include <deal.II/lac/trilinos_precondition.h>

#include <deal.II/lac/sparse_direct.h>
#include <deal.II/lac/sparse_ilu.h>
#include <deal.II/lac/constraint_matrix.h>

#include <fstream>
#include <sstream>
#include <string>
#include <typeinfo>

#include "parameters.h"
#include "boundary.h"
#include "coefficient.h"
#include "exact.h"
#include "preconditioner.h"
#include "rhs.h"

#ifndef ELASTIC_H
#define ELASTIC_H

#define ZERO 1.0e-8

using namespace std;
using namespace dealii;
namespace Elastic
{
// Solver structure
struct Solver
{
    // typedef SolverCG<>		inner;
    // typedef SolverCG<>		schur;
    typedef SolverGMRES<TrilinosWrappers::Vector>	inner;
    typedef SolverGMRES<>	schur;
};

template <int dim>
class ElasticProblem
{
public:
    ElasticProblem (const unsigned int degree, parameters *_par);
    void run ();

private:
    // pointer to parameter object
    parameters *par;
    // conditional outputs
    ConditionalOStream info_0, info_1, info_2;
    // timer
    TimerOutput computing_timer;
    // Change string to uppercase
    std::string to_upper(const std::string str);

    // Create matlab code
    void generate_matlab_study();
    // Write matrix to data file
    void write_matrix(const FullMatrix<double> &M, string filename );
    // Write matrix to data file
    void write_matrix(const TrilinosWrappers::SparseMatrix &M, string filename );
    // Write vector to data file
    void write_vector(const TrilinosWrappers::BlockVector &V, string filename );
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
    const unsigned int						n_blocks;
    Triangulation<dim>						triangulation;

    FESystem<dim>							fe;
    DoFHandler<dim>							dof_handler;

    ConstraintMatrix                        constraints;

    BlockSparsityPattern					sparsity_pattern;
    TrilinosWrappers::BlockSparseMatrix		system_matrix;
    TrilinosWrappers::BlockSparseMatrix		system_preconditioner; 		// preconditioner [A 0;Bt S]

    TrilinosWrappers::BlockVector			solution;
    TrilinosWrappers::BlockVector			system_rhs, load, body_force, precond_rhs;

    std_cxx1x::shared_ptr<typename Preconditioner::inner> A0_preconditioner;
    std_cxx1x::shared_ptr<typename Preconditioner::inner> A1_preconditioner;
    std_cxx1x::shared_ptr<typename Preconditioner::schur> S_preconditioner;
    //std_cxx1x::shared_ptr<TrilinosWrappers::PreconditionAMG> A_preconditioner;
};
}

/*
 ------------- IMPLEMENTATION --------------
 */
template <int dim>
Elastic::ElasticProblem<dim>::ElasticProblem (const unsigned int degree, parameters *_par)
    :
      par(_par),
      info_0(std::cout, _par->info == 0),
      info_1(std::cout, _par->info == 1),
      info_2(std::cout, _par->info == 2),
      computing_timer (info_0,
                       TimerOutput::summary,
                       TimerOutput::wall_times),
      degree (degree),
      n_blocks(dim+1),
      triangulation (Triangulation<dim>::maximum_smoothing),
      fe (FE_Q<dim>(degree+1), dim,
          FE_Q<dim>(degree), 1),
      dof_handler (triangulation)
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
    /**
     * TODO:
     * Can we use this
     * Fastens the algorithm by a small factor.
     * In ref 7 excluding this runtime takes 1100 secs
     * and with this we have 917 secs.
     */
    //    DoFRenumbering::Cuthill_McKee (dof_handler); // PROBLEM?? can we use it?

    std::vector<unsigned int> block_component (n_blocks,0);
    for(int i=0; i<n_blocks; ++i)
        block_component[i] = i;

    // DOF renumbering
    DoFRenumbering::component_wise (dof_handler, block_component);

    {
        std::vector<bool> ns_mask (dim+1, true); // NO_SLIP
        std::vector<bool> vs_mask (dim+1, true); // V_SLIP

        ns_mask[0] = true;
        ns_mask[1] = true;
        ns_mask[2] = false;

        vs_mask[0] = true;
        vs_mask[1] = false;
        vs_mask[2] = false;

        constraints.clear();
        VectorTools::interpolate_boundary_values (dof_handler,
                                                  NO_SLIP,
                                                  ZeroFunction<dim>(3),
                                                  constraints,
                                                  ns_mask);

        VectorTools::interpolate_boundary_values (dof_handler,
                                                  V_SLIP,
                                                  ZeroFunction<dim>(3),
                                                  constraints,
                                                  vs_mask);
    }
    constraints.close();

    A0_preconditioner.reset ();
    A1_preconditioner.reset ();
    S_preconditioner.reset ();

    system_matrix.clear ();
    system_preconditioner.clear ();

    std::vector<unsigned int> dofs_per_block (n_blocks);
    DoFTools::count_dofs_per_block (dof_handler, dofs_per_block, block_component);

    info_0 << "   Number of active cells: "
           << triangulation.n_active_cells()
           << ", Number of degrees of freedom: "
           << dof_handler.n_dofs()
           << " (" << dofs_per_block[0];
    for(int i=1; i<n_blocks;++i)
        info_0 << " + " << dofs_per_block[i];
    info_0 << ")" << std::endl;

    par->dofs << triangulation.n_active_cells() << "\t(" << dofs_per_block[0];
    for(int i=1; i<n_blocks;++i)
        par->dofs << " + " << dofs_per_block[i];
    par->dofs << ")" << std::endl;

    {
        BlockCompressedSimpleSparsityPattern bcsp(n_blocks, n_blocks);
        for(int i=0; i<n_blocks; ++i){
            for(int j=0; j<n_blocks; ++j){
                bcsp.block(i,j).reinit (dofs_per_block[i],
                                        dofs_per_block[j]);
            }
        }
        bcsp.collect_sizes();
        DoFTools::make_sparsity_pattern (dof_handler, bcsp, constraints, true);

        sparsity_pattern.copy_from(bcsp);
    }

    // Print sparsity pattern of the matrix
    if(par->print_matrices){
        std::ofstream out ("sparsity_pattern.1");
        sparsity_pattern.print_gnuplot (out);
    }


    system_matrix.reinit (sparsity_pattern);
    system_preconditioner.reinit (sparsity_pattern);

    solution.reinit (n_blocks);
    system_rhs.reinit (n_blocks);
    precond_rhs.reinit (n_blocks);
    load.reinit (n_blocks);
    body_force.reinit (n_blocks);

    for(int i=0; i<n_blocks; ++i){
        solution.block(i).reinit (dofs_per_block[i]);
        system_rhs.block(i).reinit (dofs_per_block[i]);
        precond_rhs.block(i).reinit (dofs_per_block[i]);
        load.block(i).reinit (dofs_per_block[i]);
        body_force.block(i).reinit (dofs_per_block[i]);
    }

    solution.collect_sizes ();
    system_rhs.collect_sizes ();
    precond_rhs.collect_sizes ();
    load.collect_sizes ();
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

    // Dummy cell matrix for preconditioner, it is allways zero
    Vector<double>      cell_rhs (dofs_per_cell),
                        cell_pre_rhs(dofs_per_cell);

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
        cell_pre_rhs    = 0;
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

        if(par->one_schur_it)
            l_S.gauss_jordan();

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
        // end assembly A preconditioner

        // printing local matrices
        if(par->print_local && first){ // print_local, first
            write_matrix(cell_matrix,"l_m");
            write_matrix(cell_precond,"l_p");
        }

        // local-to-global
        cell->get_dof_indices (local_dof_indices);

        constraints.distribute_local_to_global(cell_matrix, cell_rhs,
                                               local_dof_indices,
                                               system_matrix, system_rhs);
        constraints.distribute_local_to_global(cell_precond, cell_pre_rhs,
                                               local_dof_indices,
                                               system_preconditioner, precond_rhs);
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
    A0_preconditioner
            = std_cxx1x::shared_ptr<typename Preconditioner::inner>(new typename Preconditioner::inner());

    A1_preconditioner
            = std_cxx1x::shared_ptr<typename Preconditioner::inner>(new typename Preconditioner::inner());

    S_preconditioner
            = std_cxx1x::shared_ptr<typename Preconditioner::schur>(new typename Preconditioner::schur());

    std::vector<std::vector<bool> > constant_modes;
    std::vector<bool>  displacement_components (dim+1,false);

    // A00
    displacement_components[0] = true;
    DoFTools::extract_constant_modes (dof_handler,
                                      displacement_components,
                                      constant_modes);

    TrilinosWrappers::PreconditionAMG::AdditionalData amg_A0;
    amg_A0.constant_modes = constant_modes;

    amg_A0.elliptic = true;
    amg_A0.higher_order_elements = false;
    amg_A0.smoother_sweeps = 2;
    amg_A0.aggregation_threshold = par->threshold;


    A0_preconditioner->initialize( system_preconditioner.block(0,0),amg_A0);
    //

    //A11
    displacement_components[0] = false;
    displacement_components[1] = true;
    DoFTools::extract_constant_modes (dof_handler,
                                      displacement_components,
                                      constant_modes);

    TrilinosWrappers::PreconditionAMG::AdditionalData amg_A1;
    amg_A1.constant_modes = constant_modes;

    amg_A1.elliptic = true;
    amg_A1.higher_order_elements = false;
    amg_A1.smoother_sweeps = 2;
    amg_A1.aggregation_threshold = par->threshold;

    A1_preconditioner->initialize( system_preconditioner.block(1,1),amg_A1);
    //

    TrilinosWrappers::PreconditionAMG::AdditionalData amg_S;

    amg_S.elliptic = true;
    amg_S.higher_order_elements = false;
    amg_S.smoother_sweeps = 2;
    amg_S.aggregation_threshold = par->threshold;

    // elem-by-elem Schur
    S_preconditioner->initialize( system_preconditioner.block(2,2), amg_S);
}

template <int dim>
void
Elastic::ElasticProblem<dim>::solve ()
{

    const BlockSchurPreconditioner<typename Preconditioner::inner, // A0, schur
            typename Preconditioner::inner, // A1, schur
            typename Preconditioner::schur>
            preconditioner( system_preconditioner, *A0_preconditioner, *A1_preconditioner, *S_preconditioner, par); // system_matrix

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

    computing_timer.enter_section("DOF setup");
    setup_dofs ();
    computing_timer.exit_section("DOF setup");

    info_0 << "   Assembling ... ";
    computing_timer.enter_section("Assembling");
    assemble_system ();
    computing_timer.exit_section("Assembling");
    info_0 << " DONE" << std::endl;

    info_0 << "   AMG preconditioners ... ";
    computing_timer.enter_section("AMG preconditioners");
    setup_AMG ();
    computing_timer.exit_section("AMG preconditioners");
    info_0 << " DONE" << std::endl;

    if(par->print_matrices ){ // define print_data
        std::cout << "   ** Printing matrices in Matlab form **" << std::endl << std::flush;
        int n_blocks = dim+1;
        for(int i=0; i<n_blocks; ++i){
            for(int j=0; j<n_blocks; ++j){
                string a = "a" + std::to_string(i) + std::to_string(j);
                string p = "p" + std::to_string(i) + std::to_string(j);
                write_matrix(system_matrix.block(i,j),a);
                write_matrix(system_preconditioner.block(i,j),p);
            }
        }
        write_vector(system_rhs,"rhs");
    }

    int inv_iter = 0, schur_iter = 0;
    if(par->solve){

        info_0 << "   system solver ... ";
        computing_timer.enter_section("System solver");
        solve ();
        computing_timer.exit_section("System solver");
        info_0 << " DONE" << std::endl;

        // printing solver info
        for(unsigned int i = 0; i < par->inv_iterations.size(); i++){
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
    //    if(par->POISSON == 0.2 ){
    generate_matlab_study();
    //    }
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
    filename << "surface_values" << par->str_poisson << ".dat";
    std::ofstream detector_data(filename.str().c_str());
    std::vector<Point<dim> > detector_locations;

    Vector<double> value(3);

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

    for (unsigned int i=0 ; i<detector_locations.size(); ++i){

        VectorTools::point_value (dof_handler,
                                  solution,
                                  detector_locations[i],
                                  value);
        detector_data << (detector_locations[i](0)*par->L/1e3) << ", "
                          << (detector_locations[i](1)*par->L/1e3) << ", "
                              << value(0) << ", " << value(1) << ", " << (value(2)/par->L) << std::endl;
    }

                         info_0 << "... extracting surface values ..."
                      << " dx = " << dx << ", n = " << n
                      << std::endl;
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

    // Create matlab code
    template <int dim>
    void
            Elastic::ElasticProblem<dim>::generate_matlab_study(){
        int figure = 1;

        // Generate Matlab filename
        ostringstream tempOS;
        tempOS << "matrices" << par->str_poisson << ".m";
        string extension = tempOS.str().c_str();

        // Open file
        ofstream myfile(extension.c_str());
        if(!myfile.is_open()){
            cout << "Print_matlab: Unable to open file...";
            return;
        }

        myfile << "close all;clear;" << endl;

        if(par->print_matrices){
            int n_blocks = dim+1;
            myfile << "%loading data" << std::endl;
            for(int i=0; i<n_blocks; ++i){
                for(int j=0; j<n_blocks; ++j){
                    string a = "data_a" + std::to_string(i) + std::to_string(j)
                            + " = load_sparse('data_a" + std::to_string(i) + std::to_string(j) + ".dat');";
                    string p = "data_p" + std::to_string(i) + std::to_string(j)
                            + " = load_sparse('data_p" + std::to_string(i) + std::to_string(j) + ".dat');";
                    myfile << a << endl;
                    myfile << p << endl;
                }
            }
            myfile << "load('data_l_m.dat');"  << endl
                   << "load('data_l_p.dat');"  << endl;
            //               << "load('data_rhs.dat');"   << endl;

            myfile << "% Creating A and P matrices" << std::endl
                   << "A = [";
            for(int i=0; i<n_blocks; ++i){
                for(int j=0; j<n_blocks; ++j){
                    string a = "data_a" + std::to_string(i) + std::to_string(j) + " ";
                    myfile << a;
                }
                myfile << ";" << endl;
            }
            myfile << "];" <<endl;


            myfile << "P = [";
            for(int i=0; i<n_blocks; ++i){
                for(int j=0; j<n_blocks; ++j){
                    string p = "data_p" + std::to_string(i) + std::to_string(j) + " ";
                    myfile << p;
                }
                myfile << ";" << endl;
            }
            myfile << "];";
        }
        myfile.close();

        // Create another file
        tempOS.clear();
        tempOS.str("");
        tempOS << "SV" << par->str_poisson << ".m";
        extension = tempOS.str().c_str();

        // Open file
        myfile.open(extension.c_str());
        if(!myfile.is_open()){
            cout << "Print_matlab: Unable to open file...";
            return;
        }

        std::stringstream xplot;
        xplot << "set(gca,'XTick',";

        if(par->cases() == 2){
            xplot << "-4000:500:0)" << endl
                  << "xlabel('y(xt) ";
        }
        else if(par->cases() == 1){
            xplot << "0:5000:"<<par->x2*1e4<<")" << endl
                  << "xlabel('x(y=0) ";
        }
        else if(par->cases() == 0){
            xplot << "0:2000:"<<par->x2*1e4<<")" << endl
                  << "xlabel('x(y=0) ";
        }

        xplot << "in (Km)');" << endl;

        myfile	<< "sol = " << "load('surface_values" << par->str_poisson << ".dat');" << endl << endl
                << "x = sol(:,1);" << endl
                << "y = sol(:,2);" << endl
                << "domain = "
                << ((par->cases() == 2)?"y(:,1)":"x(:,1)") << ";" << endl;

        myfile	<< "horizontal = sol(:,3);" << endl
                << "vertical   = sol(:,4);" << endl
                << "pressure   = sol(:,5);" << endl << endl;

        myfile	<< "% Change default axes fonts." << endl
                << "set(0,'DefaultAxesFontName', 'Times New Roman');" << endl
                << "set(0,'DefaultAxesFontSize', 16);" << endl << endl;

        myfile	<< "% Change default text fonts." << endl
                << "set(0,'DefaultTextFontname', 'Times New Roman');" << endl
                << "set(0,'DefaultTextFontSize', 16);" << endl << endl;

        myfile	<< "figure("<<figure<<"); plot(domain, horizontal); %title('horizontal');" << endl
                << xplot.str()
                << "ylabel('Horizontal displacement in (m)');" << endl
                << "grid on" << endl
                << "%set(gca, 'GridLineStyle', '-');" << endl;
        figure++;

        myfile	<< "figure("<< figure <<"); plot(domain,  vertical); %title('vertical');" << endl
                << xplot.str()
                << "ylabel('Vertical displacement in (m)');" << endl
                << "grid on" << endl
                << "%set(gca, 'GridLineStyle', '-');" << endl;
        figure++;

        myfile	<< "% analytical constants, alpha = rho_i*h/rho_r, delta = (1+v)(1-2v)/(E(1-v)), 1/(2mu + lambda) = (1+v)(1-2v)/(E(1-v))" << endl
                << "v = 0.2; E = " << (par->S*par->YOUNG) << "; g = " << par->g0 << "; h = "<< par->h <<";" << endl
                << "rho_i = "<< par->rho_i << "; rho_r = " << par->rho_r << ";" << endl
                << "alpha = "<< par->alpha <<";" << endl
                << "beta  = "<< par->beta <<";" << endl
                << "delta = (1+v)*(1-2*v)/(E*(1-v));" << endl
                << "gamma = "<< par->gamma <<";" << endl << endl;

        myfile	<< "yb = " << (par->L*par->y1) << ";" << endl
                << "A1 = (1+v)*(1-2*v)/(E*(1-v))*g*rho_i*h; % rho_i" << endl
                << "mu = E/(2*(1+v));" << endl
                << "lambda = E*v/((1+v)*(1-2*v));" << endl << endl;

        myfile	<< "y = 1000*domain; %yb : abs(yb)/1000:0.0;" << endl
                << "v1 = A1.*(yb-y);" << endl
                << "v2 = alpha.*(exp(rho_r*g*delta.*yb)-exp(rho_r*g*delta.*y));" << endl
                << "p = -rho_i*g*h*gamma;" << endl << endl;

        myfile	<< "%figure("<<figure<<"); plot(y,v1,y,v2);" << endl;
        figure++;

        myfile	<< "p1 = sol(:,4);" << endl
                << "%figure("<<figure<<"); plot(domain,p1,domain,p);" << endl
                << endl;

        myfile.close();

        // Open file
        myfile.open("load_sparse.m");
        if(!myfile.is_open()){
            cout << "Print_matlab: Unable to open file...";
            return;
        }

        myfile << "function mat = load_sparse(f)" << endl
               << "% Reading output matrix from deal.II into matlab" << endl
               << "fid = fopen(f, 'rt');" << endl
               << "[m n] = fscanf(fid, '(%d,%d) %e\\n'); % has the particular format" << endl
               << "% reorganizes the data" << endl
               << "count = 1;" << endl
               << "for k = [1:3:n]" << endl
               << "    I(count) = m(k)+1;" << endl
               << "    J(count) = m(k+1)+1;" << endl
               << "    v(count) = m(k+2);" << endl
               << "    count = count + 1;" << endl
               << "end" << endl
               << "fclose(fid);" << endl
               << "% create a sparse stiffness matrix with the input" << endl
               << "mat=sparse(I,J,v);" << endl
               << "end" << endl;
        myfile.close();
    }

    // Create a matlab file with matrices
    template <int dim>
    void
            Elastic::ElasticProblem<dim>::write_matrix(const TrilinosWrappers::SparseMatrix &M, string filename ){
        int PSC = 13;
        //Printing in matlab form
        string name = "data_" + filename + ".dat";
        std::ofstream matrix (name.c_str());
        matrix << setprecision(PSC);//std::setprecision(std::numeric_limits<double>::digits10);
        M.print(matrix);
        matrix.close();
    }

    // Create a matlab file with vectors
    template <int dim>
    void
            Elastic::ElasticProblem<dim>::write_vector(const TrilinosWrappers::BlockVector &V, string filename ){
        int PSC = 13;
        //Printing in matlab form
        string name = "data_" + filename + ".dat";
        std::ofstream vector (name.c_str());
        vector << setprecision(PSC);//std::setprecision(std::numeric_limits<double>::digits10);
        V.print(vector);
        vector.close();
    }

    template <int dim>
    void
            Elastic::ElasticProblem<dim>::write_matrix(const FullMatrix<double> &M, string filename ){
        //Printing in matlab form
        string name = "data_" + filename + ".dat";
        std::ofstream matrix (name.c_str());
        M.print(matrix,10);
        matrix.close();
    }

    /*
 **************************************
 */

#endif

#include "elastic.h"
#include "parameters.h"

int main (int argc, char** argv)
{
	using namespace std;
	try
    {
		using namespace dealii;
		
		// MPI_Init (&argc, &argv);
		Utilities::MPI::MPI_InitFinalize mpi_initialization (argc, argv);
		
		// get parameters
		parameters *par = new parameters(argc,argv); 
		
		ostringstream filename;
		filename << "iterations" << par->str_poisson << ".log";
		ofstream pout(filename.str().c_str());
		
		// Attach deallog to process output
		deallog.attach(pout);
		deallog.depth_console (0);
		
		
		if(par->dimension == 2){
			Elastic::ElasticProblem<2> elastic_problem(par->degree, par);
			elastic_problem.run ();
		// }else if(par.dimension == 3){
			//ElasticProblem<3> elastic_problem(par.degree);
			//elastic_problem.run ();
		}else printf("wrong dimension");
        
    }
	catch (std::exception &exc)
    {
		std::cerr << std::endl << std::endl
		<< "----------------------------------------------------"
		<< std::endl;
		std::cerr << "Exception on processing: " << std::endl
		<< exc.what() << std::endl
		<< "Aborting!" << std::endl
		<< "----------------------------------------------------"
		<< std::endl;
		
		return 1;
    }
	catch (...)
    {
		std::cerr << std::endl << std::endl
		<< "----------------------------------------------------"
		<< std::endl;
		std::cerr << "Unknown exception!" << std::endl
		<< "Aborting!" << std::endl
		<< "----------------------------------------------------"
		<< std::endl;
		return 1;
    }
	
	return 0;
}

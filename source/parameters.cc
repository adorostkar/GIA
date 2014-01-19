#include "parameters.h"

// Just defaults, no write is allowed since it will rewrite any parameters added in the file.
parameters::parameters(const std::string f){
    paramFile = f;
	set_default();

	inv_iterations   = std::vector<unsigned int>();
	schur_iterations = std::vector<unsigned int>();

	compute_additionals();
}

// Constructor that gets the application parameters
parameters::parameters(int argc, char* argv[], const std::string f){
	// Set default file name
	paramFile = f;
	
	find_filename(argc,argv);

	set_default();

	if(fexists())
		read_Parameters();

	// If any option is added (including file name) read it.
	if(argc > 1)
		parse_command_line(argc,argv);


	inv_iterations   = std::vector<unsigned int>();
	schur_iterations = std::vector<unsigned int>();

	// if options have to be write back in the file.
	if(writeback)
		write_Parameters();

	compute_additionals();

	// If verbose, print all the parameters out
	if (verbose)
		print_values();
}

// copy constructor.
parameters::parameters(const parameters &pm){
	dimension = pm.dimension;
	degree = pm.degree;
	L = pm.L;
	U = pm.U;
	S = pm.S;
	T = pm.T;
	x1 = pm.x1;
	y1 = pm.y1;
	x2 = pm.x2;
	y2 = pm.y2;
	Ix = pm.Ix;
	h = pm.h;
    b_ice = pm.b_ice;
	b_up = pm.b_up;
	b_left = pm.b_left;
	b_right = pm.b_right;
	b_bottom = pm.b_bottom;
	refinements = pm.refinements;
	xdivisions = pm.xdivisions;
	ydivisions = pm.ydivisions;
	surf_samples = pm.surf_samples;
	YOUNG = pm.YOUNG;
	POISSON = pm.POISSON;
	ETA = pm.ETA;
	rho_i = pm.rho_i;
	rho_r = pm.rho_r;
	gravity = pm.gravity;
	load_enabled = pm.load_enabled;
	weight_enabled = pm.weight_enabled;
	adv = pm.adv;
	div = pm.div;
	precond = pm.precond;
	solve = pm.solve;
	InvMatPreTOL = pm.InvMatPreTOL;
	SchurTOL = pm.SchurTOL;
	TOL = pm.TOL;
	threshold = pm.threshold;
	one_schur_it = pm.one_schur_it;
	solver = pm.solver;
	info = pm.info;
	print_local = pm.print_local;
	print_matrices = pm.print_matrices;
	cas = pm.cas;
	paramFile = pm.paramFile;
	system_iter = pm.system_iter;
	load = pm.load;
	weight = pm.weight;
	g0 = pm.g0;
	scale1 = pm.scale1;
	scale2 = pm.scale2;
	scale3 = pm.scale3;
	alpha = pm.alpha;
	beta = pm.beta;
	delta = pm.delta;
	gamma = pm.gamma;
	adv_enabled = pm.adv_enabled;
	div_enabled = pm.div_enabled;
    inv_iterations = pm. inv_iterations;
    schur_iterations = pm.schur_iterations;
	str_poisson = pm.str_poisson;
	writeback = pm.writeback;
	verbose = pm.verbose;
}

// equality operator
parameters & parameters::operator=(const parameters &pm){
	if(this == &pm) return *this;

	dimension = pm.dimension;
	degree = pm.degree;
	L = pm.L;
	U = pm.U;
	S = pm.S;
	T = pm.T;
	x1 = pm.x1;
	y1 = pm.y1;
	x2 = pm.x2;
	y2 = pm.y2;
	Ix = pm.Ix;
	h = pm.h;
    b_ice = pm.b_ice;
	b_up = pm.b_up;
	b_left = pm.b_left;
	b_right = pm.b_right;
	b_bottom = pm.b_bottom;
	refinements = pm.refinements;
	xdivisions = pm.xdivisions;
	ydivisions = pm.ydivisions;
	surf_samples = pm.surf_samples;
	YOUNG = pm.YOUNG;
	POISSON = pm.POISSON;
	ETA = pm.ETA;
	rho_i = pm.rho_i;
	rho_r = pm.rho_r;
	gravity = pm.gravity;
	load_enabled = pm.load_enabled;
	weight_enabled = pm.weight_enabled;
	adv = pm.adv;
	div = pm.div;
	precond = pm.precond;
	solve = pm.solve;
	InvMatPreTOL = pm.InvMatPreTOL;
	SchurTOL = pm.SchurTOL;
	TOL = pm.TOL;
	threshold = pm.threshold;
	one_schur_it = pm.one_schur_it;
	solver = pm.solver;
	info = pm.info;
	print_local = pm.print_local;
	print_matrices = pm.print_matrices;
	cas = pm.cas;
	paramFile = pm.paramFile;
	system_iter = pm.system_iter;
	load = pm.load;
	weight = pm.weight;
	g0 = pm.g0;
	scale1 = pm.scale1;
	scale2 = pm.scale2;
	scale3 = pm.scale3;
	alpha = pm.alpha;
	beta = pm.beta;
	delta = pm.delta;
	gamma = pm.gamma;
	adv_enabled = pm.adv_enabled;
	div_enabled = pm.div_enabled;
    inv_iterations = pm. inv_iterations;
    schur_iterations = pm.schur_iterations;
	str_poisson = pm.str_poisson;
	writeback = pm.writeback;
	verbose = pm.verbose;
	return *this;
}

// set parameters to their defaults
void parameters::set_default(){
	dimension		= 2;
	degree			= 1;


	L				= 1.0e7;
	U				= 1.0;
	S				= 4.0e11;
	T				= 1.0;


	x1				= 0.0;
	x2				= 1.0e7;
	y1				= -4e6;
	y2				= 0.0;
	Ix				= 1e6;
	h				= 2170;

    b_ice           = LOAD;
    b_up			= FREE;
	b_left			= V_SLIP;
	b_right			= NEUMANN;
	b_bottom		= NO_SLIP;


	refinements		= 3;
	xdivisions		= 10;
	ydivisions		= 4;
	surf_samples	= 2000;


	YOUNG			= 4e11;
	POISSON			= 0.2;
	ETA 			= 100;
	rho_i 			= 917;
	rho_r 			= 3300;
	gravity 		= -9.8;


	load_enabled 	= true;
	weight_enabled 	= false;
	adv 			= 1;
	div 			= 1;

	precond 		= 0;
	solve 			= true;
	InvMatPreTOL 	= 1e-2;
	SchurTOL 		= 0.1;
	TOL 			= 1e-7;
	threshold 		= 0.02;
	one_schur_it 	= false;

	solver 			= 1;
	info 			= false;

	print_local 	= false;
	print_matrices 	= false;

	cas = 0;

	writeback = false;
	verbose = false;
}

void parameters::read_Parameters(){
	// create and open an archive for input
    std::ifstream ifs(paramFile.c_str());
    // Read parameters from file
    std::string iptSt;
    if(ifs.is_open()){
    	while(!ifs.eof()){
        	ifs >> iptSt;
    		if(iptSt == "dimension")
    			ifs >> dimension;
    		else if(iptSt == "degree")
    			ifs >> degree;
    		else if(iptSt == "L")
    			ifs >> L;
			else if(iptSt == "U")
				ifs >> U;
			else if(iptSt == "S")
				ifs >> S;
			else if(iptSt == "T")
				ifs >> T;
			else if(iptSt == "x1")
				ifs >> x1;
			else if(iptSt == "y1")
				ifs >> y1;
			else if(iptSt == "x2")
				ifs >> x2;
			else if(iptSt == "y2")
				ifs >> y2;
			else if(iptSt == "Ix")
				ifs >> Ix;
			else if(iptSt == "h")
				ifs >> h;
            else if(iptSt == "b_ice"){
                std::string s;
                ifs >> s;
                b_ice = str2boundary(s);
            }else if(iptSt == "b_up"){
				std::string s;
				ifs >> s;
				b_up = str2boundary(s);
			}else if(iptSt == "b_left"){
				std::string s;
				ifs >> s;
				b_left = str2boundary(s);
			}else if(iptSt == "b_right"){
				std::string s;
				ifs >> s;
				b_right = str2boundary(s);
			}else if(iptSt == "b_bottom"){
				std::string s;
				ifs >> s;
				b_bottom = str2boundary(s);
			}else if(iptSt == "refinements")
				ifs >> refinements;
			else if(iptSt == "xdivisions")
				ifs >> xdivisions;
			else if(iptSt == "ydivisions")
				ifs >> ydivisions;
			else if(iptSt == "surf_samples")
				ifs >> surf_samples;
			else if(iptSt == "YOUNG")
				ifs >> YOUNG;
			else if(iptSt == "POISSON")
				ifs >> POISSON;
			else if(iptSt == "ETA")
				ifs >> ETA;
			else if(iptSt == "rho_i")
				ifs >> rho_i;
			else if(iptSt == "rho_r")
				ifs >> rho_r;
			else if(iptSt == "gravity")
				ifs >> gravity;
			else if(iptSt == "load_enabled")
				ifs >> load_enabled;
			else if(iptSt == "weight_enabled")
				ifs >> weight_enabled;
			else if(iptSt == "adv")
				ifs >> adv;
			else if(iptSt == "div")
				ifs >> div;
			else if(iptSt == "precond")
				ifs >> precond;
			else if(iptSt == "solve")
				ifs >> solve;
			else if(iptSt == "InvMatPreTOL")
				ifs >> InvMatPreTOL;
			else if(iptSt == "SchurTOL")
				ifs >> SchurTOL;
			else if(iptSt == "TOL")
				ifs >> TOL;
			else if(iptSt == "threshold")
				ifs >> threshold;
			else if(iptSt == "one_schur_it")
				ifs >> one_schur_it;
			else{
                std::string temp;
                std::getline(ifs, temp);
            }
        }
        ifs.close();
    }

}

void parameters::write_Parameters(){
	using namespace std;
	// create and open an archive for input
    ofstream ofs;
    ofs.open(paramFile.c_str(), ios::trunc);
    // write parameters to file
    if(ofs.is_open()){
    	ofs << "## Problem dimension." << endl;
    	ofs << "dimension "		<< dimension << endl;

    	ofs << "\n## Degree of the polynomial basis functions." << endl;
    	ofs << "degree "		<< degree << endl;

    	ofs << "\n## Scaling value for length." << endl;
    	ofs << "L "				<< L 					<< endl; 

    	ofs << "\n## Scaling value for displacemnet." << endl;
    	ofs << "U "				<< U 					<< endl;

    	ofs << "\n## Scaling value for stress." << endl;
    	ofs << "S "				<< S 					<< endl;

    	ofs << "\n## Scaling value for time." << endl;
    	ofs << "T "				<< T 					<< endl;

        ofs << "\n## X coordinate of the bottom left side of the domain." << endl;
    	ofs << "x1 "			<< x1 					<< endl;

        ofs << "\n## Y coordinate of the bottom left side of the domain." << endl;
    	ofs << "y1 "			<< y1 					<< endl;

        ofs << "\n## X coordinate of the top right side of the domain." << endl;
    	ofs << "x2 "			<< x2 					<< endl;

        ofs << "\n## Y coordinate of the top right side of the domain." << endl;
    	ofs << "y2 "			<< y2					<< endl;

        ofs << "\n## Width of the load(Ice)." << endl;
    	ofs << "Ix "			<< Ix					<< endl;

    	ofs << "\n## Thickness of ice used as load." << endl;
    	ofs << "h "				<< h 					<< endl;

        ofs << "\n## Boundary type on the ice." << endl;
        ofs << "b_ice "			<< boudary2str(b_ice)	<< endl;

    	ofs << "\n## Boundary type on the top wall." << endl;
    	ofs << "b_up "			<< boudary2str(b_up)	<< endl;

    	ofs << "\n## Boundary type on the left wall." << endl;
    	ofs << "b_left "		<< boudary2str(b_left)	<< endl;

    	ofs << "\n## Boundary type on the write wall." << endl;
    	ofs << "b_right "		<< boudary2str(b_right)	<< endl;

    	ofs << "\n## Boundary type on the bottom wall." << endl;
    	ofs << "b_bottom "		<< boudary2str(b_bottom)<< endl;

    	ofs << "\n## Number of refinement levels of the mesh." << endl;
    	ofs << "refinements "	<< refinements			<< endl;

        ofs << "\n## Number of initial divisions in x direction." << endl;
    	ofs << "xdivisions "	<< xdivisions			<< endl;

        ofs << "\n## Number of initial divisions in y direction." << endl;
    	ofs << "ydivisions "	<< ydivisions			<< endl;

    	ofs << "\n## The number of samples at the surface." << endl;
    	ofs << "surf_samples " 	<< surf_samples			<< endl;

    	ofs << "\n## Young's modulus." << endl;
    	ofs << "YOUNG "			<< YOUNG 				<< endl;

    	ofs << "\n## POISSON ratio." << endl;
    	ofs << "POISSON "		<< POISSON 				<< endl;

    	ofs << "\n## Description." << endl;
    	ofs << "ETA "			<< ETA 					<< endl;

    	ofs << "\n## Density of Ice used as load." << endl;
    	ofs << "rho_i "			<< rho_i				<< endl;

    	ofs << "\n## Density of rock (inner earth)." << endl;
    	ofs << "rho_r "			<< rho_r				<< endl;

    	ofs << "\n## Gravity force of the earth." << endl;
    	ofs << "gravity "		<< gravity				<< endl;

    	ofs << "\n## Description." << endl;
    	ofs << "load_enabled " 	<< load_enabled			<< endl;

    	ofs << "\n## Description." << endl;
    	ofs << "adv "			<< adv					<< endl;

    	ofs << "\n## Description." << endl;
    	ofs << "div "			<< div					<< endl;

    	ofs << "\n## {0 = diag(A), 1 = A} Whether to rewrite block A_11 in the preconditioner." << endl;
    	ofs << "precond "		<< precond 				<< endl;

    	ofs << "\n## {0,1} specifying whether to solve the system or not." << endl;
    	ofs << "solve "			<< solve 				<< endl;

    	ofs << "\n## The tolerances of the inner solvers for computing the effect of inv(A)." << endl;
    	ofs << "InvMatPreTOL "	<< InvMatPreTOL			<< endl;

    	ofs << "\n## The tolerances of the inner solvers for computing the effect of \tilde(A)." << endl;
    	ofs << "SchurTOL "		<< SchurTOL 			<< endl;

    	ofs << "\n## The system tolerance controlling the relative L2-error in the global scheme." << endl;
    	ofs << "TOL "			<< TOL 					<< endl;

    	ofs << "\n## AMG preconditioner threshold." << endl;
    	ofs << "threshold "		<< threshold			<< endl;

    	ofs << "\n## {0,1} Whether to do only one Schure iteration on the matrix." << endl;
    	ofs << "one_schur_it " << one_schur_it			<< endl;

    }
    ofs.close();
}

void parameters::find_filename(int argc, char* argv[]){
	for (int i = 0; i < argc; ++i)
	{
		if(strcmp(argv[i],"--file") == 0){
			paramFile = argv[i+1];
			return;
		}

	}
}

void parameters::parse_command_line(int argc, char* argv[]){
	int next_option;
	bool flag = false;
	verbose = false;
    do {
		next_option = getopt_long (argc, argv, short_options,
								   long_options, NULL);
        switch (next_option){
			case '0': // --footing
				x1			= 0.0;
				y1			= -4e6;
				x2			= 1e7;
				y2			= 0.0;
				Ix			= 1e6;
				xdivisions	= 10;
				ydivisions	= 4;
				cas = 0;
				break;
			case '1': // --Large, GIA footing, Large domain
				x1			= 0.0;
				y1			= -4e6;
				x2			= 25e6;
				y2			= 0.0;
				Ix			= 1e6;
				xdivisions	= 25;
				ydivisions	= 4;
				flag = true;
				cas = 1;
				break;
			case '2': // --uniform, GIA uniform load
				x1			= 0.0;
				y1			=-4e6;
				x2			= 4e6;
				y2			= 0.0;
				Ix			= 4e6;
				xdivisions	= 1;
				ydivisions	= 1;
				
                b_ice				= LOAD;
                b_up				= FREE;
				b_left				= V_SLIP;
				b_right				= V_SLIP;
				b_bottom			= NO_SLIP;

				cas = 2;
				break;
			case '3':
				paramFile=argv[optind];
				break;
			case '4':
				writeback = true;
				break;
            case '5':
                dimension = atoi(optarg);
                break;
			case 'a':   // -a or --adv
				adv = atof(optarg);
				break;	
			case 'd':   // -d or --div
				div = atof(optarg);
				break;
			case 'e':   // -e or --elastic
				adv = 0;
				div = 0;
				break;
			case 'f':   // -f or --surf_samples
				surf_samples = atoi(optarg);
				break;
			case 'g':   // -g or --precond
				precond = atoi(optarg);
				break;
			case 'h':   /* -h or --help */
				/* User has requested usage information.  Print it to standard
				 output, and exit with exit code zero (normal termination).  */
				printf("help - file not written.\n");
				print_usage (0);
				exit(0);
				break;
			case 'i':   // -i or --inv_tol
				InvMatPreTOL = atof(optarg);
				break;
			case 'k':   // -k or --info
				info = atoi(optarg);
				break;
			case 'm':   // -m or --matlab
				print_local    = 1;
				print_matrices = 1;
				break;
			case 'p':   // -p or --poisson
				POISSON = atof(optarg);
				break;
			case 'r':   // -r or --refinements
				refinements = atoi(optarg);
				break;
			case 's':   // -s or --schur_tol
				SchurTOL = atof(optarg);
				break;
			case 't':   // -t or --tol
				TOL = atof(optarg);
				break;
			case 'u':   // -u or --unique_schur one_schur_it
				one_schur_it = 1;
				break;	
			case 'v':   // -v or --verbose
				verbose = true;
				break;
			case 'w':   // -w or --weight
				weight_enabled = 1;
				break;
			case 'y':   // -y or --young
				YOUNG = atof(optarg);
				break;
			case 'z':   // -z or --threshold
				threshold = atof(optarg);
				break;
			case '?':   /* The user specified an invalid option.  */
				/* Print usage information to standard error, and exit with exit
				 code one (indicating abnormal termination).  */
				printf("?");
				print_usage (1);
				exit(1);
				break;
			case -1:    /* Done with options.  */
				break;	
			default:
			    /* Something else: unexpected. Don't process  */
				// printf("-- unexpected --\n");
				// abort ();
				break;
		}
	}
	while (next_option != -1);
}

void parameters::print_usage(int exitNum){
	using namespace std;
	stringstream msg;
	int firstCol = 6;
	int secondCol = 20;
	int thirdCol = 20;
	msg << left << setw(firstCol)	<< setfill(' ') << "-#" 
		<< left << setw(secondCol)	<< setfill(' ') << "--[####]" 
		<< left << setw(thirdCol)	<< setfill(' ') << "[arguments]" 
		<< "[ Description ]" << endl;

	msg << left << setw(firstCol)	<< setfill(' ') << "" 
		<< left << setw(secondCol)  << setfill(' ') << "--file" 
		<< left << setw(thirdCol)   << setfill(' ') << "paramFile" 
		<< "Option file"  << endl;

	msg << left << setw(firstCol)	<< setfill(' ') << "" 
		<< left << setw(secondCol)	<< setfill(' ') << "--footing" 
		<< left << setw(thirdCol)	<< setfill(' ') << "void" 
		<< "GIA footing"  << endl;

	msg << left << setw(firstCol)	<< setfill(' ') << "" 
		<< left << setw(secondCol)	<< setfill(' ') << "--large" 
		<< left << setw(thirdCol)	<< setfill(' ') << "void" 
		<< "GIA footing(Large domain)"  << endl;

	msg << left << setw(firstCol)	<< setfill(' ') << "" 
		<< left << setw(secondCol)	<< setfill(' ') << "--uniform" 
		<< left << setw(thirdCol)	<< setfill(' ') << "void" 
		<< "GIA uniform load"  << endl;

    msg << left << setw(firstCol)	<< setfill(' ') << ""
        << left << setw(secondCol)	<< setfill(' ') << "--dim"
        << left << setw(thirdCol)	<< setfill(' ') << "integer"
        << "Problem dimension"  << endl;

    msg << left << setw(firstCol)	<< setfill(' ') << ""
        << left << setw(secondCol)	<< setfill(' ') << "--writeback"
        << left << setw(thirdCol)	<< setfill(' ') << "void"
        << "Write parameters back to the file or in the default file."  << endl;

	msg << left << setw(firstCol)	<< setfill(' ') << "-a" 
		<< left << setw(secondCol)	<< setfill(' ') << "--adv" 
		<< left << setw(thirdCol)	<< setfill(' ') << "double" 
		<< "Overwrite adv"  << endl;

    msg << left << setw(firstCol)	<< setfill(' ') << "-b"
        << left << setw(secondCol)	<< setfill(' ') << "--boundary_cond"
        << left << setw(thirdCol)	<< setfill(' ') << "boundary_Type"
        << "The boundary condition at the domain surface."  << endl;

	msg << left << setw(firstCol)	<< setfill(' ') << "-d" 
		<< left << setw(secondCol)	<< setfill(' ') << "--div" 
		<< left << setw(thirdCol)	<< setfill(' ') << "double" 
		<< "Overwrite div"  << endl;

    msg << left << setw(firstCol)	<< setfill(' ') << "-e"
        << left << setw(secondCol)	<< setfill(' ') << "--elastic"
        << left << setw(thirdCol)	<< setfill(' ') << "void"
        << "Set adv and div to zero."  << endl;

    msg << left << setw(firstCol)	<< setfill(' ') << "-f"
        << left << setw(secondCol)	<< setfill(' ') << "--surf_samples"
        << left << setw(thirdCol)	<< setfill(' ') << "int"
        << "The number of samples at the surface."  << endl;

    msg << left << setw(firstCol)	<< setfill(' ') << "-g"
        << left << setw(secondCol)	<< setfill(' ') << "--precond"
        << left << setw(thirdCol)	<< setfill(' ') << "int{0,1}"
        << "Overwrite block 11 in the preconditioner(0=diag(A),1=A)"  << endl;

    msg << left << setw(firstCol) 	<< setfill(' ') << "-h"
        << left << setw(secondCol)	<< setfill(' ') << "--help"
        << left << setw(thirdCol)	<< setfill(' ')  << "void"
        << "Dispalay usage information"  << endl;

    msg << left << setw(firstCol)	<< setfill(' ') << "-i"
        << left << setw(secondCol)	<< setfill(' ') << "--inv_tol"
        << left << setw(thirdCol)	<< setfill(' ') << "double"
        << "Overwrite InvMatPreTOL"	<< endl;

    msg << left << setw(firstCol)	<< setfill(' ') << "-k"
        << left << setw(secondCol)	<< setfill(' ') << "--info"
        << left << setw(thirdCol)	<< setfill(' ') << "void"
        << "Print out GIA parameters."  << endl;

    msg << left << setw(firstCol)	<< setfill(' ') << "-m"
        << left << setw(secondCol)	<< setfill(' ') << "--matlab"
        << left << setw(thirdCol)	<< setfill(' ') << "void"
        << "Enable Matlab printing of matrices."  << endl;

    msg << left << setw(firstCol)	<< setfill(' ') << "-o"
        << left << setw(secondCol)	<< setfill(' ') << "--output"
        << left << setw(thirdCol)	<< setfill(' ') << "paramFile"
        << "Write output to file"	<< endl;

	msg << left << setw(firstCol)	<< setfill(' ') << "-p" 
		<< left << setw(secondCol)	<< setfill(' ') << "--poisson" 
		<< left << setw(thirdCol)	<< setfill(' ') << "double" 
		<< "poisson ratio"  << endl;

    msg << left << setw(firstCol)	<< setfill(' ') << "-r"
        << left << setw(secondCol)	<< setfill(' ') << "--refinements"
        << left << setw(thirdCol)	<< setfill(' ') << "int"
        << "Overwrite refinements"  << endl;

    msg << left << setw(firstCol)	<< setfill(' ') << "-s"
        << left << setw(secondCol)	<< setfill(' ') << "--schur_tol"
        << left << setw(thirdCol)	<< setfill(' ') << "double"
        << "Overwrite SchurTOL" 	<< endl;

    msg << left << setw(firstCol)	<< setfill(' ') << "-t"
        << left << setw(secondCol)	<< setfill(' ') << "--system_tol"
        << left << setw(thirdCol)	<< setfill(' ') << "double"
        << "Overwrite TOL"  << endl;

    msg << left << setw(firstCol)	<< setfill(' ') << "-u"
        << left << setw(secondCol)	<< setfill(' ') << "--unique_schur"
        << left << setw(thirdCol)	<< setfill(' ') << "void"
        << "Enable only one Schure iteration."  << endl;

    msg << left << setw(firstCol)	<< setfill(' ') << "-v"
        << left << setw(secondCol)	<< setfill(' ') << "--verbose"
        << left << setw(thirdCol)	<< setfill(' ') << "void"
        << "Verbose runtime with details."  << endl;

    msg << left << setw(firstCol)	<< setfill(' ') << "-w"
        << left << setw(secondCol)	<< setfill(' ') << "--weight"
        << left << setw(thirdCol)	<< setfill(' ') << "void"
        << "Enable weight as body force."  << endl;

	msg << left << setw(firstCol)	<< setfill(' ') << "-y" 
		<< left << setw(secondCol)	<< setfill(' ') << "--young" 
		<< left << setw(thirdCol)	<< setfill(' ') << "double" 
		<< "yound modulus"  << endl;

	msg << left << setw(firstCol)	<< setfill(' ') << "-z" 
		<< left << setw(secondCol)	<< setfill(' ') << "--threshold" 
		<< left << setw(thirdCol)	<< setfill(' ') << "double" 
		<< "AMG threshold"  << endl;

    if(!exitNum){
	    std::cout << msg.str();
	}else{
		std::cerr << msg.str();
	}
}

void parameters::print_variables(){
	using namespace std;
	int emptySpace = 5;
	int firstCol = 35;
	int secondCol = 25;
	int thirdCol = 15;
	int forthCol = 15;
	std::ostringstream sTOL;
	one_schur_it? (sTOL << "x") : (sTOL << SchurTOL);
	
	// string Title
	std::ostringstream title, s_precond, outStr;
	switch(cas){
		case 0:
			title << "Case: Footing";
			s_precond << "P_00: diag(A)";
			break;
		case 1:
			title << "Case: Footing, Large domain";
			s_precond << "P_00: A";
			break;
		case 2:
			title << "Case: Uniform load";
			s_precond << "P_00: A";
			break;
	}
	
	if(info == 0){
		outStr << title.str() << ", " << s_precond.str() << std::endl;
		outStr << left << setw(emptySpace)<< setfill(' ') << "" 
				 << left << setw(firstCol)	<< setfill(' ') << "YOUNG = " 
				 << left << setw(secondCol) << setfill(' ') << YOUNG*S 
				 << left << setw(thirdCol)	<< setfill(' ') << "POISSON = " 
				 << left << setw(forthCol)	<< setfill(' ') << POISSON 
				 << endl;

		outStr << left << setw(emptySpace)<< setfill(' ') << "" 
				 << left << setw(firstCol)	<< setfill(' ') << "ETA = " 
				 << left << setw(secondCol) << setfill(' ') << ETA 
				 << left << setw(thirdCol)	<< setfill(' ') << "" 
				 << left << setw(forthCol)	<< setfill(' ') << "" 
				 << endl;

		outStr << left << setw(emptySpace)<< setfill(' ') << "" 
				 << left << setw(firstCol)	<< setfill(' ') << "scale1(L^2/(SU)) = " 
				 << left << setw(secondCol) << setfill(' ') << scale1 
				 << left << setw(thirdCol)	<< setfill(' ') << "scale2(L/SU) = " 
				 << left << setw(forthCol)	<< setfill(' ') << scale2 
				 << endl;

		outStr << left << setw(emptySpace)<< setfill(' ') << "" 
				 << left << setw(firstCol)	<< setfill(' ') << "scale3 (L/S)*rho_r*g = " 
				 << left << setw(secondCol) << setfill(' ') << scale3 
				 << left << setw(thirdCol)	<< setfill(' ') << "" 
				 << left << setw(forthCol)	<< setfill(' ') << "" 
				 << endl;

		outStr << left << setw(emptySpace)<< setfill(' ') << "" 
				 << left << setw(firstCol)	<< setfill(' ') << "load" 
				 << left << setw(secondCol) << setfill(' ') << load/scale2 
				 << left << setw(thirdCol)	<< setfill(' ') << "load_enabled" 
				 << left << setw(forthCol)	<< setfill(' ') << load_enabled
				 << endl;

		outStr << left << setw(emptySpace)<< setfill(' ') << "" 
				 << left << setw(firstCol)	<< setfill(' ') << "weight" 
				 << left << setw(secondCol) << setfill(' ') << weight 
				 << left << setw(thirdCol)	<< setfill(' ') << "weight_enabled" 
				 << left << setw(forthCol)	<< setfill(' ') << weight_enabled
				 << endl;

		outStr << left << setw(emptySpace)<< setfill(' ') << "" 
				 << left << setw(firstCol)	<< setfill(' ') << "adv" 
				 << left << setw(secondCol) << setfill(' ') << adv*scale3 
				 << left << setw(thirdCol)	<< setfill(' ') << "adv_enabled" 
				 << left << setw(forthCol)	<< setfill(' ') << adv_enabled
				 << endl;

		outStr << left << setw(emptySpace)<< setfill(' ') << "" 
				 << left << setw(firstCol)	<< setfill(' ') << "div" 
				 << left << setw(secondCol) << setfill(' ') << div*scale3 
				 << left << setw(thirdCol)	<< setfill(' ') << "div_enabled" 
				 << left << setw(forthCol)	<< setfill(' ') << div_enabled
				 << endl;

		outStr << left << setw(emptySpace)<< setfill(' ') << "" 
				 << left << setw(firstCol)	<< setfill(' ') << "TOL(Inner,Schur)" 
				 << left << setw(secondCol) << setfill(' ') << TOL << "(" << InvMatPreTOL << "," <<sTOL.str() << ")" 
				 << left << setw(thirdCol)	<< setfill(' ') << "" 
				 << left << setw(forthCol)	<< setfill(' ') << "" 
				 << endl;

		outStr << left << setw(emptySpace)<< setfill(' ') << "" 
				 << left << setw(firstCol)	<< setfill(' ') << "AMG threshold" 
				 << left << setw(secondCol) << setfill(' ') << threshold 
				 << left << setw(thirdCol)	<< setfill(' ') << "" 
				 << left << setw(forthCol)	<< setfill(' ') << "" 
				 << endl;

		outStr << left << setw(emptySpace)<< setfill(' ') << "" 
				 << left << setw(firstCol)	<< setfill(' ') << "refinements" 
				 << left << setw(secondCol) << setfill(' ') << refinements 
				 << left << setw(thirdCol)	<< setfill(' ') << "" 
				 << left << setw(forthCol)	<< setfill(' ') << "" 
				 << endl;

		outStr << left << setw(emptySpace)<< setfill(' ') << "" 
				 << "Scaled dimensions (x1, y1, x2, y2), (Ix,h) = (" << x1 << ", " << y1 << ", " << x2 << ", " << y2 << "), (" << Ix << ", " << h/L << ")" << std::endl;
	}else if(info == 1){
		firstCol = 15;
		outStr << left << setw(firstCol)	<< setfill(' ') << "case" 
				 << left << setw(firstCol)	<< setfill(' ') << "precond"
				 << left << setw(firstCol)	<< setfill(' ') << "adv_enabled"
				 << left << setw(firstCol)	<< setfill(' ') << "div_enabled"
				 << left << setw(firstCol)	<< setfill(' ') << "refinements"
				 << left << setw(firstCol)	<< setfill(' ') << "POISSON" 
				 << endl;

		outStr << left << setw(firstCol)	<< setfill(' ') << cas
				 << left << setw(firstCol)	<< setfill(' ') << precond
				 << left << setw(firstCol)	<< setfill(' ') << adv_enabled
				 << left << setw(firstCol)	<< setfill(' ') << div_enabled
				 << left << setw(firstCol)	<< setfill(' ') << refinements
				 << left << setw(firstCol)	<< setfill(' ') << POISSON 
				 << endl;
	}else if(info == 2){ // invTOL test
		firstCol = 15;
		outStr << left << setw(firstCol)	<< setfill(' ') << "case" 
				 << left << setw(firstCol)	<< setfill(' ') << "precond"
				 << left << setw(firstCol)	<< setfill(' ') << "adv_enabled"
				 << left << setw(firstCol)	<< setfill(' ') << "div_enabled"
				 << left << setw(firstCol)	<< setfill(' ') << "refinements"
				 << left << setw(firstCol)	<< setfill(' ') << "POISSON" 
				 << left << setw(firstCol)	<< setfill(' ') << "InvMatPreTOL"
				 << left << setw(firstCol)	<< setfill(' ') << "SchurTOL"
				 << std::endl;

		outStr << left << setw(firstCol)	<< setfill(' ') << cas
				 << left << setw(firstCol)	<< setfill(' ') << precond
				 << left << setw(firstCol)	<< setfill(' ') << adv_enabled
				 << left << setw(firstCol)	<< setfill(' ') << div_enabled
				 << left << setw(firstCol)	<< setfill(' ') << refinements
				 << left << setw(firstCol)	<< setfill(' ') << POISSON 
				 << left << setw(firstCol)	<< setfill(' ') << InvMatPreTOL 
				 << left << setw(firstCol)	<< setfill(' ') << SchurTOL 
				 << endl;
	}

	std::cout << outStr.str();
}

void parameters::print_values(){
	int shift_num = 30;
	using namespace std;

	cout << "====================== Parameters =====================" << endl;
	cout << left << setw(shift_num) << setfill('.') << "Problem Dimension " 
			<< dimension << endl;

	cout << left << setw(shift_num) << setfill('.') << "Degree "
			<< degree << endl;


	cout << left << setw(shift_num) << setfill('.') << "L "
			<< L << endl;
	cout << left << setw(shift_num) << setfill('.') << "U "
			<< U << endl;
	cout << left << setw(shift_num) << setfill('.') << "S "
			<< S << endl;
	cout << left << setw(shift_num) << setfill('.') << "T "
			<< T << endl;

	cout << left << setw(shift_num) << setfill('.') << "x1 "
			<< x1 << endl;
	cout << left << setw(shift_num) << setfill('.') << "x2 "
			<< x2 << endl;
	cout << left << setw(shift_num) << setfill('.') << "y1 "
			<< y1 << endl;
	cout << left << setw(shift_num) << setfill('.') << "y2 "
			<< y2 << endl;
	cout << left << setw(shift_num) << setfill('.') << "Ix "
			<< Ix << endl;
	cout << left << setw(shift_num) << setfill('.') << "h "
			<< h << endl;


    cout << left << setw(shift_num) << setfill('.') << "Under ice boundary boundary "
            << boudary2str(b_ice) << endl;
	cout << left << setw(shift_num) << setfill('.') << "Upper wall boundary "
			<< boudary2str(b_up) << endl;
	cout << left << setw(shift_num) << setfill('.') << "Left wall boundary "
			<< boudary2str(b_left) << endl;
	cout << left << setw(shift_num) << setfill('.') << "Right wall boundary "
			<< boudary2str(b_right) << endl;
	cout << left << setw(shift_num) << setfill('.') << "Bottom wall boundary "
			<< boudary2str(b_bottom) << endl;

	cout << left << setw(shift_num) << setfill('.') << "refinements "
			<< refinements << endl;
	cout << left << setw(shift_num) << setfill('.') << "xdivisions "
			<< xdivisions << endl;
	cout << left << setw(shift_num) << setfill('.') << "ydivisions "
			<< ydivisions << endl;
	cout << left << setw(shift_num) << setfill('.') << "surf_samples "
			<< surf_samples << endl;


	cout << left << setw(shift_num) << setfill('.') << "YOUNG "
			<< YOUNG << endl;
	cout << left << setw(shift_num) << setfill('.') << "POISSON "
			<< POISSON << endl;
	cout << left << setw(shift_num) << setfill('.') << "ETA "
			<< ETA << endl;
	cout << left << setw(shift_num) << setfill('.') << "rho_i "
			<< rho_i << endl;
	cout << left << setw(shift_num) << setfill('.') << "rho_r "
			<< rho_r << endl;
	cout << left << setw(shift_num) << setfill('.') << "gravity "
			<< gravity << endl;


	cout << left << setw(shift_num) << setfill('.') << "load_enabled "
			<< load_enabled << endl;
	cout << left << setw(shift_num) << setfill('.') << "weight_enabled "
			<< weight_enabled << endl;
	cout << left << setw(shift_num) << setfill('.') << "adv "
			<< adv << endl;
	cout << left << setw(shift_num) << setfill('.') << "div "
			<< div << endl;

	cout << left << setw(shift_num) << setfill('.') << "precond "
			<< precond << endl;
	cout << left << setw(shift_num) << setfill('.') << "solve "
			<< solve << endl;
	cout << left << setw(shift_num) << setfill('.') << "InvMatPreTOL "
			<< InvMatPreTOL << endl;
	cout << left << setw(shift_num) << setfill('.') << "SchurTOL "
			<< SchurTOL << endl;
	cout << left << setw(shift_num) << setfill('.') << "TOL "
			<< TOL << endl;
	cout << left << setw(shift_num) << setfill('.') << "threshold "
			<< threshold << endl;
	cout << left << setw(shift_num) << setfill('.') << "one_schur_it "
			<< one_schur_it << endl;

	cout << left << setw(shift_num) << setfill('.') << "solver "
			<< solver << endl;
	cout << left << setw(shift_num) << setfill('.') << "info "
			<< info << endl;

	cout << left << setw(shift_num) << setfill('.') << "print_local "
			<< print_local << endl;
	cout << left << setw(shift_num) << setfill('.') << "print_matrices "
			<< print_matrices << endl;

	cout << "====================================================" << endl;
}

boundary_Type parameters::str2boundary(std::string tempSt){
	boundary_Type bt = NEUMANN;
	if(tempSt == std::string("NEUMANN"))
		bt = NEUMANN;
	else if(tempSt == std::string("NO_SLIP"))
    	bt = NO_SLIP;
	else if(tempSt == std::string("V_SLIP"))
		bt = V_SLIP;
	else if(tempSt == std::string("LOAD"))
    	bt = LOAD;
    else if(tempSt == std::string("FREE"))
        bt = FREE;

	return bt;
}

std::string parameters::boudary2str(boundary_Type bt){
	std::string tempSt;
	switch (bt){
		case NEUMANN:
	    	tempSt = "NEUMANN";
	    	break;
		case NO_SLIP:
	    	tempSt = "NO_SLIP";
	    	break;
	    case V_SLIP:
	    	tempSt = "V_SLIP";
	    	break;
	    case LOAD:
	    	tempSt = "LOAD";
	    	break;
        case FREE:
            tempSt = "FREE";
            break;
	}
	return tempSt;
}

bool parameters::fexists(){
	struct stat buf;
	if (stat(paramFile.c_str(), &buf) != -1){
        return true;
    }
    return false;
}

void parameters::compute_additionals(){

	// Scale parameters
	YOUNG = YOUNG/S;
	x1 = x1/L;
	x2 = x2/L;
	y1 = y1/L;
	y2 = y2/L;
	Ix = Ix/L;
	
	// converting par.POISSON in a string and removing the '.'
	std::ostringstream strs;
	strs << POISSON;
	str_poisson = strs.str();
	std::replace(str_poisson.begin(), str_poisson.end(), '.', '_');
	// end converting
	
	// computing constants
	double E = YOUNG*S, v = POISSON, g = fabs(gravity);
	
	alpha = rho_i*h/rho_r;
	beta  = E*(1.0-2*v)/(4*v*(1+v));
	delta = (1.0+v)*(1.0-2*v)/(E*(1-v));
	gamma = 2*v*(1+v)/(E*(1-v));
	g0	  = fabs(gravity);
	
	// scaling
	scale1 = L*L/(S*U);
	scale2 = L/(S*U);
	scale3 = scale2*rho_r*gravity;
	
	weight = (weight_enabled)?(scale1*rho_r*gravity):0.0;
	
	load = (load_enabled)?(scale2*rho_i*gravity*h):0.0;
	//par.adv  = (par.adv_enabled )?(par.scale2*par.rho_r*par.gravity) :0.0;
	//par.div  = (par.div_enabled )?(par.scale2*par.rho_r*par.gravity) :0.0;
	
	adv_enabled = (int)fabs(adv); // to be removed
	div_enabled = (int)fabs(div); // to be removed
	
	// printf("**adv = %f, div = %f",par.adv,par.div);
}






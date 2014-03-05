
#include "parameters.h"
bool parameters::instanceFlag = false;
parameters* parameters::singlton = NULL;

using namespace std;
namespace po = boost::program_options;
parameters::~parameters() {
    instanceFlag = false;
    singlton = NULL;
}

parameters* parameters::getInstance() {
    if(!instanceFlag){
        singlton = new parameters();
        instanceFlag = true;
    }

    return singlton;
}

parameters* parameters::getInstance(int _argc, char* _argv[]) {
    if(!instanceFlag){
        singlton = new parameters(_argc,_argv);
        instanceFlag = true;
    }
    return singlton;
}

// The sequence of this method is important
// First default file is red, then other file and at last
// the command line options
parameters::parameters(int argc, char* argv[]):
    general("General configurations"),
    vars("Problem variables"),
    file_options("All file options"),
    cmdLine_options("Allowed options")
{
    create_options();

    default_file = "default.cfg";
    if(!fexists(default_file)){
        std::cerr << "default.cfg file must always be available!\nSample file is created.\n";
        write_sample_file();
        exit(1);
    }

    // Command line options setup
    po::store(po::parse_command_line(argc, argv, cmdLine_options), vm);

    if(vm.count("file"))
        param_file = vm["file"].as<std::string>();

    ifstream ifsd;
    if(fexists(param_file)){
        ifsd.open(param_file);
        po::store(po::parse_config_file(ifsd, file_options), vm);
        ifsd.close();
    }

    ifsd.open(default_file);
    po::store(po::parse_config_file(ifsd, file_options), vm);
    ifsd.close();

    po::notify(vm);

    setup_variables(vm);
    validate_options();
    compute_additionals();

    inv_iterations   = std::vector<unsigned int>();
    schur_iterations = std::vector<unsigned int>();
}

parameters::parameters():
    general("General configurations"),
    vars("Problem variables"),
    file_options("All file options"),
    cmdLine_options("Allowed options")
{
    create_options();
    // Other setups.
    default_file = "default.cfg";
    if(!fexists(default_file)){
        std::cerr << "Default file must always be available!";
        write_sample_file();
        exit(1);
    }

    ifstream ifsd(default_file);
    po::store(po::parse_config_file(ifsd, file_options), vm);
    po::notify(vm);
    ifsd.close();

    setup_variables(vm);
    validate_options();
    compute_additionals();

    inv_iterations   = std::vector<unsigned int>();
    schur_iterations = std::vector<unsigned int>();
}

void parameters::create_options() {
    // Environment options setup
    general.add_options()
            ("help,h", "Produce help message")
            ("file,f", po::value<std::string>(&param_file), "Set input file")
            ("info,k", po::value<int>(&info), "Information level{0,1,2}")
            ("matrix,m", "Print matrices and matlab files")
            ("samplefile", "Write a sample configuration file");

    vars.add_options()
            ("dim,d", po::value<int>(), "Problem dimension")
            ("adv", po::value<bool>(), "Enable/Disable advection term {1|0}")
            ("div", po::value<bool>(), "Enable/Disable divergance term {1|0}")
            ("elastic,e", "Solve only elastic")
            ("precondition,c", po::value<bool>(),
             "Precondition A using whole/diagBlock of A {1|0}")
            ("inv_tol,i", po::value<double>(), "Tolerance for inverse calculation")
            ("poisson,p", po::value<double>(&POISSON), "Poisson ratio")
            ("refinement,r", po::value<int>(&refinements), "Number of refinements")
            ("schur_tol,s", po::value<double>(), "Tolerance to compute Schur complement")
            ("system_tol,t", po::value<double>(), "System solver tollerance")
            ("one_schur,o", "One schur iteration")
            ("young,y", po::value<double>(&YOUNG), "Set Young's modulus")
            ("threshold,z", po::value<double>(), "Application threashold");

    file_options.add_options()
            ("dimension",po::value<int>(&dimension), "Set Problem dimension")
            ("degree",po::value<int>(&degree), "Set degree of polynomial")
            ("refinement", po::value<int>(&refinements), "Number of refinements")
            ("young", po::value<double>(&YOUNG), "Set Young's modulus")
            ("poisson", po::value<double>(&POISSON), "Poisson ratio")
            ("eta", po::value<double>(&ETA), "ETA value")
            ("divisions.x", po::value<int>(&xdivisions), "Number of initial divisions in width")
            ("divisions.y", po::value<int>(&ydivisions), "Number of initial divisions in depth")
            ("earth.width", po::value<double>(&x2), "Set earth width")
            ("earth.depth", po::value<double>(&y1), "Set earth depth")
            ("earth.density", po::value<double>(&rho_r), "Density of earth")
            ("earth.gravity", po::value<double>(&gravity), "Gravity of earth")
            ("ice.width", po::value<double>(&Ix), "Set Ice width")
            ("ice.depth",po::value<double>(&h), "Set Ice depth")
            ("ice.density", po::value<double>(&rho_i), "Density of ice")
            ("boundaries.ice", po::value<string>(), "Under ice boundary")
            ("boundaries.up", po::value<string>(), "Top wall boundary")
            ("boundaries.right", po::value<string>(), "Right wall boundary")
            ("boundaries.left", po::value<string>(), "Left wall boundary")
            ("boundaries.bottom", po::value<string>(), "Bottom wall boundary")
            ("enable.load", po::value<bool>(&load_enabled), "Enable/Disable load")
            ("enable.weight", po::value<bool>(&weight_enabled), "Enable/Disable weight")
            ("enable.advection", po::value<bool>(&adv_enabled), "Enable/Disable advection")
            ("enable.divergance", po::value<bool>(&div_enabled), "Enable/Disable divergance")
            ("enable.precondition", po::value<bool>(&precond),
             "Precondition A using whole/diagBlock of A {1|0}")
            ("enable.one_schur_it", po::value<bool>(&one_schur_it), "One Schure iteration")
            ("tolerance.inverse",po::value<double>(&InvMatPreTOL), "Tolerance for inverse calculation")
            ("tolerance.schur", po::value<double>(&SchurTOL), "Tolerance to compute Schur complement")
            ("tolerance.system", po::value<double>(&TOL), "System solver tolerance")
            ("amg.threshold", po::value<double>(&threshold), "AMG preconditioner threshold");

    cmdLine_options.add(general).add(vars);
}

void parameters::setup_variables(po::variables_map& vm) {
    if(vm.count("help")){
        cout << cmdLine_options << "\n";
        exit(0);
    }
    if(vm.count("samplefile")){
        write_sample_file();
        exit(0);
    }
    if(vm.count("matrix")){
        print_matrices = true;
    }
    if(vm.count("dim")){
        dimension = vm["dim"].as<int>();
    }
    if(vm.count("adv")){
        adv_enabled = vm["adv"].as<bool>();
    }
    if(vm.count("div")){
        div_enabled = vm["div"].as<bool>();
    }
    if(vm.count("elastic")){
        adv_enabled = false;
        div_enabled = false;
    }
    if(vm.count("precondition")){
        precond = vm["precondition"].as<bool>();
    }
    if(vm.count("inv_tol")){
        InvMatPreTOL = vm["inv_tol"].as<double>();
    }
    if(vm.count("schur_tol")){
        SchurTOL = vm["schur_tol"].as<double>();
    }
    if(vm.count("system_tol")){
        TOL = vm["system_tol"].as<double>();
    }
    if(vm.count("one_schur"))
        one_schur_it = true;
    if(vm.count("threshold")){
        threshold = vm["threshold"].as<double>();
    }

    if(vm.count("boundaries.ice")){
        b_ice =  str2boundary(vm["boundaries.ice"].as<string>());
    }
    if(vm.count("boundaries.up")){
        b_up =  str2boundary(vm["boundaries.up"].as<string>());
    }
    if(vm.count("boundaries.right")){
        b_right =  str2boundary(vm["boundaries.right"].as<string>());
    }
    if(vm.count("boundaries.left")){
        b_left =  str2boundary(vm["boundaries.left"].as<string>());
    }
    if(vm.count("boundaries.bottom")){
        b_bottom =  str2boundary(vm["boundaries.bottom"].as<string>());
    }
}

void parameters::compute_additionals() {
    x1 = 0.0;
    y2 = 0.0;
    L  = x2;
    U  = 1.0;
    S  = YOUNG;
    T  = 1.0;

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
    double E = YOUNG*S, v = POISSON;

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
}

void parameters::validate_options(){
    using namespace std;
    bool is_correct = true;
    if(dimension > 3 || dimension < 2){
        cerr << "Problem dimension is not supported.\n";
        is_correct = false;
    }
    if(degree <= 0 || xdivisions <= 0 || ydivisions <= 0 ||
            YOUNG <= 0 || ETA <=0 || refinements <= 0 || TOL <= 0 ||
            SchurTOL <= 0 || InvMatPreTOL <= 0 || rho_i <= 0 || rho_r <= 0){
        cerr << "Invalid or negative value for one of the variables below:\n";
        cerr << "degree, xdivisions, ydivisions, young, eta, refinements, TOL, schure_tol, InvMatPreTOL, rho_r, rho_i\n";
        is_correct = false;
    }
    if(Ix > x2){
        cerr << "Ice width is either too small or biger than earth width\n";
        is_correct = false;
    }
    if(POISSON < 0 || POISSON > 0.5){
        cerr << "Poisson's ratio should be in the interval 0-0.5\n";
        is_correct = false;
    }
    if(!is_correct)
        exit(1);
}

// Helper methods
bFlags::boundary_Type parameters::str2boundary(std::string tempSt){
    bFlags::boundary_Type bt = bFlags::NEUMANN;
    if(tempSt == std::string("NEUMANN"))
        bt = bFlags::NEUMANN;
    else if(tempSt == std::string("NO_SLIP"))
        bt = bFlags::NO_SLIP;
    else if(tempSt == std::string("V_SLIP"))
        bt = bFlags::V_SLIP;
    else if(tempSt == std::string("LOAD"))
        bt = bFlags::LOAD;
    else if(tempSt == std::string("FREE"))
        bt = bFlags::FREE;

    return bt;
}

std::string parameters::boudary2str(bFlags::boundary_Type bt){
    std::string tempSt;
    switch (bt){
    case bFlags::NEUMANN:
        tempSt = "NEUMANN";
        break;
    case bFlags::NO_SLIP:
        tempSt = "NO_SLIP";
        break;
    case bFlags::V_SLIP:
        tempSt = "V_SLIP";
        break;
    case bFlags::LOAD:
        tempSt = "LOAD";
        break;
    case bFlags::FREE:
        tempSt = "FREE";
        break;
    }
    return tempSt;
}

bool parameters::fexists(std::string filename){
    struct stat buf;
    if (stat(filename.c_str(), &buf) != -1){
        return true;
    }
    return false;
}

std::vector<std::string>& parameters::split(const std::string &s, char delim, std::vector<std::string> &elems) {
    stringstream ss(s);
    string item;
    while (getline(ss, item, delim)) {
        elems.push_back(item);
    }
    return elems;
}


std::vector<std::string> parameters::split(const std::string &s, char delim) {
    vector<string> elems;
    split(s, delim, elems); return elems;
}

std::ostream & parameters::print_values(std::ostream &ostr){
    int c1=20,ch=25,cf=35;
    ios::fmtflags f(ostr.flags());
    ostr << setw(ch) << setfill('=') << "Parameters" << endl;
    ostr << left << setfill('.');

    ostr<< setw(c1) << "Dimension=" << dimension << endl;
    ostr<< setw(c1) << "Degree=" << degree << endl;
    ostr<< setw(c1) << "L=" << L << endl;
    ostr<< setw(c1) << "U=" << U << endl;
    ostr<< setw(c1) << "S=" << S << endl;
    ostr<< setw(c1) << "T=" << T << endl;
    ostr<< setw(c1) << "x1=" << x1*L << endl;
    ostr<< setw(c1) << "x2=" << x2*L << endl;
    ostr<< setw(c1) << "y1=" << y1*L << endl;
    ostr<< setw(c1) << "y2="<< y2*L << endl;
    ostr<< setw(c1) << "Ix=" << Ix*L << endl;
    ostr<< setw(c1) << "h=" << h << endl;
    ostr<< setw(c1) << "b_ice=" << boudary2str(b_ice) << endl;
    ostr<< setw(c1) << "b_up=" << boudary2str(b_up) << endl;
    ostr<< setw(c1) << "b_left=" << boudary2str(b_left) << endl;
    ostr<< setw(c1) << "b_right=" << boudary2str(b_right) << endl;
    ostr<< setw(c1) << "b_down=" << boudary2str(b_bottom) << endl;
    ostr<< setw(c1) << "refinements=" << refinements << endl;
    ostr<< setw(c1) << "xdivisions=" << xdivisions << endl;
    ostr<< setw(c1) << "ydivisions=" << ydivisions << endl;
    ostr<< setw(c1) << "YOUNG=" << YOUNG*S << endl;
    ostr<< setw(c1) << "POISSON=" << POISSON << endl;
    ostr<< setw(c1) << "ETA=" << ETA << endl;
    ostr<< setw(c1) << "rho_i=" << rho_i << endl;
    ostr<< setw(c1) << "rho_r=" << rho_r << endl;
    ostr<< setw(c1) << "gravity=" << gravity << endl;
    ostr<< setw(c1) << "load_enabled=" << load_enabled << endl;
    ostr<< setw(c1) << "weight_enabled=" << weight_enabled << endl;
    ostr<< setw(c1) << "adv_enabled=" << adv_enabled << endl;
    ostr<< setw(c1) << "div_enabled=" << div_enabled << endl;
    ostr<< setw(c1) << "precond=" << precond << endl;
    ostr<< setw(c1) << "InvMatPreTOL=" << InvMatPreTOL << endl;
    ostr<< setw(c1) << "SchurTOL=" << SchurTOL << endl;
    ostr<< setw(c1) << "TOL=" << TOL << endl;
    ostr<< setw(c1) << "threshold=" << threshold << endl;
    ostr<< setw(c1) << "one_schur_it=" << one_schur_it << endl;
    ostr<< setw(c1) << "info=" << info << endl;
    ostr<< setw(c1) << "print_matrices=" << print_matrices << endl;
    ostr<< setw(c1) << "system_iter=" << system_iter << endl;
    ostr<< setw(c1) << "load=" << load << endl;
    ostr<< setw(c1) << "weight " << weight << endl;
    ostr<< setw(c1) << "g0=" << g0 << endl;
    ostr<< setw(c1) << "scale1=" << scale1 << endl;
    ostr<< setw(c1) << "scale2=" << scale2 << endl;
    ostr<< setw(c1) << "scale3=" << scale3 << endl;
    ostr<< setw(c1) << "alpha=" << alpha  << endl;
    ostr<< setw(c1) << "beta=" << beta  << endl;
    ostr<< setw(c1) << "delta=" << delta << endl;
    ostr<< setw(c1) << "gamma" << gamma << endl;
    ostr<< setw(c1) << "str_poisson=" << str_poisson <<endl;

    ostr << setw(cf) << setfill('=') << "" << endl;

    ostr.flags(f);
    return ostr;
}

std::ostream & parameters::print_variables(std::ostream & outStr){
    using namespace std;
    ios::fmtflags f(outStr.flags());
    int c1 = (one_schur_it) ? 15 : 22;

    outStr << left;
    outStr << "Problem:" << endl;
    outStr << setw(c1) << "\tP_00: ";
    (precond) ? outStr << "A_00" << endl : outStr << "diag(A_00)" << endl;

    outStr << setw(c1) << "\tRefinements: " << refinements << endl;


    if(one_schur_it)
        outStr << setw(c1) << "\tSystemTol: " << TOL << endl;
    else
        outStr << setw(c1) << "\tTOL(P_00,Schur): "
               << TOL << "(" << InvMatPreTOL << ", " << SchurTOL << ")"
           << endl;

    outStr << setw(c1) << "\tAdv/Div: ";
    (adv_enabled) ? outStr << "enabled, " : outStr << RED << "disabled, " << RESET;
    (div_enabled) ? outStr << "enabled" << endl : outStr << RED << "disabled" << RESET<< endl;

    outStr << setw(c1) << "\tLoad/Weight: ";
    (load_enabled) ? outStr << "enabled, " : outStr << RED << "disabled, " << RESET;
    (weight_enabled) ? outStr << "enabled" << endl : outStr << RED << "disabled" << RESET << endl;

    outStr << setw(c1) << "\tPoisson: " << POISSON << endl;
    outStr << setw(c1) << "\tYoung: " << YOUNG*S << endl;
    outStr << setw(c1) << "\tEta: " << ETA << endl;
    outStr << setw(c1) << "\tEarth (Km)" << x2*L*1e-3 << ", " << y1*L*1e-3 << endl;
    outStr << setw(c1) << "\tIce (Km)" << Ix*L*1e-3 << ", " << h*1e-3 << endl << endl;

    outStr.flags(f);
    return outStr;
}

void parameters::write_sample_file(){
    using namespace std;
    // create and open an archive for input
    ofstream ofs;
    string filename = "default.cfg";
    ofs.open(filename.c_str(), ios::trunc);

    if(!ofs.is_open()){
        cerr << "Could not open sample file\n";
        exit(1);
    }
    ofs << "## Problem dimension.\n" <<
           "dimension=2\n" <<
           "## Degree of the polynomial basis functions.\n"
           "degree=1\n" <<
           "## Number of refinement\n" <<
           "refinement=3\n" <<
           "## Young's modulus.\n" <<
           "young=4e11\n" <<
           "## POISSON ratio.\n" <<
           "poisson=0.2\n" <<
           "## Description.\n" <<
           "eta=100\n" <<
           "## Divisions in x and y\n"<<
           "[divisions]\n" <<
           "\tx=10\n" <<
           "\ty=4\n" <<
           "## Earth properties.\n" <<
           "[earth]\n" <<
           "\twidth=1.0e7\n" <<
           "\tdepth=-4e6\n" <<
           "\tdensity=3300\n" <<
           "\tgravity=-9.8\n" <<
           "## Ice properties.\n" <<
           "[ice]\n" <<
           "\twidth=1e6\n" <<
           "\tdepth=2e3\n" <<
           "\tdensity=917\n" <<
           "## Boundary types.\n" <<
           "[boundaries]\n" <<
           "\tice=LOAD\n" <<
           "\tup=FREE\n" <<
           "\tleft=V_SLIP\n" <<
           "\tright=NEUMANN\n" <<
           "\tbottom=NO_SLIP\n" <<
           "## Toggles to enable/disable\n" <<
           "[enable]\n"<<
           "\tload=1\n"
           "\tweight=0\n" <<
           "\tadvection=1\n" <<
           "\tdivergance=1\n" <<
           "\tone_schur_it=0\n" <<
           "## Use the whole block A_11 in the preconditioner.\n" <<
           "\tprecondition=1\n" <<
           "## Solvers tolerances\n"<<
           "[tolerance]\n"
           "\tinverse=1e-2\n" <<
           "\tschur=1e-1\n" <<
           "\tsystem=1e-7\n" <<
           "# AMG options\n" <<
           "[amg]\n"
           "\tthreshold=0.02\n";

    ofs.close();
}

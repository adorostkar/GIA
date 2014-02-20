
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
    cmdLine_options("Allowed options")
{
    // Other setups.
    default_file = "default.config";
    if(!fexists(default_file)){
        std::cerr << "Default file must always be available!\nSample file is created.\n";
        write_sample_file();
        exit(1);
    }

    setup_variables_from_file(default_file);

    // Environment options setup
    set_env_vars();

    po::variables_map vm;
    po::store(po::parse_command_line(argc, argv, cmdLine_options), vm);
    po::notify(vm);

    if(vm.count("file"))
        param_file = vm["file"].as<std::string>();

    if(fexists(param_file)){
        setup_variables_from_file(param_file);
    }

    setup_variables_from_cmd(vm);

    validate_options();

    compute_additionals();

    inv_iterations   = std::vector<unsigned int>();
    schur_iterations = std::vector<unsigned int>();
}

parameters::parameters():
    general("General configurations"),
    vars("Problem variables"),
    cmdLine_options("Allowed options")
{
    // Other setups.
    default_file = "default.config";
    if(!fexists(default_file)){
        std::cerr << "Default file must always be available!";
        write_sample_file();
        exit(1);
    }

    setup_variables_from_file(default_file);

    validate_options();

    compute_additionals();

    inv_iterations   = std::vector<unsigned int>();
    schur_iterations = std::vector<unsigned int>();
}

void parameters::set_env_vars() {
    // Environment options setup
    general.add_options()
            ("help,h", "Produce help message")
            ("file,f", po::value<std::string>(&param_file), "Set input file")
            ("info,k", po::value<int>(&info), "Information level{0,1,2}")
            ("matrix,m", "Print matrices and matlab files")
            ("samplefile", "Write a sample configuration file");

    vars.add_options()
            ("dim", po::value<int>(&dimension), "Problem dimension")
            ("adv,a", po::value<bool>(&adv_enabled), "Enable/Disable advection term {1|0}")
            ("div,d", po::value<bool>(&div_enabled), "Enable/Disable divergance term {1|0}")
            ("elastic,e", "Solve only elastic")
            ("precondition,c", po::value<bool>(&precond),
             "Precondition A using whole/diagBlock of A {1|0}")
            ("inv_tol,i", po::value<double>(&InvMatPreTOL), "Tolerance for inverse calculation")
            ("poisson,p", po::value<double>(&POISSON), "Poisson ratio")
            ("refinement,r", po::value<int>(&refinements), "Number of refinements")
            ("schur_tol,s", po::value<double>(&SchurTOL), "Tolerance to compute Schur complement")
            ("system_tol,t", po::value<double>(&TOL), "System solver tollerance")
            ("one_schur,o", "One schur iteration")
            ("weight,w", po::value<bool>(&weight_enabled), "Enable/disable weight computation{1|0}")
            ("young,y", po::value<double>(&YOUNG), "Set Young's modulus")
            ("threshold,z", po::value<double>(&threshold), "Application threashold");

    cmdLine_options.add(general).add(vars);
}

void parameters::setup_variables_from_cmd(po::variables_map& vm) {
    if(vm.count("help")){
        cout << cmdLine_options << "\n";
        exit(0);
    }
    if(vm.count("samplefile")){
        write_sample_file();
        exit(0);
    }
    if(vm.count("matrix")){
        print_local = true;
        print_matrices = true;
    }
    if(vm.count("elastic")){
        adv_enabled = false;
        div_enabled = false;
    }
    if(vm.count("one_schur"))
        one_schur_it = true;
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

void parameters::setup_variables_from_file(std::string filename) {
    using namespace std;
    // create and open an archive for input
    ifstream ifs(filename.c_str());
    // Read parameters from file
    string token;
    vector<string> tokens;
    if(!ifs.is_open()){
        std::cerr << "Could not open " << filename << ", aborting execution" << std::endl;
    }

    while(!ifs.eof() && ifs.good()){
        ifs >> token;
        tokens = split(token, '=');
        set_values(tokens);
    }
    ifs.close();
}

void parameters::set_values(std::vector<std::string>& tokens){
    using namespace std;
    if(tokens.size() != 2)
        return;
    string token = tokens[0], value = tokens[1];

    if(token == "dimension"){
        dimension = atoi(value.c_str());
    }else if(token == "degree"){
        degree = atoi(value.c_str());
    }else if(token == "earth_width"){
        x2 = atof(value.c_str());
    }else if(token == "earth_depth"){
        y1 = atof(value.c_str());
    }else if(token == "ice_width"){
        Ix = atof(value.c_str());
    }else if(token == "ice_depth"){
        h = atof(value.c_str());
    }else if(token == "b_ice"){
        b_ice = str2boundary(value);
    }else if(token == "b_up"){
        b_up = str2boundary(value);
    }else if(token == "b_left"){
        b_left = str2boundary(value);
    }else if(token == "b_right"){
        b_right = str2boundary(value);
    }else if(token == "b_bottom"){
        b_bottom = str2boundary(value);
    }else if(token == "refinements"){
        refinements = atoi(value.c_str());
    }else if(token == "xdivisions"){
        xdivisions = atoi(value.c_str());
    }else if(token == "ydivisions"){
        ydivisions = atoi(value.c_str());
    }else if(token == "young"){
        YOUNG = atof(value.c_str());
    }else if(token == "poisson"){
        POISSON = atof(value.c_str());
    }else if(token == "eta"){
        ETA = atof(value.c_str());
    }else if(token == "rho_i"){
        rho_i = atof(value.c_str());
    }else if(token == "rho_r"){
        rho_r = atof(value.c_str());
    }else if(token == "gravity"){
        gravity = atof(value.c_str());
    }else if(token == "load_enabled"){
        load_enabled = (value == "true" || value == "1" || value == "on");
    }else if(token == "weight_enabled"){
        weight_enabled = (value == "true" || value == "1" || value == "on");
    }else if(token == "adv_enabled"){
        adv_enabled = (value == "true" || value == "1" || value == "on");
    }else if(token == "div_enabled"){
        div_enabled = (value == "true" || value == "1" || value == "on");
    }else if(token == "precondition"){
        precond = (value == "true" || value == "1" || value == "on");
    }else if(token == "InvMatPreTOL"){
        InvMatPreTOL = atof(value.c_str());
    }else if(token == "SchurTOL"){
        SchurTOL = atof(value.c_str());
    }else if(token == "TOL"){
        TOL = atof(value.c_str());
    }else if(token == "threshold"){
        threshold = atof(value.c_str());
    }else if(token == "one_schur_it"){
        one_schur_it = (value == "true" || value == "1" || value == "on");
    }
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
    ostr<< setw(c1) << "print_local=" << print_local << endl;
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
    int c1=25, c2=20, c3=25;

    ostringstream str, wl_e, ad_e;
    wl_e << "(" << weight_enabled << ", " << load_enabled << ")";
    ad_e << "(" << adv_enabled <<", "<< div_enabled <<")";
    str << TOL << "(" << InvMatPreTOL << ",";
    if(one_schur_it)
        str << "x" << ")";
    else
        str<< SchurTOL << ")";
    string title = (x2 == Ix) ? "Case: Footing" : "Case: Uniform load";
    (precond == 0) ? title.append(", P_00: diag(A)") : title.append(", P_00: A");

    outStr << left;
    outStr << title << endl;
    outStr << "\t"
           << setw(c1) << "YOUNG = " << setw(c2) << YOUNG*S
           << setw(c3) << "POISSON = " << POISSON
           << endl;

    outStr << "\t"
           << setw(c1) << "scale1(L^2/(SU)) = " << setw(c2) << scale1
           << setw(c3) << "ETA = " << ETA
           << endl;

    outStr << "\t"
           << setw(c1) << "scale3 (L/S)*rho_r*g = " << setw(c2) << scale3
           << setw(c3) << "scale2(L/SU) = " << scale2
           << endl;

    outStr << "\t"
           << setw(c1) << "weight/load enabled=" << setw(c2) << wl_e.str()
           << setw(c3) << "adv/div enabled" << ad_e.str()
           << endl;

    outStr << "\t"
           << setw(c1) << "refinements" << setw(c2) << refinements
           << setw(c3) << "AMG threshold" << threshold
           << endl;

    outStr << "\t"
           << setw(c1) << "TOL(Inner,Schur)" << setw(c2) << str.str()
           << endl;

    outStr << "\t" << "Earth and ice (width, depth) = (" << x2 << ", " << y1 << "), (" << Ix << ", " << h/L << ")" << endl;

    outStr.flags(f);
    return outStr;
}

void parameters::write_sample_file(){
    using namespace std;
    // create and open an archive for input
    ofstream ofs;
    string filename = "sample.config";
    ofs.open(filename.c_str(), ios::trunc);

    if(!ofs.is_open()){
        cerr << "Could not open sample file\n";
        exit(1);
    }
    ofs << "## Problem dimension.\n" <<
           "dimension=2\n" <<
           "## Degree of the polynomial basis functions.\n"
           "degree=1\n" <<
           "## Earth width.\n" <<
           "earth_width=1.0e7\n" <<
           "## Earth depth\n" <<
           "earth_depth=-4e6\n" <<
           "## Load(Ice) width.\n" <<
           "ice_width=1e6\n" <<
           "## Load depth.\n" <<
           "ice_depth=2e3\n" <<
           "## Boundary type under the ice.\n" <<
           "b_ice=LOAD\n" <<
           "## Boundary type on the surface.\n" <<
           "b_up=FREE\n" <<
           "## Boundary type on the left wall.\n" <<
           "b_left=V_SLIP\n" <<
           "## Boundary type on the write wall.\n" <<
           "b_right=NEUMANN\n" <<
           "## Boundary type on the bottom wall.\n" <<
           "b_bottom=NO_SLIP\n" <<
           "## Number of refinement\n" <<
           "refinements=3\n" <<
           "## Number of initial divisions in x direction.\n" <<
           "xdivisions=10\n" <<
           "## Number of initial divisions in y direction.\n" <<
           "ydivisions=4\n" <<
           "## Young's modulus.\n" <<
           "young=4e11\n" <<
           "## POISSON ratio.\n" <<
           "poisson=0.2\n" <<
           "## Description.\n" <<
           "eta=100\n" <<
           "## Ice density\n" <<
           "rho_i=917\n" <<
           "## Earth density.\n" <<
           "rho_r=3300\n" <<
           "## Gravity\n" <<
           "gravity=-9.8\n" <<
           "## if load is incorporated in the formulation or is it zero.{0,1}\n" <<
           "load_enabled=1\n"
           "## if weight is incorporated in the formulation or is it zero.{0,1}\n" <<
           "weight_enabled=0\n" <<
           "## If advection term is enabled.{0,1}\n" <<
           "adv_enabled=1\n" <<
           "## If divergance term is enabled.{0,1}\n" <<
           "div_enabled=1\n" <<
           "## {0 = diag(A), 1 = A} Whether to rewrite block A_11 in the preconditioner.\n" <<
           "precondition=0\n" <<
           "## The tolerances of the inner solvers for computing the effect of inv(A).\n" <<
           "InvMatPreTOL=1e-2\n" <<
           "## The tolerances of the inner solvers for computing the effect of \tilde(A).\n" <<
           "SchurTOL=1e-1\n" <<
           "## The system tolerance controlling the relative L2-error in the global scheme.\n" <<
           "TOL=1e-7\n" <<
           "## AMG preconditioner threshold.\n" <<
           "threshold=0.02\n" <<
           "## {0,1} Whether to do only one Schure iteration on the matrix.\n" <<
           "one_schur_it=0\n";

    ofs.close();
}

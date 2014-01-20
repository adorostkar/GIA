/* TODO:
 - Add argument check. i.e boundaries for the arguments that get int or double
 - Find a way to add variables added to struct long_options automatically to
    print_ methods and help method.
    + Idea: create a discription variable and iterate through that.
 */

// include headers that implement a archive in simple text format
#include <sys/stat.h>

#include <list>
#include <iostream>
#include <cmath>
#include <vector>
#include <getopt.h>
#include <cstdlib>
#include <cstdio>
#include <iomanip>
#include <fstream>
#include <sstream>
#include <string>
#include <string.h>
#include <algorithm>

#ifndef PARAMETERS_H
#define PARAMETERS_H
// short options
const char* const short_options = "ha:d:p:y:z:i:s:t:r:f:k:g:ewumv";
// Long options
const struct option long_options[] = { // An array describing valid long options.
	{ "footing",		no_argument,		0,		'0' },
	{ "large",			no_argument,		0,		'1' },
	{ "uniform",		no_argument,		0,		'2' },
	{ "file",			no_argument,		0,		'3' },
	{ "writeback",		no_argument,		0,		'4' },
    { "dim",    		required_argument,	0,		'5' },
	{ "adv",			required_argument,	0,		'a' },
	{ "div",			required_argument,	0,		'd' },
	{ "elastic",		no_argument,		0,		'e' },
	{ "surf_samples",	required_argument,	0,		'f' },
	{ "precond",		required_argument,	0,		'g' },
	{ "help",			no_argument,		0,		'h' },
	{ "inv_tol",		required_argument,	0,		'i' },
	{ "info",			required_argument,	0,		'k' },
	{ "matlab_print",	no_argument,		0,		'm' },
	{ "poisson",		required_argument,	0,		'p' },
	{ "refinements",	required_argument,	NULL,	'r' },
	{ "schur_tol",		required_argument,	0,		's' },
	{ "system_tol",		required_argument,	0,		't' },
	{ "unique_schur",	no_argument,		0,		'u' },
	{ "verbose",		0,					NULL,	'v' },
	{ "weight",			no_argument,		0,		'w' },
	{ "young",			required_argument,	0,		'y' },
	{ "threshold",		required_argument,	0,		'z' },
	{ NULL,				0,					NULL,	0   }   /* Required at end of array.  */
	};

// Boundary types enumerator,
// Different available boundaries.
enum boundary_Type {
	NEUMANN = 0,
	NO_SLIP = 1,
    V_SLIP  = 2,
    LOAD    = 3,
    FREE    = 4
	};

// Parameter class,
// This class contains all application options including the ones only seen by
// application and the ones available to the user to change.
// User can alter options by command line or by using an option file.
class parameters
{
public:
    ~parameters();
    static parameters* getInstance();
    static parameters* getInstance(int _argc, char *_argv[]);

	// This method reads the parameters
	// from the file(if no value is given, defaults will be used)
	void read_Parameters();
	// Write back parameters to file
	void write_Parameters();
	// This method reads the command line options
	void parse_command_line(int argc, char** argv);

	// findout the filename fo parameter file
	void find_filename(int argc, char** argv);
	// Print argument usage
	void print_usage(int exitNum);

	// Print parameter values
	void print_values();

	// Print variable settings, usefull information for the problem
	void print_variables();

	// Set parameter values to default
    void set_default();

    // get case number
    int cases(){return cas;}
	// File name of the option file
	std::string					paramFile;
	/* problem variables */

	int							dimension, degree;
	int							refinements;
	int							xdivisions, ydivisions;
	// info = {0,1,2}
	int							solver,precond, info;
	int							system_iter,surf_samples;

	double						adv, div;
	double						load, weight;
	double						gravity;
	double						InvMatPreTOL,SchurTOL,TOL,threshold;
	double						YOUNG,POISSON,ETA;
	double						rho_i,rho_r,g0;
	double						x1,y1,x2,y2,Ix,h;
	double						L,U,S,T; // Scaling parameters for length, displacement, Stress and time
	double						scale1,scale2,scale3;
	double						alpha, beta, delta, gamma;
    // alpha = rho_i*h/rho_r
    // beta  = E*(1-2*v)/( 4*v*(1+v) )
    // delta = (1+v)(1-2v)/(E(1-v)) = 1/(2*mu + lambda)
    // gamma = 2v(1+v)/E(1-v) = 1/(2*beta + mu)

	bool						print_local,print_matrices, one_schur_it;
	bool						load_enabled,weight_enabled,adv_enabled,div_enabled;
	bool						solve;
    // writeback parameters to file.
	bool						writeback;
	
    boundary_Type				b_ice, b_up, b_left, b_right, b_bottom;
    std::vector<unsigned int>		inv_iterations, schur_iterations;
	
	std::string					str_poisson;
	std::stringstream			dofs;
private:
	bool						verbose; // Verbose parameter set
    // Problem case number. This can only be retrieved and can't be changed
    // during runtime.
    int                         cas;

    // Instance flag
    static bool instanceFlag;
    // Instance
    static parameters *singlton;

    //default constructor
    parameters(const std::string f = "vardata.conf");
    // If the application does not have any command line options
    // this works the same as default constructor
    // If no --file options is passed, the program assumes default file is used.
    parameters(int argc, char* argv[], const std::string f = "vardata.conf");

	// Used to convert string boundary type to enum boundary types
	boundary_Type str2boundary(std::string tempSt);
	// convert boundaries to text
	std::string boudary2str(boundary_Type bt);
    
    bool fexists();

    // Compute additional parameters i.e. alpha, beta, gamma, delta, g0
    // scale1, scale2, weight, load, adv_enabled, div_enabled, scale2
    // inv_iteration, schure_iterations
    void compute_additionals();
};
#endif

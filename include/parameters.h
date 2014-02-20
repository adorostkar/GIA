/** TODO
 */
#include <boost/program_options.hpp>
#include <cstdlib>
#include <cstdio>
#include <fstream>
#include <iostream>
#include <iomanip>
#include <vector>
#include <string>
#include <string.h>
#include <sstream>
#include <sys/stat.h>

#ifndef PARAMETERS_H_
#define PARAMETERS_H_

/*! Boundary types enumerator,
 * Different available boundaries.
 * Used power of two to be able to
 * use multiple flags.
 */
struct bFlags{
    enum boundary_Type {
        NEUMANN = 1<<0,
        NO_SLIP = 1<<1,
        V_SLIP  = 1<<2,
        LOAD    = 1<<3,
        FREE    = 1<<4
    };
};

class parameters {
public:
    // Variables
    std::string					param_file, default_file;

    int							dimension, degree,
                                refinements,
                                xdivisions, ydivisions,
                                info, // {0,1,2}
                                system_iter;

    double						load, weight,
                                gravity,
                                InvMatPreTOL, SchurTOL, TOL, threshold,
                                YOUNG, POISSON, ETA,
                                rho_i, rho_r, g0,
                                x1, x2, y1, y2, Ix, h,
                                L, U, S, T;

    /*!
     * \brief scale1 = L^2/(SU)
     * \brief scale2 = L/SU
     * \brief scale3 = (L/S)*rho_r*g
     */
    double						scale1, scale2, scale3;

    /*!
     * \brief alpha = rho_i*h/rho_r
     * \brief beta  = E*(1-2*v)/( 4*v*(1+v) )
     * \brief delta = (1+v)(1-2v)/(E(1-v)) = 1/(2*mu + lambda)
     * \brief gamma = 2v(1+v)/E(1-v) = 1/(2*beta + mu)
     */
    double						alpha, beta, delta, gamma;

    /*!
     * \brief print_local prints local matrices to output
     * \brief print_matrices prints global matrices to output
     * \brief one_schur_it uses one iteration to compute Schure
     */
    bool						precond, print_local,print_matrices, one_schur_it;

    /*!
     * \brief load_enabled is load enabled on the surface.
     * \brief weight_enabled is body force is enabled.
     * \brief adv_enabled is advection term is enabled.
     * \brief div_enabled id divergance term is enabled.
     */
    bool						load_enabled,weight_enabled,adv_enabled,div_enabled;

    bFlags::boundary_Type		b_ice, b_up, b_left, b_right, b_bottom;
    std::vector<unsigned int>	inv_iterations, schur_iterations;

    std::string					str_poisson;
    std::stringstream			dofs;


    // Methods
    ~parameters();
    static parameters* getInstance();
    static parameters* getInstance(int _argc, char *_argv[]);

    void write_sample_file();
    std::ostream &print_variables(std::ostream & str);
    std::ostream &print_values(std::ostream &ostr);
private:
    // Variables
    boost::program_options::variables_map vm;
    boost::program_options::options_description general;
    boost::program_options::options_description vars;
    boost::program_options::options_description file_options;
    boost::program_options::options_description cmdLine_options;
    // Instance flag
    static bool instanceFlag;
    // Instance
    static parameters *singlton;

    // Methods
    parameters();
    parameters(int argc, char* argv[]);
    void create_options();
    void compute_additionals();
    void setup_variables(boost::program_options::variables_map & vm);

    bFlags::boundary_Type str2boundary(std::string tempSt);
    // convert boundaries to text
    std::string boudary2str(bFlags::boundary_Type bt);

    bool fexists(std::string filename);
    std::vector<std::string>& split(const std::string &s, char delim, std::vector<std::string> &elems);
    std::vector<std::string> split(const std::string &s, char delim);
    void validate_options();
};

#endif /* PARAMETERS_H_ */

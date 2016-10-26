// ESW Fri Oct 21 21:37:18 EDT 2016
// Structs and functions that manage all input parameters.

#ifndef MG_INPUTS
#define MG_INPUTS

#include <vector>

#include "operators.h"
#include "null_gen.h"
#include "mg.h"
#include "mg_complex.h"
#include "u1_utils.h"

// Print usage. 
void display_usage();

// Forward declaration.
struct mg_input_params;

// Parse inputs, put into mg_input_params.
int parse_inputs(int argc, char** argv, mg_input_params *params);

struct mg_input_params
{
    
    double mass;
    op_type opt;
    int lattice_size_x;
    int lattice_size_y; 
    
    // Information about the outer solver.
    outer_solver out_solve; 
    double outer_precision; 
    int outer_max_iter; 
    bool outer_restart;
    int outer_restart_freq; 
    
    // Information about the null space generation.
    null_vector_params nvec_params; 
    bool nvec_use_eigen; 
    
    // Information about MG.
    std::vector<int> blocksizes; 
    int n_refine; 
    mg_multilevel_type mlevel_type; 
    bool normal_eqn_mg; 
    
    // Information about MG smoother.
    inner_solver in_smooth; 
    double omega_smooth;
    std::vector<int> pre_smooths;
    std::vector<int> post_smooths; 
    bool normal_eqn_smooth; 
    
    // Information about the coarse solve. 
    inner_solver in_solve; 
    std::vector<double> inner_precisions; 
    int inner_restart; 
    int inner_max; 
    
    // Information about gauge field.
    gauge_create_type gauge_load; 
    bool do_gauge_transform; 
    double beta; 
    char* load_cfg; 
    bool do_load; 
    
    // Test types
    bool do_free; 
//#ifdef EIGEN_TEST
    bool do_eigentest;
    int set_eigen;
    int set_cv;
//#endif

    // Constructor with defaults. 
    mg_input_params();
    
};

#endif


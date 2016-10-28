// This is a sample code that constructs general kinetic terms!

#include <iostream>
#include <iomanip> // to set output precision.
#include <fstream>
#include <cmath>
#include <string>
#include <sstream>
#include <complex>
#include <random>
#include <vector>
#include <cstring> // should be replaced by using sstream


#include "generic_vector.h"
#include "verbosity.h"

#include "generic_inverters.h"
#include "generic_inverters_precond.h"

#include "generic_bicgstab.h"
#include "generic_cg.h"
#include "generic_cr.h"
#include "generic_gcr.h"
#include "generic_minres.h"
#include "generic_cg_flex_precond.h"
#include "generic_gcr_var_precond.h"
#include "generic_bicgstab_precond.h"
#include "generic_bicgstab_l.h"

#include "mg.h"
#include "mg_complex.h"
#include "null_gen.h"
#include "u1_utils.h"
#include "lattice.h"
#include "operators.h"
#include "operators_stencil.h"
#include "coarse_stencil.h"
#include "tests.h"
#include "input_params.h"

// Are we checking eigenvalues?
#define EIGEN_TEST

#ifdef EIGEN_TEST
#include "arpack_interface.h"
#endif

// Stencil constructing golden test! Compare constructed to P^\dag A P.
// Currently doesn't work for lattices that are too small (2x2 for D, 4x4 for D^\dag D)
//#define STENCIL_CONSTRUCT_TEST

// Stencil piece test! Compare using the full stencil to 
// summing over all pieces (+x, +y, etc) of the stencil.
//#define STENCIL_PIECE_TEST

// Efficient stencil test! Test building a coarse stencil from a fine stencil.
//#define STENCIL_EFFICIENT_TEST

// Dagger stencil test! Make sure we correctly build the dagger of the stencil.
//#define STENCIL_DAGGER_TEST

// Test stencils with the shift built into the clover term vs. using the shift functionality in the stencil.
//#define STENCIL_SHIFT_TEST

// Do restrict/prolong test?
//#define PDAGP_TEST

// Do two vector restrict/prolong?
//#define PDAGP_2TEST

// Try solving just the coarse solver. 
//#define COARSE_ONLY

using namespace std; 

// Define pi.
#define PI 3.141592653589793

// Print null vectors?
//#define PRINT_NULL_VECTOR


// What type of test should we do?
enum mg_test_types
{
    TOP_LEVEL_ONLY = 0,       // only test the top level solver.
    SMOOTHER_ONLY = 1,        // Top level + smoother
    TWO_LEVEL = 2,            // Two level MG 
    THREE_LEVEL = 3           // Three level MG
};


// Custom routine to load gauge field.
void internal_load_gauge_u1(complex<double>* lattice, int x_fine, int y_fine, double BETA);

enum src_type
{
    POINT = 0,
    RANDOM_GAUSSIAN = 1,
    ORIGIN_POINT = 2 // for correlator test. 
};

int main(int argc, char** argv)
{  
    // Declare some variables.
    cout << setiosflags(ios::scientific) << setprecision(6);
    int i, j, k;
    complex<double> *lattice; // Holds the gauge field.
    complex<double> *lhs, *rhs, *check; // For some Kinetic terms.
    complex<double> *tmp, *tmp2; // For temporary space. 
    double tmp_mass; // holds temporary mass values.
    double explicit_resid = 0.0;
    double bnorm = 0.0;
    std::mt19937 generator (1337u); // RNG, 1337u is the seed. 
    inversion_info invif;
    staggered_u1_op stagif;
    
    // Introducing coordinate functions slowly.
    int nd = 2;
    int* coord = new int[nd];
    for (i = 0; i < nd; i++)
    {
        coord[i] = 0;
    }
    int lattice_size[nd];
    
    // Set parameters. This input parameters class automatically
    // gets populated by a set of default parameters. Look in 
    // "input_params.cpp" to see the defaults.
    mg_input_params params; // multigrid params. 
    
    
    // Other initialization-type deals. 
    // Describe the source type.
    src_type source = RANDOM_GAUSSIAN; // POINT, RANDOM_GAUSSIAN, or ORIGIN_POINT (for correlator test).
    
    // What test are we performing?
    mg_test_types my_test = TWO_LEVEL; //THREE_LEVEL; // TWO_LEVEL is the default which won't override anything.
    
    if (my_test == THREE_LEVEL) // FOR TEST ONLY
    {
        params.n_refine = 2;
    }
    
    if (my_test == SMOOTHER_ONLY)
    {
        params.in_solve = NONE; 
    }
    
    
    
    /////////////////////////////////////////////
    // Get a few parameters from command line. //
    /////////////////////////////////////////////
    if (!parse_inputs(argc, argv, &params)) // returns 0 on failure.
    {
        return 0;
    }
    
    //printf("Mass %.8e Blocksize %d %d Null Vectors %d\n", MASS, X_BLOCKSIZE, Y_BLOCKSIZE, n_null_vector);
    //return 0;
    
    ///////////////////////////////////////
    // End of human-readable parameters! //
    ///////////////////////////////////////
    
    
    string op_name;
    void (*op)(complex<double>*, complex<double>*, void*) = square_staggered_u1;
    void (*op_dagger)(complex<double>*, complex<double>*, void*) = square_staggered_dagger_u1; 
    switch (params.opt)
    {
        case STAGGERED:
            op = square_staggered_u1;
            op_dagger = square_staggered_dagger_u1;
            op_name = "Staggered U(1)";
            break;
        case LAPLACE:
            op_name = "Free Laplace";
            op = op_dagger = square_laplace;
            break;
        case LAPLACE_NC2:
            op_name = "Free Laplace Nc = 2";
            op = op_dagger = square_laplace;
            break;
        case G5_STAGGERED:
            op_name = "Gamma_5 Staggered U(1)";
            op = op_dagger = square_staggered_gamma5_u1;
            break;
        case STAGGERED_NORMAL:
            op_name = "Staggered U(1) Normal";
            op = op_dagger = square_staggered_normal_u1;
            break;
        case STAGGERED_INDEX:
            op_name = "Staggered U(1) Index Operator";
            op = op_dagger = staggered_index_operator;
            break;
    }
    cout << "[OP]: Operator " << op_name << " Mass " << params.mass << "\n";
    
    // Only relevant for free laplace test.
    int Nc = 1;  // Only value that matters for staggered
    if (params.opt == LAPLACE_NC2)
    {
        Nc = 2;
    }
    
    // Unset eo for Laplace. 
    if (params.opt == LAPLACE || params.opt == LAPLACE_NC2) // FOR TEST ONLY
    {
        params.nvec_params.bstrat = BLOCK_NONE;
    }
    
    // Describe the fine lattice. 
    lattice_size[0] = params.lattice_size_x;
    lattice_size[1] = params. lattice_size_y;
    
    // Create a lattice object.
    Lattice Lat(nd, lattice_size, Nc);
    
    cout << "[VOL]: X " << lattice_size[0] << " Y " << lattice_size[1] << " Volume " << Lat.get_volume();
    if (params.opt == LAPLACE || params.opt == LAPLACE_NC2) // FOR TEST ONLY
    {
        cout << " Nc " << Nc;
    }
    cout << "\n";
    
    // Do some allocation.
    // Initialize the lattice. Indexing: index = y*N + x.
    lattice = new complex<double>[2*Lat.get_lattice_size()];
    lhs = new complex<double>[Lat.get_lattice_size()];
    rhs = new complex<double>[Lat.get_lattice_size()];   
    check = new complex<double>[Lat.get_lattice_size()]; 
    tmp = new complex<double>[Lat.get_lattice_size()];
    tmp2 = new complex<double>[Lat.get_lattice_size()];
    
    // In case we need a gauge transform.
    complex<double>* gauge_trans = new complex<double>[Lat.get_lattice_size()];
    
    // Zero it out.
    zero<double>(lattice, 2*Lat.get_lattice_size());
    zero<double>(rhs, Lat.get_lattice_size());
    zero<double>(lhs, Lat.get_lattice_size());
    zero<double>(check, Lat.get_lattice_size());
    zero<double>(gauge_trans, Lat.get_lattice_size());
    zero<double>(tmp, Lat.get_lattice_size());
    zero<double>(tmp2, Lat.get_lattice_size());
    //
    
    // Fill stagif.
    stagif.lattice = lattice;
    stagif.mass = params.mass; 
    stagif.x_fine = Lat.get_lattice_dimension(0); // DEPRECIATE
    stagif.y_fine = Lat.get_lattice_dimension(1); // DEPRECIATE
    //stagif.Lat = &Lat; 
    stagif.Nc = Nc; // Only relevant for laplace test only.
    
    // Create the verbosity structure.
    inversion_verbose_struct verb;
    verb.verbosity = VERB_DETAIL;
    verb.verb_prefix = "[L1]: ";
    verb.precond_verbosity = VERB_NONE; //VERB_DETAIL;
    verb.precond_verb_prefix = "Prec ";
    
    // Describe the gauge field. 
    cout << "[GAUGE]: Creating a gauge field.\n" << flush; 
    
    switch (params.gauge_load)
    {
        case GAUGE_UNIT:
            unit_gauge_u1(lattice, Lat.get_lattice_dimension(0), Lat.get_lattice_dimension(1));
            break;
        case GAUGE_RANDOM:
            gauss_gauge_u1(lattice, Lat.get_lattice_dimension(0), Lat.get_lattice_dimension(1), generator, params.beta);
            cout << "[GAUGE]: Created a U(1) gauge field with angle standard deviation " << 1.0/sqrt(params.beta) << "\n";
            break;
        case GAUGE_LOAD:
            // Unit first in case loading fails.
            unit_gauge_u1(lattice, Lat.get_lattice_dimension(0), Lat.get_lattice_dimension(1));
            // Load the gauge field.
            if (params.do_load)
            {
                read_gauge_u1(lattice, Lat.get_lattice_dimension(0), Lat.get_lattice_dimension(1), params.load_cfg);
                cout << "[GAUGE]: Loaded a U(1) gauge field from " << params.load_cfg << "\n";
            }
            else // various predefined cfgs. 
            {
                internal_load_gauge_u1(lattice, Lat.get_lattice_dimension(0), Lat.get_lattice_dimension(1), params.beta); // defined at end of file.
            }
            break;
    }
    
    if (params.do_gauge_transform)
    {
        // Generate and perform a random gauge transformation.
        rand_trans_u1(gauge_trans, Lat.get_lattice_dimension(0), Lat.get_lattice_dimension(1), generator);
        apply_gauge_trans_u1(lattice, gauge_trans, Lat.get_lattice_dimension(0), Lat.get_lattice_dimension(1));
        cout << "[GAUGE]: Performed a random gauge rotation.\n";
    }
    
    cout << "[GAUGE]: The average plaquette is " << get_plaquette_u1(lattice, Lat.get_lattice_dimension(0), Lat.get_lattice_dimension(1)) << ".\n";
    
    // Sanity check if we're doing a multigrid solve.
    if (params. n_refine == 0)
    {
        my_test = TOP_LEVEL_ONLY;
    }
    
    string op_null_name;
    void (*op_null)(complex<double>*, complex<double>*, void*) = square_staggered_u1; 
    switch (params.nvec_params.opt_null)
    {
        case STAGGERED:
            op_null = square_staggered_u1;
            op_null_name = "Staggered U(1)";
            break;
        case LAPLACE:
            op_null_name = "Free Laplace";
            op_null = square_laplace;
            break;
        case LAPLACE_NC2:
            op_null_name = "Free Laplace Nc = 2";
            op_null = square_laplace;
            break;
        case G5_STAGGERED:
            op_null_name = "Gamma_5 Staggered U(1)";
            op_null = square_staggered_gamma5_u1;
            break;
        case STAGGERED_NORMAL:
            op_null_name = "Staggered U(1) Normal";
            op_null = square_staggered_normal_u1;
            break;
        case STAGGERED_INDEX:
            op_name = "Staggered U(1) Index Operator";
            op = staggered_index_operator;
            break;
    }
    
    cout << "[OP]: Null Gen Operator " << op_null_name << "\n";
    
    // Build an mg_struct.
    mg_operator_struct_complex mgstruct;
    mgstruct.x_fine = Lat.get_lattice_dimension(0);
    mgstruct.y_fine = Lat.get_lattice_dimension(1); 
    mgstruct.Nc = Nc; // only matters for square laplace.
    mgstruct.n_refine = params.n_refine; 
    mgstruct.dslash_count = new dslash_tracker(params.n_refine); 
    
    if ((int)params.blocksizes.size() != params.n_refine && params.blocksizes.size() != 1 && params.n_refine > 0)
    {
        cout << "[ERROR]: Incorrect number of block sizes supplied. Needs to be either 1 or nrefine.\n";
        return 0;
    }
    
    mgstruct.blocksize_x = new int[params.n_refine];
    mgstruct.blocksize_y = new int[params.n_refine];
    if (params.blocksizes.size() == 1)
    {
        for (i = 0; i < params.n_refine; i++)
        {
            mgstruct.blocksize_x[i] = params.blocksizes[0];
            mgstruct.blocksize_y[i] = params.blocksizes[0];
        }
    }
    else // there are unique blocksizes for each level.
    {
        for (i = 0; i < params.n_refine; i++)
        {
            mgstruct.blocksize_x[i] = params.blocksizes[i];
            mgstruct.blocksize_y[i] = params.blocksizes[i];
        }
    }
    
    //mgstruct.blocksize_x[0] = 8;
    //mgstruct.blocksize_y[0] = 8;
    //for (i = 1; i < n_refine; i++)
    
    // If we do the free test, we construct a constant vector.
    if (params.do_free)
    {
        params.nvec_params.n_null_vector = 1;
    }
    
    switch (params.nvec_params.bstrat)
    {
        case BLOCK_NONE:
            mgstruct.n_vector = params.nvec_params.n_null_vector;
            params.nvec_params.null_partitions = 1;
            break;
        case BLOCK_EO:
        case BLOCK_TOPO:
            mgstruct.n_vector = 2*params.nvec_params.n_null_vector;
            params.nvec_params.null_partitions = 2;
            break;
        case BLOCK_CORNER:
            mgstruct.n_vector = 4*params.nvec_params.n_null_vector;
            params.nvec_params.null_partitions = 4; 
            break;
    }
    
    mgstruct.matrix_vector = op; //square_staggered_u1;
    mgstruct.matrix_vector_dagger = op_dagger; 
    mgstruct.matrix_extra_data = (void*)&stagif; 
    
    
    cout << "[MG]: X_Block ";
    for (i = 0; i < params.n_refine; i++)
    {
        cout << mgstruct.blocksize_x[i] << " ";
    }
    cout << "Y_Block ";
    for (i = 0; i < params.n_refine; i++)
    {
        cout << mgstruct.blocksize_y[i] << " "; 
    }
    cout << "NullVectors " << params.nvec_params.n_null_vector << "\n";
    
    // Build lattice objects for each level. 
    mgstruct.latt = new Lattice*[params.n_refine+1];
    mgstruct.latt[0] = new Lattice(Lat); // The fine level already exists.
    
    if (params.n_refine > 0)
    {
        for (i = 1; i <= params.n_refine; i++)
        {
            lattice_size[0] = mgstruct.latt[i-1]->get_lattice_dimension(0)/mgstruct.blocksize_x[i-1];
            lattice_size[1] = mgstruct.latt[i-1]->get_lattice_dimension(1)/mgstruct.blocksize_y[i-1];
            mgstruct.latt[i] = new Lattice(nd, lattice_size, mgstruct.n_vector);
        }
    }
    
    // Print info about lattices.
    for (i = 0; i <= params.n_refine; i++)
    {
        cout << "[LATTICE_L" << i+1 << "]: X " << mgstruct.latt[i]->get_lattice_dimension(0) <<
                " Y: " << mgstruct.latt[i]->get_lattice_dimension(1) <<
                " Vol: " << mgstruct.latt[i]->get_volume() <<
                " Nc: " << mgstruct.latt[i]->get_nc() << 
                " Size: " << mgstruct.latt[i]->get_lattice_size() << "\n";
    }
    
    // Set the starting mg_struct state.
    mgstruct.curr_level = 0; // Ready to do top level -> second level.
    mgstruct.curr_dof_fine = mgstruct.latt[0]->get_nc(); // Top level has only one d.o.f. per site. 
    mgstruct.curr_x_fine = mgstruct.latt[0]->get_lattice_dimension(0);
    mgstruct.curr_y_fine = mgstruct.latt[0]->get_lattice_dimension(1);
    mgstruct.curr_fine_size = mgstruct.latt[0]->get_lattice_size();
    
    if (params.n_refine > 0)
    {
        mgstruct.curr_dof_coarse = mgstruct.latt[1]->get_nc();
        mgstruct.curr_x_coarse = mgstruct.latt[1]->get_lattice_dimension(0);
        mgstruct.curr_y_coarse = mgstruct.latt[1]->get_lattice_dimension(1);
        mgstruct.curr_coarse_size = mgstruct.latt[1]->get_lattice_size();
    }
    else // give it some garbage so it doesn't complain.
    {
        mgstruct.curr_dof_coarse = mgstruct.n_vector; 
        mgstruct.curr_x_coarse = mgstruct.x_fine/mgstruct.blocksize_x[0];
        mgstruct.curr_y_coarse = mgstruct.y_fine/mgstruct.blocksize_y[0];
        mgstruct.curr_coarse_size = mgstruct.curr_y_coarse*mgstruct.curr_x_coarse*mgstruct.curr_dof_coarse;
    }
    
    
    // Build the mg inverter structure.
    // Set up the MG preconditioner. 
    mg_precond_struct_complex mgprecond;

    mgprecond.in_smooth_type = params.in_smooth; // What inner smoother? MinRes or GCR.
    mgprecond.omega_smooth = params.omega_smooth; // What relaxation parameter should we use (MinRes only!)
    mgprecond.mlevel_type = params.mlevel_type; // Do we smooth then go down, or smooth then start a new Krylov?
    mgprecond.in_solve_type = params.in_solve; // What inner solver? NONE, MINRES, CG, GCR, BICGSTAB. Should also be used for recursive solve!!
    mgprecond.n_max = params.inner_max; // max number of steps to use for inner solver.
    mgprecond.n_restart = params.inner_restart; // frequency of restart (relevant for CG, GCR).
    mgprecond.mgstruct = &mgstruct; // Contains null vectors, fine operator. (Since we don't construct the fine op.)
    mgprecond.matrix_extra_data = (void*)&mgstruct; // What extra_data the coarse operator expects. 
    
    // Inner precision. Default 1e-2. Maximum relative residual for coarse solves.
    if ((int)params.inner_precisions.size() != params.n_refine && params.inner_precisions.size() != 1 && params.n_refine > 0)
    {
        cout << "[ERROR]: Incorrect number of inner precisions supplied. Needs to be either 1 or nrefine.\n";
        return 0;
    }
    
    mgprecond.rel_res = new double[params.n_refine];
    
    if (params.inner_precisions.size() == 1) { for (i = 0; i < params.n_refine; i++) { mgprecond.rel_res[i] = params.inner_precisions[0]; } }
    else // there are unique pre smooth counts for each level.
    { for (i = 0; i < params.n_refine; i++) { mgprecond.rel_res[i] = params.inner_precisions[i]; } }
    
    
    // Presmooth and postsmooth. Default 6 MinRes pre, post. 
    if ((int)params.pre_smooths.size() != params.n_refine && params.pre_smooths.size() != 1 && params.n_refine > 0)
    {
        cout << "[ERROR]: Incorrect number of presmoother iterations supplied. Needs to be either 1 or nrefine.\n";
        return 0;
    }
    if ((int)params.post_smooths.size() != params.n_refine && params.post_smooths.size() != 1 && params.n_refine > 0)
    {
        cout << "[ERROR]: Incorrect number of postsmoother iterations supplied. Needs to be either 1 or nrefine.\n";
        return 0;
    }
    
    mgprecond.n_pre_smooth = new int[params.n_refine];
    mgprecond.n_post_smooth = new int[params.n_refine];
    
    if (params.pre_smooths.size() == 1) { for (i = 0; i < params.n_refine; i++) { mgprecond.n_pre_smooth[i] = params.pre_smooths[0]; } }
    else // there are unique pre smooth counts for each level.
    { for (i = 0; i < params.n_refine; i++) { mgprecond.n_pre_smooth[i] = params.pre_smooths[i]; } }
    
    if (params.post_smooths.size() == 1) { for (i = 0; i < params.n_refine; i++) { mgprecond.n_post_smooth[i] = params.post_smooths[0]; } }
    else // there are unique post smooth counts for each level.
    { for (i = 0; i < params.n_refine; i++) { mgprecond.n_post_smooth[i] = params.post_smooths[i]; } }
    
    mgprecond.normal_eqn_smooth = params.normal_eqn_smooth;
    mgprecond.normal_eqn_mg = params.normal_eqn_mg; 
    if (params.normal_eqn_mg) // MG where we use D^\dagger_coarse D_coarse
    {
        mgprecond.coarse_matrix_vector = coarse_square_staggered_normal; // Function which applies the coarse operator. 
        mgprecond.fine_matrix_vector = fine_square_staggered_normal; // Function which applies the fine operator. 
    
        if (params.normal_eqn_smooth) // Sort of pointless. 
        {
            mgprecond.coarse_matrix_vector_dagger = coarse_square_staggered_normal; // function which applies coarse dagger.
            mgprecond.fine_matrix_vector_dagger = fine_square_staggered_normal; // function which applies fine dagger.
            mgprecond.coarse_matrix_vector_normal = coarse_square_staggered_normal; // function which applies normal of coarse op
            mgprecond.fine_matrix_vector_normal = fine_square_staggered_normal; // function which applies normal of fine op.
        }
    }
    else // Standard MG.
    {
        mgprecond.coarse_matrix_vector = coarse_square_staggered; // Function which applies the coarse operator. 
        mgprecond.fine_matrix_vector = fine_square_staggered; // Function which applies the fine operator. 
    
        if (params.normal_eqn_smooth)
        {
            mgprecond.coarse_matrix_vector_dagger = coarse_square_staggered_dagger; // function which applies coarse dagger.
            mgprecond.fine_matrix_vector_dagger = fine_square_staggered_dagger; // function which applies fine dagger.
            mgprecond.coarse_matrix_vector_normal = coarse_square_staggered_normal; // function which applies normal of coarse op
            mgprecond.fine_matrix_vector_normal = fine_square_staggered_normal; // function which applies normal of fine op.
        }
    }

    // End set up MG preconditioners

    // Set right hand side.
    switch (source)
    {
        case POINT: // Set a point.
            for (i = 0; i < Nc; i++)
            {
                rhs[(Lat.get_lattice_dimension(0)/2+(Lat.get_lattice_dimension(1)/2)*Lat.get_lattice_dimension(0))*Nc+i] = 1.0;
            }
            break;
        case RANDOM_GAUSSIAN: // Random rhs.
            gaussian<double>(rhs, Lat.get_lattice_size(), generator);
            break;
        case ORIGIN_POINT: // Set a point for correlator computation.
            for (i = 0; i < Nc; i++)
            {
                rhs[i] = 1.0;
            }
    }

    // Get norm for rhs.
    bnorm = sqrt(norm2sq<double>(rhs, Lat.get_lattice_size()));

    // Set a point on the lhs.
    //lhs[x_fine/2+(y_fine/2)*x_fine+1] = 1.0;
    
    // Allocate null vectors.
    mgstruct.null_vectors = new complex<double>**[mgstruct.n_refine];
    for (i = 0; i < mgstruct.n_refine; i++)
    {
        mgstruct.null_vectors[i] = new complex<double>*[mgstruct.n_vector];
        for (j = 0; j < mgstruct.n_vector; j++)
        {
            mgstruct.null_vectors[i][j] = new complex<double>[mgstruct.latt[i]->get_lattice_size()];
            zero<double>(mgstruct.null_vectors[i][j], mgstruct.latt[i]->get_lattice_size());
        }
    }
    
    // Repeat null precisions if needed.
    if (params.nvec_params.null_precisions.size() == 1 && params.n_refine > 1)
    {
        for (i = 1; i < params.n_refine; i++)
        {
            params.nvec_params.null_precisions[i] = params.nvec_params.null_precisions[0];
        }
    }
    
    // Allocate null max iters
    if (params.nvec_params.null_max_iters.size() == 1 && params.n_refine > 1)
    {
        for (i = 1; i < params.n_refine; i++)
        {
            params.nvec_params.null_max_iters[i] = params.nvec_params.null_max_iters[0];
        }
    }
    
    // Allocate stencils.
    mgstruct.stencils = new stencil_2d*[params.n_refine+1];
    for (i = 0; i <= mgstruct.n_refine; i++)
    {
        // We can't allocate a stencil if the lattice gets too small.
        if (get_stencil_size(params.opt) == 2 && (mgstruct.latt[i]->get_lattice_dimension(0) < 4 || mgstruct.latt[i]->get_lattice_dimension(1) < 4))
        {
            mgstruct.stencils[i] = 0; // We can't generate a stencil. It's small enough that we just explicitly project.
        }
        else if (mgstruct.latt[i]->get_lattice_dimension(0) < 2 || mgstruct.latt[i]->get_lattice_dimension(1) < 2) // it's a stencil size 1 or a big enough stencil size 2.
        {
            mgstruct.stencils[i] = 0;
        }
        else // We can build a stencil.
        {
            mgstruct.stencils[i] = new stencil_2d(mgstruct.latt[i], get_stencil_size(params.opt));
        }
    }
    
    // Allocate stencils for daggered operator, if needed.
    if (params.normal_eqn_smooth || params.normal_eqn_mg) 
    {
        mgstruct.have_dagger_stencil = true;
        
        mgstruct.dagger_stencils = new stencil_2d*[params.n_refine+1];
        for (i = 0; i <= mgstruct.n_refine; i++)
        {
            // We can't allocate a stencil if the lattice gets too small.
            if (get_stencil_size(params.opt) == 2 && (mgstruct.latt[i]->get_lattice_dimension(0) < 4 || mgstruct.latt[i]->get_lattice_dimension(1) < 4))
            {
                mgstruct.dagger_stencils[i] = 0; // We can't generate a stencil. It's small enough that we just explicitly project.
            }
            else if (mgstruct.latt[i]->get_lattice_dimension(0) < 2 || mgstruct.latt[i]->get_lattice_dimension(1) < 2) // it's a stencil size 1 or a big enough stencil size 2.
            {
                mgstruct.dagger_stencils[i] = 0;
            }
            else // We can build a stencil.
            {
                mgstruct.dagger_stencils[i] = new stencil_2d(mgstruct.latt[i], get_stencil_size(params.opt));
            }
        }
    }
    else
    {
        mgstruct.have_dagger_stencil = false;
        mgstruct.dagger_stencils = 0;
    }
    
    if (!(my_test == TOP_LEVEL_ONLY || my_test == SMOOTHER_ONLY))
    {
        cout << "[MG]: Creating " << mgstruct.n_vector << " null vectors.\n" << flush; 
        
        // Generate top level!
        cout << "[NULLVEC]: About to generate null vectors.\n" << flush; 
        
        mgstruct.matrix_vector = op_null; // Trick op to null gen op.
        
        // Back up verbosity string.
        std::string verb_string = verb.verb_prefix;
        
        // Back up stencils.
        stencil_2d** save_stencil = mgstruct.stencils;
        
        // Allocate stencils for null vector generation. 
        mgstruct.stencils = new stencil_2d*[params.n_refine];
        for (i = 0; i < mgstruct.n_refine; i++)
        {
            // We can't allocate a stencil if the lattice gets too small.
            if (get_stencil_size(params.nvec_params.opt_null) == 2 && (mgstruct.latt[i]->get_lattice_dimension(0) < 4 || mgstruct.latt[i]->get_lattice_dimension(1) < 4))
            {
                mgstruct.stencils[i] = 0; // We can't generate a stencil. It's small enough that we just explicitly project.
            }
            else if (mgstruct.latt[i]->get_lattice_dimension(0) < 2 || mgstruct.latt[i]->get_lattice_dimension(1) < 2) // it's a stencil size 1 or a big enough stencil size 2.
            {
                mgstruct.stencils[i] = 0;
            }
            else // We can build a stencil.
            {
                mgstruct.stencils[i] = new stencil_2d(mgstruct.latt[i], get_stencil_size(params.nvec_params.opt_null));
            }
        }
        
        // Fine level is handled specially. In some case, we have a special function to
        // populate the fine staggered stencil without a bunch of matrix multiplies. We don't have it for all
        // operators.
        switch (params.nvec_params.opt_null)
        {
            case STAGGERED: 
                // Encode mass in shift, not in the stencil.
                tmp_mass = stagif.mass; stagif.mass = 0.0;
                get_square_staggered_u1_stencil(mgstruct.stencils[0], &stagif);
                stagif.mass = tmp_mass;
                mgstruct.stencils[0]->shift = params.nvec_params.null_mass;
                break;
            case G5_STAGGERED:
                // Encode mass in a phase shift, not in the stencil.
                tmp_mass = stagif.mass; stagif.mass = 0.0;
                get_square_staggered_gamma5_u1_stencil(mgstruct.stencils[0], &stagif);
                stagif.mass = tmp_mass;
                mgstruct.stencils[0]->eo_shift = params.nvec_params.null_mass; 
                break;
            case LAPLACE: // not implemented yet...
            case LAPLACE_NC2:
                // Encode mass in shift, not in the stencil.
                tmp_mass = stagif.mass; stagif.mass = 0.0;
                generate_stencil_2d(mgstruct.stencils[0], fine_square_staggered, (void*)&mgstruct);
                stagif.mass = tmp_mass;
                mgstruct.stencils[0]->shift = params.nvec_params.null_mass; 
                break;
            case STAGGERED_NORMAL:
            case STAGGERED_INDEX:
                // Use the null generation mass. 
                tmp_mass = stagif.mass; stagif.mass = params.nvec_params.null_mass;
                generate_stencil_2d(mgstruct.stencils[0], fine_square_staggered, (void*)&mgstruct);
                stagif.mass = tmp_mass;
                break;
        }
        cout << "[STENCIL_L1]: Built stencil.\n" << flush;
        
        for (int n = 0; n < mgstruct.n_refine; n++)
        {   
            cout << "[L" << mgstruct.curr_level+1 << "_NULLVEC]: Current fine size is " << mgstruct.curr_fine_size << "\n";
            
            // Update verb_prefix temporarily.
            verb.verb_prefix = "[L" + to_string(mgstruct.curr_level+1) + "_NULLVEC]: ";
            
            if (params.do_free)
            {
                null_generate_free(&mgstruct, &params.nvec_params, params.do_gauge_transform, gauge_trans);
            }
            else if (params.nvec_use_eigen) // use ARPACK. It can't get here if we don't have EIGEN_TEST
            {
                #ifdef EIGEN_TEST 
                int n_eigen = mgstruct.n_vector/params.nvec_params.null_partitions;
                int n_cv = min(10*mgstruct.n_vector/params.nvec_params.null_partitions, mgstruct.curr_fine_size);
                arpack_dcn_t* ar_strc = arpack_dcn_init(mgstruct.curr_fine_size, n_eigen, n_cv); 
                char eigtype[3]; strcpy(eigtype, "SM"); // Smallest magnitude eigenvalues.
                complex<double>* eigs_tmp = new complex<double>[mgstruct.n_vector/params.nvec_params.null_partitions];
                arpack_solve_t info_solve = arpack_dcn_getev(ar_strc, eigs_tmp, mgstruct.null_vectors[mgstruct.curr_level], mgstruct.curr_fine_size, n_eigen, n_cv, params.nvec_params.null_max_iters[n], eigtype, params.nvec_params.null_precisions[n], 0.0, fine_square_staggered, (void*)&mgstruct); 
                delete[] eigs_tmp; 
                arpack_dcn_free(&ar_strc);
                
                 // Print info about the eigensolve.
                cout << "[L" << mgstruct.curr_level+1 << "_NULLVEC_ARPACK]: Number of converged eigenvalues: " << info_solve.nconv << "\n";
                cout << "[L" << mgstruct.curr_level+1 << "_NULLVEC_ARPACK]: Number of iteration steps: " << info_solve.niter << "\n";
                cout << "[L" << mgstruct.curr_level+1 << "_NULLVEC_ARPACK]: Number of matrix multiplies: " << info_solve.nops << "\n";
                
                for (i = 0; i < mgstruct.n_vector/params.nvec_params.null_partitions; i++)
                {
                    // Aggregate in chirality (or corners) as needed.  
                    // This is handled differently if we're on the top level or further down. 
                    if (mgstruct.curr_level == 0) // top level routines
                    {
                        null_partition_staggered(&mgstruct, i, params.nvec_params.bstrat, &Lat);
                    }
                    else // not on the top level
                    {
                        null_partition_coarse(&mgstruct, i, params.nvec_params.bstrat);
                    }

                    for (k = 0; k < params.nvec_params.null_partitions; k++)
                    {
                        normalize(mgstruct.null_vectors[mgstruct.curr_level][i+k*params.nvec_params.n_null_vector], mgstruct.curr_fine_size);
                    }
                }
#endif // EIGEN_TEST

            }
            else // Generate null vectors. 
            {
                null_generate_random_smooth(&mgstruct, &params.nvec_params, &verb, &generator);
            }



    #ifdef PRINT_NULL_VECTOR
            // Print vector.
            
            cout << "\n\nPrinting null vectors:\n"; 
            for (int n = 0; n < mgstruct.n_vector; n++)
            {
                cout << "\nVector " << n << "\n";
                for (int y = 0; y < mgstruct.curr_y_fine; y++)
                {
                    for (int x = 0; x < mgstruct.curr_x_fine; x++)
                    {
                        cout << "(";
                        for (int c = 0; c < mgstruct.curr_dof_fine; c++)
                        {
                            cout << mgstruct.null_vectors[mgstruct.curr_level][n][y*mgstruct.curr_x_fine*mgstruct.curr_dof_fine+x*mgstruct.curr_dof_fine+c] << ",";
                        }
                        cout << ") ";
                    }
                    cout << "\n";
                }
            }
    #endif // PRINT_NULL_VECTOR

            cout << "[MG]: Performing block orthonormalize of null vectors...\n";
            block_orthonormalize(&mgstruct); 
            
            // Only need to build stencil for null op.
            if (n != mgstruct.n_refine-1)
            {
                
                // Build a temporary coarse stencil if we can! We'll have to do a rebuild later
                // to use the correct mass.

                // First, make sure we CAN build the coarse stencil. When we tried to allocate stencil objects, 
                // it failed if the volume wasn't big enough.
                if (mgstruct.stencils[mgstruct.curr_level+1] != 0)
                {
                    if (get_stencil_size(params.nvec_params.opt_null) == 2)
                    {
                        // Generate the old, inefficient way. This method builds the shifts directly in.
                        generate_stencil_2d(mgstruct.stencils[mgstruct.curr_level+1], coarse_square_staggered, (void*)&mgstruct);
                        cout << "[STENCIL_L" << mgstruct.curr_level+2 << "]: Null vector generation inefficient stencil build.\n" << flush;
                    }
                    else
                    {
                        // Generate the new, shiny way!
                        if (params.nvec_params.bstrat == BLOCK_TOPO || (params.nvec_params.opt_null == G5_STAGGERED && (params.nvec_params.bstrat == BLOCK_NONE || params.nvec_params.bstrat == BLOCK_TOPO)))
                        {
                            // Need to block mass in since we've destroyed good chirality.
                            generate_coarse_from_fine_stencil(mgstruct.stencils[mgstruct.curr_level+1], mgstruct.stencils[mgstruct.curr_level], &mgstruct, false);  // build shifts in.
                        }
                        else // we're good!
                        {
                            generate_coarse_from_fine_stencil(mgstruct.stencils[mgstruct.curr_level+1], mgstruct.stencils[mgstruct.curr_level], &mgstruct, true);  // true -> ignore shifts, don't build in.
                            // Set shifted masses.
                            switch (params.nvec_params.opt_null)
                            {
                                case STAGGERED: 
                                case LAPLACE:
                                case LAPLACE_NC2:
                                    mgstruct.stencils[mgstruct.curr_level+1]->shift = mgstruct.stencils[mgstruct.curr_level]->shift;
                                    break;
                                case G5_STAGGERED:
                                    mgstruct.stencils[mgstruct.curr_level+1]->dof_shift = (mgstruct.curr_level == 0 ? mgstruct.stencils[mgstruct.curr_level]->eo_shift : mgstruct.stencils[mgstruct.curr_level]->dof_shift);
                                    break;
                                case STAGGERED_NORMAL:
                                case STAGGERED_INDEX:
                                    // Can't get here anyway.
                                    break;
                            }
                        }
                        cout << "[STENCIL_L" << mgstruct.curr_level+2 << "]: Null vector generation efficient stencil build.\n" << flush;
                    }
                }
                else
                {
                    cout << "[STENCIL_L" << mgstruct.curr_level+2 << "]: Lattice too small to build null vector generation stencil.\n" << flush;
                }
                
                level_down(&mgstruct);
            }
            
        }
        
        // Reset the tricks we've done. 
        
        // Clear out temporary null vector stencils.
        for (i = 1; i < mgstruct.n_refine; i++)
        {
            if (mgstruct.stencils[i] != 0) { delete mgstruct.stencils[i]; }
        }
        delete mgstruct.stencils[0];
        delete[] mgstruct.stencils; 
        
        mgstruct.matrix_vector = op; // Reset op to solved op.
        verb.verb_prefix = verb_string; 
        mgstruct.stencils = save_stencil;
        
        // Un-pop to the finest level.
        for (int n = 1; n < mgstruct.n_refine; n++)
        {
            level_up(&mgstruct);
        }
        
        // We built the previous stencils with the generation mass! Rebuild with the correct mass.
        
        // Clear and rebuild.
        
        // First, build the top level.
        mgstruct.stencils[0]->clear_stencils();
        switch (params.opt)
        {
            case STAGGERED:
                // Encode mass in shift, not in the stencil.
                tmp_mass = stagif.mass; stagif.mass = 0.0;
                get_square_staggered_u1_stencil(mgstruct.stencils[0], &stagif);
                if (mgstruct.have_dagger_stencil)
                {
                    get_square_staggered_dagger_u1_stencil(mgstruct.dagger_stencils[0], &stagif);
                }
                stagif.mass = tmp_mass;
                mgstruct.stencils[0]->shift = stagif.mass;
                if (mgstruct.have_dagger_stencil)
                {
                    mgstruct.dagger_stencils[0]->shift = stagif.mass;
                }
                cout << "[STENCIL_L1]: Final efficient stencil build.\n" << flush;
                break;
            case G5_STAGGERED:
                // Encode mass in shift, not in the stencil.
                tmp_mass = stagif.mass; stagif.mass = 0.0;
                get_square_staggered_gamma5_u1_stencil(mgstruct.stencils[0], &stagif);
                if (mgstruct.have_dagger_stencil)
                {
                    // G5_STAGGERED is hermitian indefinite.
                    get_square_staggered_gamma5_u1_stencil(mgstruct.dagger_stencils[0], &stagif);
                }
                stagif.mass = tmp_mass;
                mgstruct.stencils[0]->eo_shift = stagif.mass;
                if (mgstruct.have_dagger_stencil)
                {
                    mgstruct.dagger_stencils[0]->eo_shift = stagif.mass;
                }
                cout << "[STENCIL_L1]: Final efficient stencil build.\n" << flush;
                break;
            case LAPLACE: // not implemented yet...
            case LAPLACE_NC2:
                // Encode mass in shift, not in the stencil.
                tmp_mass = stagif.mass; stagif.mass = 0.0;
                generate_stencil_2d(mgstruct.stencils[0], fine_square_staggered, (void*)&mgstruct);
                if (mgstruct.have_dagger_stencil)
                {
                    // These operators are all hermitian (pos def or indef)
                    generate_stencil_2d(mgstruct.dagger_stencils[0], fine_square_staggered, (void*)&mgstruct);
                }
                stagif.mass = tmp_mass;
                mgstruct.stencils[0]->shift = stagif.mass;
                if (mgstruct.have_dagger_stencil)
                {
                    mgstruct.dagger_stencils[0]->shift = stagif.mass;
                }
                cout << "[STENCIL_L1]: Final inefficient stencil build.\n" << flush;
                break; 
            case STAGGERED_NORMAL:
            case STAGGERED_INDEX:
                generate_stencil_2d(mgstruct.stencils[0], fine_square_staggered, (void*)&mgstruct);
                if (mgstruct.have_dagger_stencil)
                {
                    // These operators are all hermitian (pos def or indef)
                    generate_stencil_2d(mgstruct.dagger_stencils[0], fine_square_staggered, (void*)&mgstruct);
                }
                cout << "[STENCIL_L1]: Final inefficient stencil build.\n" << flush;
                break;
        }
        
        for (int n = 0; n < mgstruct.n_refine; n++)
        {
            if (mgstruct.stencils[n+1] != 0)
            {
                mgstruct.stencils[n+1]->clear_stencils();
                if (get_stencil_size(params.opt) == 2)
                {
                    // Generate the old, inefficient way. This builds shifts directly in.
                    generate_stencil_2d(mgstruct.stencils[n+1], coarse_square_staggered, (void*)&mgstruct);
                    if (mgstruct.have_dagger_stencil)
                    {
                        generate_stencil_2d(mgstruct.dagger_stencils[n+1], coarse_square_staggered_dagger, (void*)&mgstruct);
                    }
                    cout << "[STENCIL_L" << n+2 << "]: Final inefficient stencil build.\n" << flush;
                }
                else
                {
                    // Generate the new, shiny way!
                    if (params.nvec_params.bstrat == BLOCK_TOPO || (params.opt == G5_STAGGERED && params.nvec_params.bstrat == BLOCK_NONE))
                    {
                        // Need to block mass in since we've destroyed good chirality.
                        generate_coarse_from_fine_stencil(mgstruct.stencils[n+1], mgstruct.stencils[n], &mgstruct, false);  // build shifts in.
                    }
                    else // we're good!
                    {
                        generate_coarse_from_fine_stencil(mgstruct.stencils[n+1], mgstruct.stencils[n], &mgstruct, true);  // true -> ignore shifts, don't build in.
                        // Set shifted masses.
                        switch (params.opt)
                        {
                            case STAGGERED: 
                            case LAPLACE:
                            case LAPLACE_NC2:
                                mgstruct.stencils[n+1]->shift = mgstruct.stencils[n]->shift; 
                                break;
                            case G5_STAGGERED:
                                mgstruct.stencils[n+1]->dof_shift = (n == 0 ? mgstruct.stencils[n]->eo_shift : mgstruct.stencils[n]->dof_shift);
                                break;
                            case STAGGERED_NORMAL:
                            case STAGGERED_INDEX:
                                // Can't get here anyway.
                                break;
                        }
                    }
                    
                    
                    // Generate the new, shiny way!
                    //generate_coarse_from_fine_stencil(mgstruct.stencils[n+1], mgstruct.stencils[n], &mgstruct); 
                    if (mgstruct.have_dagger_stencil)
                    {
                        // Generate the new, shiny way!
                        if (params.nvec_params.bstrat == BLOCK_TOPO || (params.opt == G5_STAGGERED && (params.nvec_params.bstrat == BLOCK_NONE || params.nvec_params.bstrat == BLOCK_TOPO)))
                        {
                            // Need to block mass in since we've destroyed good chirality.
                            generate_coarse_from_fine_stencil(mgstruct.dagger_stencils[n+1], mgstruct.dagger_stencils[n], &mgstruct, false);  // build shifts in.
                        }
                        else // we're good!
                        {
                            generate_coarse_from_fine_stencil(mgstruct.dagger_stencils[n+1], mgstruct.dagger_stencils[n], &mgstruct, true);  // true -> ignore shifts, don't build in.
                            // Set shifted masses.
                            switch (params.opt)
                            {
                                case STAGGERED: 
                                case LAPLACE:
                                case LAPLACE_NC2:
                                    mgstruct.dagger_stencils[n+1]->shift = mgstruct.dagger_stencils[n]->shift; 
                                    break;
                                case G5_STAGGERED:
                                    mgstruct.stencils[n+1]->dof_shift = (n == 0 ? mgstruct.dagger_stencils[n]->eo_shift : mgstruct.dagger_stencils[n]->dof_shift);
                                    break;
                                case STAGGERED_NORMAL:
                                case STAGGERED_INDEX:
                                    // Can't get here anyway.
                                    break;
                            }
                        }
                        
                    }
                    cout << "[STENCIL_L" << n+2 << "]: Final efficient stencil build.\n" << flush;
                }
            }
            else
            {
                cout << "[STENCIL_L" << n+2 << "]: Lattice too small to build stencil.\n" << flush;
            }
            
            if (n != mgstruct.n_refine-1)
            {
                level_down(&mgstruct);
            }
        }
        
        // Un-pop to the finest level.
        for (int n = 1; n < mgstruct.n_refine; n++)
        {
            level_up(&mgstruct);
        }
        
    } // end skipping generation if we're only doing a top level or smoother test. 
    
    
#ifdef PRINT_NULL_VECTOR
    cout << "\nCheck projector:\n"; 
    for (int n = 0; n < mgstruct.n_vector; n++)
    {
        cout << "Vector " << n << "\n";
        for (int y = 0; y < y_fine; y++)
        {
            for (int x = 0; x < x_fine; x++)
            {
                cout << mgstruct.null_vectors[0][n][x+y*x_fine] << " ";
            }
            cout << "\n";
        }
    }
#endif // PRINT_NULL_VECTOR
    
    
#ifdef EIGEN_TEST
    
    if (params.do_eigentest)
    {   
        test_eigenvalue_overlap(&mgstruct, &stagif, params.opt, params.set_eigen, params.set_cv); 
    } // do_eigentest
#endif // EIGEN_TEST
    
#ifdef STENCIL_CONSTRUCT_TEST
    for (i = 0; i <= mgstruct.n_refine; i++)
    {
        test_stencil_construct(&mgstruct, i, get_stencil_size(params.opt)); 
    }

    return 0; 
#endif

    
#ifdef STENCIL_PIECE_TEST
    for (i = 0; i <= mgstruct.n_refine; i++)
    {
        test_stencil_piece(&mgstruct, i, get_stencil_size(opt)); 
    }

    return 0; 
#endif
    
#ifdef STENCIL_EFFICIENT_TEST
    // Do an inefficient build of a stencil for the fine operator.
    stencil_2d stenc_fine(mgstruct.latt[0], get_stencil_size(opt));
    
    generate_stencil_2d(&stenc_fine, fine_square_staggered, (void*)&mgstruct);
    
    // Good! We've built the first level.
    
    // Next, build the coarse operator.
    stencil_2d stenc_coarse(mgstruct.latt[1], get_stencil_size(opt));
    cout << "About to generate coarse stencil.\n" << flush;
    generate_coarse_from_fine_stencil(&stenc_coarse, &stenc_fine, &mgstruct); 
    cout << "Generated coarse stenci.\n" << flush; 
    
    /*for (i = 0; i < mgstruct.latt[1]->get_nc()*mgstruct.latt[1]->get_nc(); i++)
    {
        cout << stenc_coarse.hopping[3*mgstruct.latt[1]->get_lattice_size()*mgstruct.latt[1]->get_nc()+i] << " ";
        if (i % mgstruct.latt[1]->get_nc() == mgstruct.latt[1]->get_nc()-1)
        {
            cout << "\n";
        }
    }
    cout << "\n\n" << flush; */
    
    // Generate the coarse operator the old way.
    stencil_2d stenc_coarse_old(mgstruct.latt[1], get_stencil_size(opt));
    generate_stencil_2d(&stenc_coarse_old, coarse_square_staggered, (void*)&mgstruct);
    
    /*for (i = 0; i < mgstruct.latt[1]->get_nc()*mgstruct.latt[1]->get_nc(); i++)
    {
        cout << stenc_coarse_old.hopping[3*mgstruct.latt[1]->get_lattice_size()*mgstruct.latt[1]->get_nc()+i] << " ";
        if (i % mgstruct.latt[1]->get_nc() == mgstruct.latt[1]->get_nc()-1)
        {
            cout << "\n";
        }
    }
    cout << "\n" << flush; */
    
    // Let's run a comparison!
    complex<double>* tmp_rhs = new complex<double>[mgstruct.latt[1]->get_lattice_size()];
    complex<double>* tmp_lhs = new complex<double>[mgstruct.latt[1]->get_lattice_size()];
    
    gaussian<double>(tmp_rhs, mgstruct.latt[1]->get_lattice_size(), generator); 
    apply_stencil_2d(tmp_lhs, tmp_rhs, (void*)&stenc_coarse);
    cout << "Applied coarse stencil.\n" << flush; 
    
    
    complex<double>* tmp_lhs2 = new complex<double>[mgstruct.latt[1]->get_lattice_size()];
    apply_stencil_2d(tmp_lhs2, tmp_rhs, (void*)&stenc_coarse_old);
    //coarse_square_staggered(tmp_lhs2, tmp_rhs, (void*)&mgstruct);
    cout << "Applied old coarse stencil.\n" << flush; 
    //for (i = 0; i < mgstruct.latt[1]->get_lattice_size(); i++)
    //{
    //    tmp_lhs2[i] = stagif.mass*tmp_rhs[i];
    //}
    
    // Get squared difference.
    cout << "[TEST]: Level " << 2 << " Squared difference: " << diffnorm2sq<double>(tmp_lhs, tmp_lhs2, mgstruct.latt[1]->get_lattice_size()) << "\n" << flush; 
    
    return 0;  
#endif
    
#ifdef STENCIL_DAGGER_TEST
    if (!(normal_eqn_smooth || normal_eqn_mg))
    {
        cout << "The dagger test can only be performed if '--smoother-type normal_eqn' is set (which builds the dagger stencil).\n" << flush;
        
        return 0;
    }
    
    // Compare applying D^\dag D via the stencils with via explicit prolong/restrict.
    complex<double>* tmp_rhs = new complex<double>[mgstruct.latt[1]->get_lattice_size()];
    complex<double>* tmp_lhs = new complex<double>[mgstruct.latt[1]->get_lattice_size()];
    complex<double>* tmp_lhs2 = new complex<double>[mgstruct.latt[1]->get_lattice_size()];
    complex<double>* tmp_tmp = new complex<double>[mgstruct.latt[1]->get_lattice_size()];
    
    complex<double>* tmp_rhs_fine = new complex<double>[mgstruct.latt[0]->get_lattice_size()];
    complex<double>* tmp_lhs_fine = new complex<double>[mgstruct.latt[0]->get_lattice_size()];
    
    // Randomize the rhs.
    gaussian<double>(tmp_rhs, mgstruct.latt[1]->get_lattice_size(), generator);
    
    // Test stencils.
    zero<double>(tmp_tmp, mgstruct.latt[1]->get_lattice_size());
    apply_stencil_2d(tmp_tmp, tmp_rhs, mgstruct.stencils[1]);
    zero<double>(tmp_lhs, mgstruct.latt[1]->get_lattice_size());
    apply_stencil_2d(tmp_lhs, tmp_tmp, mgstruct.dagger_stencils[1]);
    
    // Test prolong.
    zero<double>(tmp_rhs_fine, mgstruct.latt[0]->get_lattice_size());
    prolong(tmp_rhs_fine, tmp_rhs, &mgstruct);
    zero<double>(tmp_lhs_fine, mgstruct.latt[0]->get_lattice_size());
    op(tmp_lhs_fine, tmp_rhs_fine, (void*)&stagif);
    zero<double>(tmp_tmp, mgstruct.latt[1]->get_lattice_size());
    restrict(tmp_tmp, tmp_lhs_fine, &mgstruct);
    
    zero<double>(tmp_rhs_fine, mgstruct.latt[0]->get_lattice_size());
    prolong(tmp_rhs_fine, tmp_tmp, &mgstruct);
    zero<double>(tmp_lhs_fine, mgstruct.latt[0]->get_lattice_size());
    op_dagger(tmp_lhs_fine, tmp_rhs_fine, (void*)&stagif);
    zero<double>(tmp_lhs2, mgstruct.latt[1]->get_lattice_size());
    restrict(tmp_lhs2, tmp_lhs_fine, &mgstruct);
    
    // Get squared difference.
    cout << "[TEST]: Level " << 2 << " Squared difference: " << diffnorm2sq<double>(tmp_lhs, tmp_lhs2, mgstruct.latt[1]->get_lattice_size()) << "\n" << flush; 
    
    // Clean up.
    delete[] tmp_rhs;
    delete[] tmp_lhs;
    delete[] tmp_lhs2;
    delete[] tmp_tmp;
    delete[] tmp_rhs_fine;
    delete[] tmp_lhs_fine; 
    
    return 0;
#endif
    
    
#ifdef STENCIL_SHIFT_TEST
    
    for (int n = 0; n <= mgstruct.n_refine; n++)
    {
        complex<double>* tmp_rhs = new complex<double>[mgstruct.latt[n]->get_lattice_size()];
        complex<double>* tmp_lhs = new complex<double>[mgstruct.latt[n]->get_lattice_size()];
        complex<double>* tmp_lhs2 = new complex<double>[mgstruct.latt[n]->get_lattice_size()];

        // Randomize the rhs.
        gaussian<double>(tmp_rhs, mgstruct.latt[n]->get_lattice_size(), generator);

        // Apply op.
        if (n == 0)
        {
            op(tmp_lhs, tmp_rhs, &stagif);
        }
        else
        {
            apply_stencil_2d(tmp_lhs, tmp_rhs, mgstruct.stencils[n]);
        }
        
        // Generate stencil with shift built in.
        stencil_2d* tmp_stencil = new stencil_2d(mgstruct.latt[n], get_stencil_size(opt));
        if (n == 0)
        {
            switch (opt)
            {
                case STAGGERED:
                    get_square_staggered_u1_stencil(tmp_stencil, &stagif);
                    break;
                case G5_STAGGERED:
                    get_square_staggered_gamma5_u1_stencil(tmp_stencil, &stagif);
                    break;
                case LAPLACE: // not implemented yet...
                case LAPLACE_NC2:
                case STAGGERED_NORMAL:
                case STAGGERED_INDEX:
                    generate_stencil_2d(tmp_stencil, fine_square_staggered, (void*)&mgstruct);
                    break;
            }
        }
        else
        {
            generate_coarse_from_fine_stencil(tmp_stencil, mgstruct.stencils[n-1], &mgstruct, false);  // build shifts in.
            level_down(&mgstruct);
        }

        // Apply stencil.
        apply_stencil_2d(tmp_lhs2, tmp_rhs, tmp_stencil);

        // Compare
        cout << "Relative difference in test is " << diffnorm2sq<double>(tmp_lhs, tmp_lhs2, mgstruct.latt[n]->get_lattice_size()) << "\n" << flush;
        
        delete[] tmp_rhs;
        delete[] tmp_lhs;
        delete[] tmp_lhs2;
    }

    return 0;
    
#endif
    
#ifdef PDAGP_TEST
    {
        // Begin PdagP test.
        
        // Describe the coarse lattice. 
        int x_coarse = x_fine/mgstruct.blocksize_x[0]; // how many coarse sites are there in the x dir?
        int y_coarse = y_fine/mgstruct.blocksize_y[0]; // how many coarse sites are there in the y dir?
        int coarse_size = x_coarse*y_coarse; 


        // Print the fine rhs.
        cout << "Check fine point src:\n"; 
        for (int y = 0; y < y_fine; y++)
        {
            for (int x = 0; x < x_fine; x++)
            {
                cout << rhs[x+y*x_fine] << " ";
            }
            cout << "\n";
        }

        cout << "Check projector:\n"; 
        for (int y = 0; y < y_fine; y++)
        {
            for (int x = 0; x < x_fine; x++)
            {
                cout << mgstruct.null_vectors[0][0][x+y*x_fine] << " ";
            }
            cout << "\n";
        }

        // Test block normalizing the projector.
        cout << "Block normalize the projector(s). X blocks: " << x_coarse << " Y blocks: " << y_coarse << "\n";
        block_normalize(&mgstruct); 

        cout << "Check block normalized vector:\n"; 
        for (int y = 0; y < y_fine; y++)
        {
            for (int x = 0; x < x_fine; x++)
            {
                cout << mgstruct.null_vectors[0][0][x+y*x_fine] << " ";
            }
            cout << "\n";
        }


        // Restrict the original source. 
        complex<double>* rhs_coarse = new complex<double>[coarse_size*mgstruct.n_vector];
        cout << "Restricting the source.\n";
        restrict(rhs_coarse, rhs, &mgstruct);

        // Check coarse source.
        cout << "Check coarse src:\n"; 
        for (int y = 0; y < y_coarse; y++)
        {
            for (int x = 0; x < x_coarse; x++)
            {
                cout << rhs_coarse[x+y*x_coarse] << " ";
            }
            cout << "\n";
        }

        // Prolong the restricted source. 
        complex<double> *rhs_PdagP = new complex<double>[Lat.get_lattice_size()];
        cout << "Prolonging the restricted source.\n";
        prolong(rhs_PdagP, rhs_coarse, &mgstruct);

        cout << "Check re-prolonged source vector:\n"; 
        for (int y = 0; y < y_fine; y++)
        {
            for (int x = 0; x < x_fine; x++)
            {
                cout << rhs_PdagP[x+y*x_fine] << " ";
            }
            cout << "\n";
        }


        // Try applying the coarse Laplace operator!
        complex<double>* rhs_A_coarse = new complex<double>[coarse_size*mgstruct.n_vector]; 
        zero<double>(rhs_A_coarse, coarse_size*mgstruct.n_vector); 

        cout << "Applying the coarse Laplace operator to coarse source.\n";
        coarse_square_laplace(rhs_A_coarse, rhs_coarse, (void*)&mgstruct);

        // Check A coarse source.
        cout << "Check A times coarse src:\n"; 
        for (int y = 0; y < y_coarse; y++)
        {
            for (int x = 0; x < x_coarse; x++)
            {
                cout << rhs_A_coarse[x+y*x_coarse] << " ";
            }
            cout << "\n";
        }

        // Prolong rhs_A_coarse.
        cout << "Prolong A times coarse source.\n";
        complex<double>* rhs_PAP_fine = new complex<double>[Lat.get_lattice_size()];
        zero<double>(rhs_PAP_fine, Lat.get_lattice_size());
        prolong(rhs_PAP_fine, rhs_A_coarse, &mgstruct);

        // Check PAP. 
        cout << "Check PAP on source.\n";
        for (int y = 0; y < y_fine; y++)
        {
            for (int x = 0; x < x_fine; x++)
            {
                cout << rhs_PAP_fine[x+y*x_fine] << " ";
            }
            cout << "\n";
        }
        
        delete[] rhs_coarse;
        delete[] rhs_A_coarse;
        delete[] rhs_PdagP;
        delete[] rhs_PAP_fine;
    }
#endif
    
    
#ifdef PDAGP_2TEST
    
    // Test adding a second projector.
    {
        mgstruct.n_vector = 2; 
        complex<double>* tmp_store = mgstruct.null_vectors[0];
        delete[] mgstruct.null_vectors;
        mgstruct.null_vectors = new complex<double>*[mgstruct.n_vector];
        mgstruct.null_vectors[0] = tmp_store; 
        mgstruct.null_vectors[1] = new complex<double>[Lat.get_lattice_size()];

        // Add an even/odd vector. 
        for (int y = 0; y < y_fine; y++)
        {
            for (int x = 0; x < x_fine; x++)
            {
                if ((x+y)%2 == 0)
                    mgstruct.null_vectors[0][1][x+y*x_fine] = complex<double>(0.0,1.0);
                else
                    mgstruct.null_vectors[0][1][x+y*x_fine] = complex<double>(0.0,-1.0);
            }
        }

        cout << "Check projector:\n"; 
        for (int n = 0; n < mgstruct.n_vector; n++)
        {
            cout << "Vector " << n << "\n";
            for (int y = 0; y < y_fine; y++)
            {
                for (int x = 0; x < x_fine; x++)
                {
                    cout << mgstruct.null_vectors[n][x+y*x_fine] << " ";
                }
                cout << "\n";
            }
        }

        // Test block normalizing the projector.
        cout << "Block normalize the projector(s). X blocks: " << x_coarse << " Y blocks: " << y_coarse << "\n";
        block_normalize(&mgstruct); 

        cout << "Check block normalized vector:\n"; 
        for (int n = 0; n < mgstruct.n_vector; n++)
        {
            cout << "Vector " << n << "\n";
            for (int y = 0; y < y_fine; y++)
            {
                for (int x = 0; x < x_fine; x++)
                {
                    cout << mgstruct.null_vectors[0][n][x+y*x_fine] << " ";
                }
                cout << "\n";
            }
        }

        // Restrict the original source. 
        complex<double>* rhs_coarse = new complex<double>[coarse_size*mgstruct.n_vector];
        cout << "Restricting the source.\n";
        restrict(rhs_coarse, rhs, &mgstruct);

        // Check coarse source.
        cout << "Check coarse src:\n"; 
        for (int y = 0; y < y_coarse; y++)
        {
            for (int x = 0; x < x_coarse; x++)
            {
                cout << "(";
                for (int n = 0; n < mgstruct.n_vector; n++)
                {
                    cout << rhs_coarse[(x+y*x_coarse)*mgstruct.n_vector+n] << ",";
                }
                cout << ") ";
            }
            cout << "\n";
        }

        // Prolong the restricted source. 
        complex<double>* rhs_PdagP = new complex<double>[Lat.get_lattice_size()];
        cout << "Prolonging the restricted source.\n";
        prolong(rhs_PdagP, rhs_coarse, &mgstruct);

        cout << "Check re-prolonged source vector:\n"; 
        for (int y = 0; y < y_fine; y++)
        {
            for (int x = 0; x < x_fine; x++)
            {
                cout << rhs_PdagP[x+y*x_fine] << " ";
            }
            cout << "\n";
        }
        
        delete[] rhs_coarse;
        delete[] rhs_PdagP;
    }
    
#endif
    
#ifdef COARSE_ONLY
    cout << "[COARSE]: Solving coarse system only.\n";
    complex<double>* rhs_coarse = new complex<double>[coarse_size*mgstruct.n_vector];
    zero<double>(rhs_coarse, coarse_size*mgstruct.n_vector);
    restrict(rhs_coarse, rhs, &mgstruct);
    
    cout << "[COARSE]: Norm of coarse rhs " << sqrt(norm2sq<double>(rhs_coarse, coarse_size*mgstruct.n_vector)) << ".\n";
    
    complex<double>* lhs_coarse = new complex<double>[coarse_size*mgstruct.n_vector];
    zero<double>(lhs_coarse, coarse_size*mgstruct.n_vector);
    
    invif = minv_vector_gcr_restart(lhs_coarse, rhs_coarse, coarse_size*mgstruct.n_vector, 10000, 1e-6, 16, coarse_square_staggered, (void*)&mgstruct);
    
    
    if (invif.success == true)
    {
     printf("[COARSE]: Iterations %d RelRes %.8e Err N Algorithm %s\n", invif.iter, sqrt(invif.resSq)/bnorm, invif.name.c_str());
    }
    else // failed, maybe.
    {
     printf("[COARSE]: Iterations %d RelRes %.8e Err Y Algorithm %s\n", invif.iter, sqrt(invif.resSq)/bnorm, invif.name.c_str());
     printf("[COARSE]: This may be because the max iterations was reached.\n");
    }

    
    

    printf("[COARSE]: Computing [check] = A [lhs] as a confirmation.\n");

    // Check and make sure we get the right answer.
    complex<double>* A_lhs_coarse = new complex<double>[coarse_size*mgstruct.n_vector];
    zero<double>(A_lhs_coarse, coarse_size*mgstruct.n_vector);
    
    coarse_square_laplace(A_lhs_coarse, lhs_coarse, (void*)&mgstruct);

    for (i = 0; i < coarse_size*mgstruct.n_vector; i++)
    {
      explicit_resid += real(conj(A_lhs_coarse[i] - rhs_coarse[i])*(A_lhs_coarse[i] - rhs_coarse[i]));
    }
    explicit_resid = sqrt(explicit_resid);

    printf("[COARSE]: [check] should equal [rhs]. The residual is %15.20e.\n", explicit_resid);
    
    complex<double>* pro_lhs_coarse = new complex<double>[N*N];
    zero<double>(pro_lhs_coarse, N*N);
    prolong(pro_lhs_coarse, lhs_coarse, &mgstruct); 
    complex<double>* pro_rhs_coarse = new complex<double>[N*N];
    zero<double>(pro_rhs_coarse, N*N);
    square_staggered_u1(pro_rhs_coarse, pro_lhs_coarse, (void*)&stagif);
    
    // Compare PAP solution to real solution. 
    cout << "\n[COARSE]: Compare solutions.\n";
    double comparison = 0;
    double resid_comparison = 0;
    for (i = 0; i < Lat.get_lattice_size(); i++)
    {
        comparison += real(conj(pro_lhs_coarse[i]-lhs[i])*(pro_lhs_coarse[i]-lhs[i]));
        resid_comparison += real(conj(pro_rhs_coarse[i]-rhs[i])*(pro_rhs_coarse[i]-rhs[i]));
    }
    comparison = sqrt(explicit_resid);
    printf("[COARSE]: The solutions deviate by %15.20e.\n", comparison);
    printf("[COARSE]: The projected residual has a rel res of %15.20e.\n", sqrt(resid_comparison)/bnorm);
    
    delete[] rhs_coarse; 
    delete[] lhs_coarse;
    delete[] A_lhs_coarse; 
    delete[] pro_lhs_coarse; 
    delete[] pro_rhs_coarse; 
    
#endif // COARSE_ONLY
    
    // What operator are we using, D^\dag D or D?
    void (*fine_op)(complex<double>*,complex<double>*,void*);
    if (params.normal_eqn_mg)
    {
        fine_op = fine_square_staggered_normal;
    }
    else
    {
        fine_op = fine_square_staggered;
    }
    
    if (my_test == TOP_LEVEL_ONLY)
    {
        // Try a direct solve.
        cout << "\n[ORIG]: Solve fine system.\n";
        
        cout << "Mgstruct fine size " << mgstruct.curr_fine_size << "\n" << flush;
        
        switch (params.out_solve)
        {
            case OUTER_GCR:
                if (params.outer_restart)
                {
                    invif = minv_vector_gcr_restart(lhs, rhs, Lat.get_lattice_size(), params.outer_max_iter, params.outer_precision, params.outer_restart_freq, fine_op, (void*)&mgstruct, &verb);
                }
                else
                {
                    invif = minv_vector_gcr(lhs, rhs, Lat.get_lattice_size(), params.outer_max_iter, params.outer_precision, fine_op, (void*)&mgstruct, &verb);
                }
                break;
            case OUTER_CG:
                if (params.outer_restart)
                {
                    invif = minv_vector_cg_restart(lhs, rhs, Lat.get_lattice_size(), params.outer_max_iter, params.outer_precision, params.outer_restart_freq, fine_op, (void*)&mgstruct, &verb);
                }
                else
                {
                    invif = minv_vector_cg(lhs, rhs, Lat.get_lattice_size(), params.outer_max_iter, params.outer_precision, fine_op, (void*)&mgstruct, &verb);
                }
                break;
            case OUTER_BICGSTAB:
                if (params.outer_restart)
                {
                    invif = minv_vector_bicgstab_restart(lhs, rhs, Lat.get_lattice_size(), params.outer_max_iter, params.outer_precision, params.outer_restart_freq, fine_op, (void*)&mgstruct, &verb);
                }
                else
                {
                    invif = minv_vector_bicgstab(lhs, rhs, Lat.get_lattice_size(), params.outer_max_iter, params.outer_precision, fine_op, (void*)&mgstruct, &verb);
                }
                break;
        }
        
        //invif = minv_vector_gcr_restart(lhs, rhs, Lat.get_lattice_size(), 100000, outer_precision, outer_restart, square_staggered_u1, (void*)&stagif);
        mgstruct.dslash_count->krylov[mgstruct.curr_level] += (params.normal_eqn_mg ? 2 : 1)*invif.ops_count; 
        

        if (invif.success == true)
        {
         printf("[ORIG]: Iterations %d RelRes %.8e Err N Algorithm %s\n", invif.iter, sqrt(invif.resSq)/bnorm, invif.name.c_str());
        }
        else // failed, maybe.
        {
         printf("[ORIG]: Iterations %d RelRes %.8e Err Y Algorithm %s\n", invif.iter, sqrt(invif.resSq)/bnorm, invif.name.c_str());
         printf("[ORIG]: This may be because the max iterations was reached.\n");
        }

        printf("[ORIG]: Computing [check] = A [lhs] as a confirmation.\n");

        // Check and make sure we get the right answer.
        zero<double>(check, Lat.get_lattice_size());
        (*fine_op)(check, lhs, (void*)&mgstruct); 
        //square_staggered_u1(check, lhs, (void*)&stagif);

        explicit_resid = 0.0;
        for (i = 0; i < Lat.get_lattice_size(); i++)
        {
          explicit_resid += real(conj(rhs[i] - check[i])*(rhs[i] - check[i]));
        }
        explicit_resid = sqrt(explicit_resid)/bnorm;

        printf("[ORIG]: [check] should equal [rhs]. The relative residual is %15.20e.\n", explicit_resid);

    } // TOP_LEVEL_ONLY

    if (my_test == SMOOTHER_ONLY || my_test == TWO_LEVEL || my_test == THREE_LEVEL)
    {
        // Let's actually test a multigrid solve!
        cout << "\n[MG]: Test MG solve.\n";

        // Block normalize the null vectors.
        block_normalize(&mgstruct); 

        

        // Well, maybe this will work?
        zero<double>(lhs, Lat.get_lattice_size());
        
        switch (params.out_solve)
        {
            case OUTER_GCR:
                if (params.outer_restart)
                {
                    invif = minv_vector_gcr_var_precond_restart(lhs, rhs, Lat.get_lattice_size(), params.outer_max_iter, params.outer_precision, params.outer_restart_freq, fine_op, (void*)&mgstruct, mg_preconditioner, (void*)&mgprecond, &verb); 
                }
                else
                {
                    invif = minv_vector_gcr_var_precond(lhs, rhs, Lat.get_lattice_size(), params.outer_max_iter, params.outer_precision, fine_op, (void*)&mgstruct, mg_preconditioner, (void*)&mgprecond, &verb); 
                }
                break;
            case OUTER_CG:
                if (params.outer_restart)
                {
                    invif = minv_vector_cg_flex_precond_restart(lhs, rhs, Lat.get_lattice_size(), params.outer_max_iter, params.outer_precision, params.outer_restart_freq, fine_op, (void*)&mgstruct, mg_preconditioner, (void*)&mgprecond, &verb); 
                }
                else
                {
                    invif = minv_vector_cg_flex_precond(lhs, rhs, Lat.get_lattice_size(), params.outer_max_iter, params.outer_precision, fine_op, (void*)&mgstruct, mg_preconditioner, (void*)&mgprecond, &verb); 
                }
                break;
            case OUTER_BICGSTAB:
                if (params.outer_restart)
                {
                    invif = minv_vector_bicgstab_precond_restart(lhs, rhs, Lat.get_lattice_size(), params.outer_max_iter, params.outer_precision, params.outer_restart_freq, fine_op, (void*)&mgstruct, mg_preconditioner, (void*)&mgprecond, &verb); 
                }
                else
                {
                    invif = minv_vector_bicgstab_precond(lhs, rhs, Lat.get_lattice_size(), params.outer_max_iter, params.outer_precision, fine_op, (void*)&mgstruct, mg_preconditioner, (void*)&mgprecond, &verb); 
                }
                break;
        }
        
        mgstruct.dslash_count->krylov[mgstruct.curr_level] += (params.normal_eqn_mg ? 2 : 1)*invif.ops_count; 
        
        

        if (invif.success == true)
        {
            printf("[L1]: Iterations %d RelRes %.8e Err N Algorithm %s\n", invif.iter, sqrt(invif.resSq)/bnorm, invif.name.c_str());
        }
        else // failed, maybe.
        {
            printf("[L1]: Iterations %d RelRes %.8e Err Y Algorithm %s\n", invif.iter, sqrt(invif.resSq)/bnorm, invif.name.c_str());
            printf("[L1]: This may be because the max iterations was reached.\n");
        }


        printf("[MG]: Computing [check] = A [lhs] as a confirmation.\n");

        // Check and make sure we get the right answer.
        (*op)(check, lhs, (void*)&stagif);
        //square_staggered_u1(check, lhs, (void*)&stagif);

        explicit_resid = diffnorm2sq(rhs, check, Lat.get_lattice_size())/bnorm;

        printf("[MG]: [check] should equal [rhs]. The relative residual is %15.20e.\n", explicit_resid);
    } // SMOOTHER_ONLY or TWO_LEVEL
    
    // Look at the two point function!
    if (source == ORIGIN_POINT)
    {
        double* corr = new double[Lat.get_lattice_dimension(1)];
        cout << "BEGIN_GOLDSTONE\n";
        for (i = 0; i < Lat.get_lattice_dimension(1); i++)
        {
            corr[i] = 0.0;
            for (j = 0; j < Lat.get_lattice_dimension(0); j++)
            {
                corr[i] += real(conj(lhs[i*Lat.get_lattice_dimension(0)+j])*lhs[i*Lat.get_lattice_dimension(0)+j]);
            }
            cout << i << " " << corr[i] << "\n";
        }
        cout << "END_GOLDSTONE\n";
        cout << "BEGIN_EFFMASS\n";
        for (i = 0; i < Lat.get_lattice_dimension(1); i++)
        {
            cout << i << " " << log(corr[i]/corr[(i+1)%Lat.get_lattice_dimension(1)]) << "\n";
        }
        cout << "END_EFFMASS\n";
    }
    
    cout << "=================\n";
    cout << "= FINAL RESULTS =\n";
    cout << "=================\n";
    
    for (i = 0; i <= mgstruct.n_refine; i++)
    {
        cout << "= [L" << i+1 << "_DSLASH]: Dslash NullVec " << mgstruct.dslash_count->nullvectors[i] << " Krylov " << mgstruct.dslash_count->krylov[i] << " PreSmooth " << mgstruct.dslash_count->presmooth[i] <<
            " PostSmooth " << mgstruct.dslash_count->postsmooth[i] << " Residual " << mgstruct.dslash_count->residual[i] << " TotalNonNull " << invif.iter << " x " << 
            ((double)mgstruct.dslash_count->krylov[i]+mgstruct.dslash_count->presmooth[i]+mgstruct.dslash_count->postsmooth[i]+mgstruct.dslash_count->residual[i])/invif.iter << " = " << 
            (mgstruct.dslash_count->krylov[i]+mgstruct.dslash_count->presmooth[i]+mgstruct.dslash_count->postsmooth[i]+mgstruct.dslash_count->residual[i]) << "\n";
    }
    
    // Free the lattice.
    delete[] lattice;
    delete[] lhs;
    delete[] rhs;
    delete[] check;
    delete[] tmp;
    delete[] tmp2;
    delete[] gauge_trans;
    
    // Clean up!
    delete[] coord;
    delete[] mgstruct.blocksize_x;
    delete[] mgstruct.blocksize_y; 
    for (i = 0; i < mgstruct.n_refine; i++)
    {
        for (j = 0; j < mgstruct.n_vector; j++)
        {
            delete[] mgstruct.null_vectors[i][j];
        }
        delete[] mgstruct.null_vectors[i];
    }
    delete mgstruct.dslash_count; 
    delete[] mgstruct.null_vectors; 
    if (mgstruct.n_refine > 0)
    {
        for (i = 0; i <= mgstruct.n_refine; i++)
        {
            delete mgstruct.latt[i];
            if (mgstruct.stencils[i] != 0) { delete mgstruct.stencils[i]; }
            if (mgstruct.have_dagger_stencil && mgstruct.dagger_stencils[i] != 0) { delete mgstruct.dagger_stencils[i]; }
        }
    }
    delete[] mgstruct.latt;
    delete[] mgstruct.stencils; 
    if (mgstruct.have_dagger_stencil)
    {
        delete mgstruct.dagger_stencils;
    }
    
    delete[] mgprecond.n_pre_smooth;
    delete[] mgprecond.n_post_smooth; 
    delete[] mgprecond.rel_res; 
    
    return 0; 
}


// Custom routine to load gauge field.
void internal_load_gauge_u1(complex<double>* lattice, int x_fine, int y_fine, double BETA)
{
    if (x_fine == 32 && y_fine == 32)
    {
        if (abs(BETA - 3.0) < 1e-8)
        {
            read_gauge_u1(lattice, x_fine, y_fine, "./cfg/l32t32b30_heatbath.dat");
            cout << "[GAUGE]: Loaded a U(1) gauge field with non-compact beta = " << 1.0/sqrt(BETA) << "\n";
        }
        else if (abs(BETA - 6.0) < 1e-8)
        {
            read_gauge_u1(lattice, x_fine, y_fine, "./cfg/l32t32b60_heatbath.dat");
            cout << "[GAUGE]: Loaded a U(1) gauge field with non-compact beta = " << 1.0/sqrt(BETA) << "\n";
        }
        else if (abs(BETA - 10.0) < 1e-8)
        {
            read_gauge_u1(lattice, x_fine, y_fine, "./cfg/l32t32b100_heatbath.dat");
            cout << "[GAUGE]: Loaded a U(1) gauge field with non-compact beta = " << 1.0/sqrt(BETA) << "\n";
        }
        else if (abs(BETA - 10000.0) < 1e-8)
        {
            read_gauge_u1(lattice, x_fine, y_fine, "./cfg/l32t32bperturb_heatbath.dat");
            cout << "[GAUGE]: Loaded a U(1) gauge field with non-compact beta = " << 1.0/sqrt(BETA) << "\n";
        }
        else
        {
            cout << "[GAUGE]: Saved U(1) gauge field with correct beta, volume does not exist.\n";
        }
    }
    else if (x_fine == 64 && y_fine == 64)
    {
        if (abs(BETA - 3.0) < 1e-8)
        {
            read_gauge_u1(lattice, x_fine, y_fine, "./cfg/l64t64b30_heatbath.dat");
            cout << "[GAUGE]: Loaded a U(1) gauge field with non-compact beta = " << 1.0/sqrt(BETA) << "\n";
        }
        else if (abs(BETA - 6.0) < 1e-8)
        {
            read_gauge_u1(lattice, x_fine, y_fine, "./cfg/l64t64b60_heatbath.dat");
            cout << "[GAUGE]: Loaded a U(1) gauge field with non-compact beta = " << 1.0/sqrt(BETA) << "\n";
        }
        else if (abs(BETA - 10.0) < 1e-8)
        {
            read_gauge_u1(lattice, x_fine, y_fine, "./cfg/l64t64b100_heatbath.dat");
            cout << "[GAUGE]: Loaded a U(1) gauge field with non-compact beta = " << 1.0/sqrt(BETA) << "\n";
        }
        else if (abs(BETA - 10000.0) < 1e-8)
        {
            read_gauge_u1(lattice, x_fine, y_fine, "./cfg/l64t64bperturb_heatbath.dat");
            cout << "[GAUGE]: Loaded a U(1) gauge field with non-compact beta = " << 1.0/sqrt(BETA) << "\n";
        }
        else
        {
            cout << "[GAUGE]: Saved U(1) gauge field with correct beta, volume does not exist.\n";
        }
    }
    else if (x_fine == 128 && y_fine == 128)
    {
        if (abs(BETA - 6.0) < 1e-8)
        {
            read_gauge_u1(lattice, x_fine, y_fine, "./cfg/l128t128b60_heatbath.dat");
            cout << "[GAUGE]: Loaded a U(1) gauge field with non-compact beta = " << 1.0/sqrt(BETA) << "\n";
        }
        else if (abs(BETA - 10.0) < 1e-8)
        {
            read_gauge_u1(lattice, x_fine, y_fine, "./cfg/l128t128b100_heatbath.dat");
            cout << "[GAUGE]: Loaded a U(1) gauge field with non-compact beta = " << 1.0/sqrt(BETA) << "\n";
        }
        else if (abs(BETA - 10000.0) < 1e-8)
        {
            read_gauge_u1(lattice, x_fine, y_fine, "./cfg/l128t128bperturb_heatbath.dat");
            cout << "[GAUGE]: Loaded a U(1) gauge field with non-compact beta = " << 1.0/sqrt(BETA) << "\n";
        }
        else
        {
            cout << "[GAUGE]: Saved U(1) gauge field with correct beta, volume does not exist.\n";
        }
    }
    else
    {
        cout << "[GAUGE]: Saved U(1) gauge field with correct beta, volume does not exist. Using unit gauge.\n";
    }
}



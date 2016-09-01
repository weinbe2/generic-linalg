// This is a sample code that constructs general kinetic terms!

#include <iostream>
#include <iomanip> // to set output precision.
#include <cmath>
#include <string>
#include <sstream>
#include <complex>
#include <random>
#include <cstring> // should be replaced by using sstream


#include "generic_vector.h"
#include "verbosity.h"

#include "generic_bicgstab.h"
#include "generic_cg.h"
#include "generic_cr.h"
#include "generic_gcr.h"
#include "generic_minres.h"
#include "generic_cg_flex_precond.h"
#include "generic_gcr_var_precond.h"

#include "mg.h"
#include "mg_complex.h"
#include "null_gen.h"
#include "u1_utils.h"
#include "lattice.h"
#include "operators.h"

// Are we checking eigenvalues?
#define EIGEN_TEST

#ifdef EIGEN_TEST
#include "arpack_interface.h"
#endif

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

// Do null vector generation? Currently uses BiCGStab
#define GEN_NULL_VECTOR



// What type of test should we do?
enum mg_test_types
{
    TOP_LEVEL_ONLY = 0,       // only test the top level solver.
    SMOOTHER_ONLY = 1,        // Top level + smoother
    TWO_LEVEL = 2,            // Two level MG 
    THREE_LEVEL = 3           // Three level MG
};



// What gauge field do we use? Load, random, unit?
enum gauge_create_type
{
    GAUGE_LOAD = 0,             // Load a gauge field.
    GAUGE_RANDOM = 1,           // Create a gauge field with deviation 1/sqrt(beta)
    GAUGE_UNIT = 2              // Use a unit gauge field.
};

// Custom routine to load gauge field.
void internal_load_gauge_u1(complex<double>* lattice, int x_fine, int y_fine, double BETA);


enum op_type
{
    STAGGERED = 0,
    LAPLACE = 1,
    LAPLACE_NC2 = 2,
    G5_STAGGERED = 3,
    STAGGERED_NORMAL = 4,
    STAGGERED_INDEX = 5
};

enum src_type
{
    POINT = 0,
    RANDOM_GAUSSIAN = 1,
    ORIGIN_POINT = 2 // for correlator test. 
};

/*struct staggered_u1_op
{
    complex<double> *lattice;
    double mass;
    int x_fine;
    int y_fine; 
    Lattice* Lat; 
    int Nc; // only relevant for square laplace. 
};*/

// Print usage. 
void display_usage()
{
    cout << "--help\n";
    cout << "--lattice-size [32, 64, 128] {##}             (default 32x32)\n";
    cout << "--operator [laplace, laplace2, staggered\n";
    cout << "       g5_staggered, normal_staggered, index] (default staggered)\n";
    cout << "--null-operator [laplace, laplace2, staggered\n";
    cout << "       g5_staggered, normal_staggered, index] (default staggered)\n";
    cout << "--null-solver [gcr, bicgstab, cg, minres]     (default bicgstab)\n";
    cout << "--null-precision [null prec]                  (default 5e-5)\n";
    cout << "--null-mass [mass]                            (default 1e-4)\n";
    cout << "--null-eo [corner, yes, no, topo]             (default yes)\n";
    cout << "--null-global-ortho-conj [yes, no]            (default no, it only helps in some weird fluke cases)\n";
    cout << "--null-ortho-eo [yes, no]                     (default no)\n"; 
    cout << "--mass [mass]                                 (default 1e-2)\n";
    cout << "--blocksize [blocksize]                       (default 4)\n";
    cout << "--nvec [nvec]                                 (default 4)\n";
    cout << "--nrefine [number coarse]                     (default 1)\n";
    cout << "--multi-strategy [smooth, recursive]          (default smooth)\n";
    cout << "--coarse-precision [coarse precision]         (default 1e-3)\n";
    cout << "--coarse-max-iter [coarse maximum iterations] (default 1024)\n";
    cout << "--gauge [unit, load, random]                  (default load)\n";
    cout << "--gauge-transform [yes, no]                   (default no)\n";
    cout << "--beta [3.0, 6.0, 10.0, 10000.0]              (default 6.0)\n";
    cout << "--npre-smooth [presmooth steps]               (default 6)\n";
    cout << "--npost-smooth [postsmooth steps]             (default 6)\n";
    cout << "--load-cfg [path]                             (default do not load, overrides beta)\n";
#ifdef EIGEN_TEST
    cout << "--do-eigentest [yes, no]                      (default no)\n";
    cout << "--n-ev [all, # smallest]                      (default all)\n";
    cout << "--n-cv [#]                                    (default min(all, 2.5*n-ev))\n";
#endif // EIGEN_TEST
}

int main(int argc, char** argv)
{  
    // Declare some variables.
    cout << setiosflags(ios::scientific) << setprecision(6);
    int i, j, k, x, y;
    complex<double> *lattice; // Holds the gauge field.
    complex<double> *lhs, *rhs, *check; // For some Kinetic terms.
    complex<double> *tmp, *tmp2; // For temporary space. 
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
    int color = 0;
    
    // Set parameters. 
    
    // How are we creating the gauge field? Load it, random, unit? Can set on command line.
    gauge_create_type gauge_load = GAUGE_LOAD; // GAUGE_LOAD, GAUGE_UNIT, GAUGE_RANDOM
    
    // Should we do a random gauge rotation? Can set on command line.
    bool do_gauge_transform = false; 
    
    // What operator are we using for the solve? (Laplace is free only.) Can set on command line.
    op_type opt = STAGGERED; // STAGGERED, LAPLACE, LAPLACE_NC2, G5_STAGGERED
    
    // What operator are we using for null vector generation? Can set on command line.
    op_type opt_null = STAGGERED;

    // What test are we performing?
    mg_test_types my_test = TWO_LEVEL; //THREE_LEVEL; // TWO_LEVEL is the default which won't override anything.
    
    // How are we generating null vectors?
    mg_null_gen_type null_gen = NULL_BICGSTAB; // NULL_BICGSTAB, NULL_GCR, NULL_CG, NULL_MINRES
    
    // L_x = L_y = Dimension for a lattice.
    int lattice_size_x = 32; // Can be set on command line with --lattice-size. 
    int lattice_size_y = 32; 
    
    // Describe the staggered fermions.
    double MASS = 0.01; // Can be overridden on command line with --mass 
    
    // Describe the source type.
    src_type source = RANDOM_GAUSSIAN; // POINT, RANDOM_GAUSSIAN, or ORIGIN_POINT (for correlator test).
    
    // Outer Inverter information.
    double outer_precision = 5e-7; 
    int outer_restart = 64; 
    
    // Multigrid information. 
    int n_refine = 1; // 1 = two level V cycle, 2 = three level V cycle, etc. 
                      // Can be set on command line with --nrefine
    if (my_test == THREE_LEVEL) // FOR TEST ONLY
    {
        n_refine = 2;
    }
    int X_BLOCKSIZE = 4; // Can be overrided with below arg on command line
    int Y_BLOCKSIZE = 4; // with --blocksize 
    blocking_strategy bstrat = BLOCK_EO; // BLOCK_NONE, BLOCK_EO, BLOCK_CORNER
    int null_partitions = 2; 
    
    // Null vector generation
    
// If GEN_NULL_VECTOR isn't defined:
//   1 for just const vector, 2 for const + even/odd vector, 4 for each corner
//    of the hypercube.
// If GEN_NULL_VECTOR is defined and eo = 0:
//   Generate "n_null_vector" null vectors which are block orthogonalized.
// IF GEN_NULL_VECTOR is defined and eo = 1:
//   Generate "n_null_vector" null vectors, partition into even and odd.
//    Total number of null vectors is 2*VECTOR_COUNT. 
// IF GEN_NULL_VECTOR is defined and eo = 3:
//   Generate "n_null_vector" null vectors, partition into four corners.
//    Total number of null vectors is 4*VECTOR_COUNT. 
    int n_null_vector = 4; // Note: Gets multiplied by 2 for LAPLACE_NC2 test.
                           // Can be overriden on command line with --nvec
    int null_max_iter = 500;
    double null_precision = 5e-5; //5e-4; // Can be overriden on command line with --null-precision
    double null_mass = 1e-4; // What mass do we use for null vector generation?
    
    // Do we globally orthogonalize null vectors both against previous null vectors and their conjugate?
    bool do_global_ortho_conj = false;
    
    // Do we split null vectors into even/odd, then orthogonalize, or do we orthogonalize first?
    bool do_ortho_eo = false; 
    
    // Inner solver.
    mg_multilevel_type mlevel_type = MLEVEL_SMOOTH; // MLEVEL_SMOOTH, MLEVEL_RECURSIVE --- do we smooth then go down, or smooth then krylov?
    inner_solver in_solve = GCR; //CR; //GCR; 
    double inner_precision = 1e-3;
    int inner_restart = 64;
    int inner_max = 1024;
    if (my_test == SMOOTHER_ONLY)
    {
        in_solve = NONE; 
    }
    
    // Smoother
    inner_solver in_smooth = GCR; //NONE; //GCR; BICGSTAB
    double omega_smooth = 0.67; // for MINRES only. 
    int pre_smooth = 6; // Can set on command line.
    int post_smooth = 6; // Can set on command line.
    
    // Gauge field information.
    double BETA = 6.0; // For random gauge field, phase angles have std.dev. 1/sqrt(beta).
                       // For heatbath gauge field, corresponds to non-compact beta.
                       // Can be set on command line with --beta.
    
    // Load an external cfg?
    char* load_cfg = NULL;
    bool do_load = false; 
    
#ifdef EIGEN_TEST
    bool do_eigentest = false;
    
    int set_eigen = -1; // default for generating all eigenvalues.
    int set_cv = -1; // default for generating all eigenvalues, or 
#endif // EIGEN_TEST 
    
    /////////////////////////////////////////////
    // Get a few parameters from command line. //
    /////////////////////////////////////////////
    for (i = 1; i < argc; i++)
    {
        if (strcmp(argv[i], "--help") == 0)
        {
            display_usage();
            return 0;
        }
        if (i+1 != argc)
        {
            if (strcmp(argv[i], "--mass") == 0)
            {
                MASS = atof(argv[i+1]);
                i++;
            }
            else if (strcmp(argv[i], "--operator") == 0)
            {
                if (strcmp(argv[i+1], "laplace") == 0)
                {
                    opt = LAPLACE; 
                }
                else if (strcmp(argv[i+1], "laplace2") == 0)
                {
                    opt = LAPLACE_NC2; 
                }
                else if (strcmp(argv[i+1], "staggered") == 0)
                {
                    opt = STAGGERED; 
                }
                else if (strcmp(argv[i+1], "g5_staggered") == 0)
                {
                    opt = G5_STAGGERED; 
                }
                else if (strcmp(argv[i+1], "normal_staggered") == 0)
                {
                    opt = STAGGERED_NORMAL;
                }
                else if (strcmp(argv[i+1], "index") == 0)
                {
                    opt = STAGGERED_INDEX;
                }
                i++;
            }
            else if (strcmp(argv[i], "--null-operator") == 0)
            {
                if (strcmp(argv[i+1], "laplace") == 0)
                {
                    opt_null = LAPLACE; 
                }
                else if (strcmp(argv[i+1], "laplace2") == 0)
                {
                    opt_null = LAPLACE_NC2; 
                }
                else if (strcmp(argv[i+1], "staggered") == 0)
                {
                    opt_null = STAGGERED; 
                }
                else if (strcmp(argv[i+1], "g5_staggered") == 0)
                {
                    opt_null = G5_STAGGERED; 
                }
                else if (strcmp(argv[i+1], "normal_staggered") == 0)
                {
                    opt_null = STAGGERED_NORMAL;
                }
                else if (strcmp(argv[i+1], "index") == 0)
                {
                    opt_null = STAGGERED_INDEX;
                }
                i++;
            }
            else if (strcmp(argv[i], "--blocksize") == 0)
            {
                X_BLOCKSIZE = atoi(argv[i+1]);
                Y_BLOCKSIZE = X_BLOCKSIZE;
                i++;
            }
            else if (strcmp(argv[i], "--nvec") == 0)
            {
                n_null_vector = atoi(argv[i+1]);
                i++;
            }
            else if (strcmp(argv[i], "--null-precision") == 0)
            {
                null_precision = atof(argv[i+1]);
                i++;
            }
            else if (strcmp(argv[i], "--null-mass") == 0)
            {
                null_mass = atof(argv[i+1]);
                i++;
            }
            else if (strcmp(argv[i], "--null-solver") == 0)
            {
                if (strcmp(argv[i+1], "gcr") == 0)
                {
                    null_gen = NULL_GCR;
                }
                else if (strcmp(argv[i+1], "bicgstab") == 0)
                {
                    null_gen = NULL_BICGSTAB;
                }
                else if (strcmp(argv[i+1], "cg") == 0)
                {
                    null_gen = NULL_CG;
                }
                else if (strcmp(argv[i+1], "minres") == 0)
                {
                    null_gen = NULL_MINRES;
                }
                i++;
            }
            else if (strcmp(argv[i], "--null-eo") == 0)
            {
                if (strcmp(argv[i+1], "yes") == 0)
                {
                    bstrat = BLOCK_EO; // even/odd
                }
                else if (strcmp(argv[i+1], "corner") == 0)
                {
                    bstrat = BLOCK_CORNER; // corners
                }
                else if (strcmp(argv[i+1], "topo") == 0)
                {
                    bstrat = BLOCK_TOPO; // chirality as defined by taste singlet. 
                }
                else // none.
                {
                    bstrat = BLOCK_NONE; // fully reduce. 
                }
                i++;
            }
            else if (strcmp(argv[i], "--null-global-ortho-conj") == 0)
            {
                if (strcmp(argv[i+1], "yes") == 0)
                {
                    do_global_ortho_conj = true; // yes, globally orthogonalize null vectors against previous and conj.
                }
                else if (strcmp(argv[i+1], "no") == 0)
                {
                    do_global_ortho_conj = false;
                }
                i++;
            }
            else if (strcmp(argv[i], "--null-ortho-eo") == 0)
            {
                if (strcmp(argv[i+1], "yes") == 0)
                {
                    do_ortho_eo = true; // yes, split null vectors into eo before orthogonalizing.
                }
                else if (strcmp(argv[i+1], "no") == 0)
                {
                    do_ortho_eo = false; 
                }
                i++;
            }
            else if (strcmp(argv[i], "--nrefine") == 0)
            {
                n_refine = atoi(argv[i+1]);
                i++;
            }
            else if (strcmp(argv[i], "--multi-strategy") == 0)
            {
                if (strcmp(argv[i+1], "smooth") == 0)
                {
                    mlevel_type = MLEVEL_SMOOTH;
                }
                else if (strcmp(argv[i+1], "recursive") == 0)
                {
                    mlevel_type = MLEVEL_RECURSIVE; 
                }
                i++;
            }
            else if (strcmp(argv[i], "--coarse-precision") == 0)
            {
                inner_precision = atof(argv[i+1]);
                i++;
            }
            else if (strcmp(argv[i], "--coarse-max-iter") == 0)
            {
                inner_max = atoi(argv[i+1]);
                i++;
            }
            else if (strcmp(argv[i], "--gauge") == 0)
            {
                if (strcmp(argv[i+1], "unit") == 0)
                {
                    gauge_load = GAUGE_UNIT;
                }
                else if (strcmp(argv[i+1], "random") == 0)
                {
                    gauge_load = GAUGE_RANDOM;
                }
                else
                {
                    gauge_load = GAUGE_LOAD;
                }
                i++;
            }
            else if (strcmp(argv[i], "--gauge-transform") == 0)
            {
                if (strcmp(argv[i+1], "yes") == 0)
                {
                    do_gauge_transform = true;
                }
                else
                {
                    do_gauge_transform = false;
                }
                i++;
            }
            else if (strcmp(argv[i], "--beta") == 0)
            {
                BETA = atof(argv[i+1]);
                i++;
            }
            else if (strcmp(argv[i], "--lattice-size") == 0)
            {
                lattice_size_x = atoi(argv[i+1]);
                
                if (i+2 != argc)
                {
                    if (argv[i+2][0] == '-' && argv[i+2][1] == '-') // look for --
                    {
                        lattice_size_y = lattice_size_x;
                    }
                    else // assume number
                    {
                        lattice_size_y = atoi(argv[i+2]);
                        i++;
                    }
                }
                else
                {
                    lattice_size_y = lattice_size_x; // At the end, don't try to grab the next element!
                }
                i++;
            }
            else if (strcmp(argv[i], "--npre-smooth") == 0)
            {
                pre_smooth = atof(argv[i+1]);
                i++;
            }
            else if (strcmp(argv[i], "--npost-smooth") == 0)
            {
                post_smooth = atof(argv[i+1]);
                i++;
            }
            else if (strcmp(argv[i], "--load-cfg") == 0)
            {
                load_cfg = argv[i+1];
                do_load = true;
                i++;
            }
#ifdef EIGEN_TEST
            else if (strcmp(argv[i], "--do-eigentest") == 0)
            {
                if (strcmp(argv[i+1], "yes") == 0)
                {
                    do_eigentest = true;
                }
                else
                {
                    do_eigentest = false;
                }
                i++;
            }
            else if (strcmp(argv[i], "--n-ev") == 0)
            {
                if (strcmp(argv[i+1], "all") == 0)
                {
                    set_eigen = -1;
                    set_cv = -1;
                }
                else
                {
                    set_eigen = atoi(argv[i+1]);
                }
                i++;
            }
            else if (strcmp(argv[i], "--n-cv") == 0)
            {
                set_cv = atoi(argv[i+1]);
                i++;
            }
#endif // EIGENTEST
            else
            {
                display_usage();
                return 0;
            }
        }
    }
    
    //printf("Mass %.8e Blocksize %d %d Null Vectors %d\n", MASS, X_BLOCKSIZE, Y_BLOCKSIZE, n_null_vector);
    //return 0;
    
    ///////////////////////////////////////
    // End of human-readable parameters! //
    ///////////////////////////////////////
    
    
    string op_name;
    void (*op)(complex<double>*, complex<double>*, void*);
    switch (opt)
    {
        case STAGGERED:
            op = square_staggered_u1;
            op_name = "Staggered U(1)";
            break;
        case LAPLACE:
            op_name = "Free Laplace";
            op = square_laplace;
            break;
        case LAPLACE_NC2:
            op_name = "Free Laplace Nc = 2";
            op = square_laplace;
            break;
        case G5_STAGGERED:
            op_name = "Gamma_5 Staggered U(1)";
            op = square_staggered_gamma5_u1;
            break;
        case STAGGERED_NORMAL:
            op_name = "Staggered U(1) Normal";
            op = square_staggered_normal_u1;
            break;
        case STAGGERED_INDEX:
            op_name = "Staggered U(1) Index Operator";
            op = staggered_index_operator;
            break;
    }
    cout << "[OP]: Operator " << op_name << " Mass " << MASS << "\n";
    
    // Only relevant for free laplace test.
    int Nc = 1;  // Only value that matters for staggered
    if (opt == LAPLACE_NC2)
    {
        Nc = 2;
    }
    
    // Unset eo for Laplace. 
    if (opt == LAPLACE || opt == LAPLACE_NC2) // FOR TEST ONLY
    {
        bstrat = BLOCK_NONE;
    }
    
    // Describe the fine lattice. 
    lattice_size[0] = lattice_size_x;
    lattice_size[1] = lattice_size_y;
    
    // Create a lattice object.
    Lattice Lat(nd, lattice_size, Nc);
    
    cout << "[VOL]: X " << lattice_size[0] << " Y " << lattice_size[1] << " Volume " << Lat.get_volume();
    if (opt == LAPLACE || opt == LAPLACE_NC2) // FOR TEST ONLY
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
    stagif.mass = MASS; 
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
    cout << "[GAUGE]: Creating a gauge field.\n";
    
    switch (gauge_load)
    {
        case GAUGE_UNIT:
            unit_gauge_u1(lattice, Lat.get_lattice_dimension(0), Lat.get_lattice_dimension(1));
            break;
        case GAUGE_RANDOM:
            gauss_gauge_u1(lattice, Lat.get_lattice_dimension(0), Lat.get_lattice_dimension(1), generator, BETA);
            cout << "[GAUGE]: Created a U(1) gauge field with angle standard deviation " << 1.0/sqrt(BETA) << "\n";
            break;
        case GAUGE_LOAD:
            // Unit first in case loading fails.
            unit_gauge_u1(lattice, Lat.get_lattice_dimension(0), Lat.get_lattice_dimension(1));
            // Load the gauge field.
            if (do_load)
            {
                read_gauge_u1(lattice, Lat.get_lattice_dimension(0), Lat.get_lattice_dimension(1), load_cfg);
                cout << "[GAUGE]: Loaded a U(1) gauge field from " << load_cfg << "\n";
            }
            else // various predefined cfgs. 
            {
                internal_load_gauge_u1(lattice, Lat.get_lattice_dimension(0), Lat.get_lattice_dimension(1), BETA); // defined at end of file.
            }
            break;
    }
    
    if (do_gauge_transform)
    {
        // Generate and perform a random gauge transformation.
        rand_trans_u1(gauge_trans, Lat.get_lattice_dimension(0), Lat.get_lattice_dimension(1), generator);
        apply_gauge_trans_u1(lattice, gauge_trans, Lat.get_lattice_dimension(0), Lat.get_lattice_dimension(1));
        cout << "[GAUGE]: Performed a random gauge rotation.\n";
    }
    
    cout << "[GAUGE]: The average plaquette is " << get_plaquette_u1(lattice, Lat.get_lattice_dimension(0), Lat.get_lattice_dimension(1)) << ".\n";
    
    // Sanity check if we're doing a multigrid solve.
    if (n_refine == 0)
    {
        my_test = TOP_LEVEL_ONLY;
    }
    
    string op_null_name;
    void (*op_null)(complex<double>*, complex<double>*, void*) = square_staggered_u1; 
    switch (opt_null)
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
    mgstruct.n_refine = n_refine; 
    mgstruct.blocksize_x = new int[n_refine];
    mgstruct.blocksize_y = new int[n_refine];
    //mgstruct.blocksize_x[0] = 8;
    //mgstruct.blocksize_y[0] = 8;
    //for (i = 1; i < n_refine; i++)
    for (i = 0; i < n_refine; i++)
    {
        mgstruct.blocksize_x[i] = X_BLOCKSIZE;
        mgstruct.blocksize_y[i] = Y_BLOCKSIZE;
    }

    switch (bstrat)
    {
        case BLOCK_NONE:
            mgstruct.n_vector = n_null_vector;
            null_partitions = 1;
            break;
        case BLOCK_EO:
        case BLOCK_TOPO:
            mgstruct.n_vector = 2*n_null_vector;
            null_partitions = 2;
            break;
        case BLOCK_CORNER:
            mgstruct.n_vector = 4*n_null_vector;
            null_partitions = 4; 
            break;
    }
    mgstruct.matrix_vector = op; //square_staggered_u1;
    mgstruct.matrix_extra_data = (void*)&stagif; 
    
    cout << "[MG]: X_Block " << X_BLOCKSIZE << " Y_Block " << Y_BLOCKSIZE << " NullVectors " << n_null_vector << "\n";
    
    // Set the starting mg_struct state.
    mgstruct.curr_level = 0; // Ready to do top level -> second level.
    mgstruct.curr_dof_fine = Nc; // Top level has only one d.o.f. per site. 
    mgstruct.curr_x_fine = mgstruct.x_fine;
    mgstruct.curr_y_fine = mgstruct.y_fine;
    mgstruct.curr_fine_size = mgstruct.curr_y_fine*mgstruct.curr_x_fine*mgstruct.curr_dof_fine;
    
    mgstruct.curr_dof_coarse = mgstruct.n_vector; 
    mgstruct.curr_x_coarse = mgstruct.x_fine/mgstruct.blocksize_x[0];
    mgstruct.curr_y_coarse = mgstruct.y_fine/mgstruct.blocksize_y[0];
    mgstruct.curr_coarse_size = mgstruct.curr_y_coarse*mgstruct.curr_x_coarse*mgstruct.curr_dof_coarse;
    
    
    // Build the mg inverter structure.
    // Set up the MG preconditioner. 
    mg_precond_struct_complex mgprecond;

    mgprecond.in_smooth_type = in_smooth; // What inner smoother? MinRes or GCR.
    mgprecond.omega_smooth = omega_smooth; // What relaxation parameter should we use (MinRes only!)
    mgprecond.n_pre_smooth = pre_smooth; // 6 MinRes smoother steps before coarsening.
    mgprecond.n_post_smooth = post_smooth; // 6 MinRes smoother steps after refining.
    mgprecond.mlevel_type = mlevel_type; // Do we smooth then go down, or smooth then start a new Krylov?
    mgprecond.in_solve_type = in_solve; // What inner solver? NONE, MINRES, CG, GCR, BICGSTAB
    mgprecond.n_max = inner_max; // max number of steps to use for inner solver.
    mgprecond.n_restart = inner_restart; // frequency of restart (relevant for CG, GCR).
    mgprecond.rel_res = inner_precision; // Maximum relative residual for inner solver.
    mgprecond.mgstruct = &mgstruct; // Contains null vectors, fine operator. (Since we don't construct the fine op.)
    mgprecond.coarse_matrix_vector = coarse_square_staggered; // Function which applies the coarse operator. 
    mgprecond.fine_matrix_vector = fine_square_staggered; // Function which applies the fine operator. 
    mgprecond.matrix_extra_data = (void*)&mgstruct; // What extra_data the coarse operator expects. 

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
    
    // Create a projector.
    mgstruct.null_vectors = new complex<double>**[mgstruct.n_refine];
    // The top level is special since there are no color indices.
    mgstruct.null_vectors[0] = new complex<double>*[mgstruct.n_vector];
    for (j = 0; j < mgstruct.n_vector; j++)
    {
        mgstruct.null_vectors[0][j] = new complex<double>[Lat.get_lattice_size()];
        zero<double>(mgstruct.null_vectors[0][j], Lat.get_lattice_size());
    }
    // Higher levels are different.
    if (mgstruct.n_refine > 1)
    {
        for (i = 1; i < mgstruct.n_refine; i++)
        {
            level_down(&mgstruct);

            mgstruct.null_vectors[i] = new complex<double>*[mgstruct.n_vector];
            for (j = 0; j < mgstruct.n_vector; j++)
            {
                mgstruct.null_vectors[i][j] = new complex<double>[mgstruct.curr_x_fine*mgstruct.curr_y_fine*mgstruct.curr_dof_fine];
                zero<double>(mgstruct.null_vectors[i][j], mgstruct.curr_x_fine*mgstruct.curr_y_fine*mgstruct.curr_dof_fine);
            }
        }
        // Come back up!
        for (i = mgstruct.n_refine; i > 1; i--)
        {
            level_up(&mgstruct);
        }
    }

    cout << "[MG]: Creating " << mgstruct.n_vector << " projector(s).\n";
#ifndef GEN_NULL_VECTOR

    // Make a constant projector.
    if (mgstruct.n_vector == 1)
    {
        cout << "[MG]: Null vector 1 is a constant.\n";
        for (i = 0; i < Lat.get_lattice_size(); i++)
        {
            mgstruct.null_vectors[0][0][i] = 1;
            if (do_gauge_transform) 
            {
                mgstruct.null_vectors[0][0][i] *= gauge_trans[i];
            }
        }
    }
    else if (mgstruct.n_vector == 2) // constant, even/odd phase. 
    {
        cout << "[MG]: Null vector 1 is a constant.\n";
        cout << "[MG]: Null vector 2 is an even/odd phase.\n";
        for (i = 0; i < Lat.get_lattice_size(); i++)
        {
            mgstruct.null_vectors[0][0][i] = 1;
            x = i % x_fine;
            y = i / x_fine;
            mgstruct.null_vectors[0][1][i] = ((x+y)%2 == 0) ? complex<double>(0.0,1.0) : complex<double>(0.0,-1.0);
            if (do_gauge_transform)
            {
                mgstruct.null_vectors[0][0][i] *= (gauge_trans[i]);
                mgstruct.null_vectors[0][1][i] *= (gauge_trans[i]);
            }
        }
    }
    else if (mgstruct.n_vector == 4) // 4 corners of hypercube.
    {
        cout << "[MG]: Null vector 1 is a constant on unit corner (0,0).\n";
        cout << "[MG]: Null vector 2 is a constant on unit corner (1,0).\n";
        cout << "[MG]: Null vector 3 is a constant on unit corner (0,1).\n";
        cout << "[MG]: Null vector 4 is a constant on unit corner (1,1).\n";
        // Generate a normal distribution.
        std::normal_distribution<> dist(1.0, 0.1);
        for (i = 0; i < Lat.get_lattice_size(); i++)
        {
            x = i % x_fine;
            y = i / x_fine;
            mgstruct.null_vectors[0][2*(y%2)+(x%2)][i] = 1.0;
            //mgstruct.null_vectors[0][2*(y%2)+(x%2)][i] = dist(generator);
            
            if (do_gauge_transform)
            {
                mgstruct.null_vectors[0][2*(y%2)+(x%2)][i] *= (gauge_trans[i]);
            }
        }
    }
    else // invalid.
    {
        cout << "Unless you are generating null vectors, you can only use 1, 2, or 4 null vectors.\n";
        return 0;
    }
    
    cout << "[MG]: Performing block orthonormalize of null vectors...\n";
    block_orthonormalize(&mgstruct); 
    
#ifdef PRINT_NULL_VECTOR
    cout << "[MG]: Check projector:\n"; 
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
    
        

#else // generate null vector
    
    // Skip this depending on our test!

    if (!(my_test == TOP_LEVEL_ONLY || my_test == SMOOTHER_ONLY))
    {
        // Generate top level!
        cout << "[NULLVEC]: About to generate null vectors.\n";
        
        mgstruct.matrix_vector = op_null; // Trick op to null gen op.
        
        // Temporarily set the mass to 1e-4 for the null vector generation. 
        stagif.mass = null_mass;
        
        // Back up verbosity string.
        std::string verb_string = verb.verb_prefix;
        
        for (int n = 0; n < mgstruct.n_refine; n++)
        {
            if (n != 0)
            {
                level_down(&mgstruct);
            }
            
            cout << "[L" << mgstruct.curr_level+1 << "_NULLVEC]: Current fine size is " << mgstruct.curr_fine_size << "\n";
            
            // Update verb_prefix temporarily.
            verb.verb_prefix = "[L" + to_string(mgstruct.curr_level+1) + "_NULLVEC]: ";
            
            // We generate null vectors by solving Ax = 0, with a
            // gaussian initial guess.
            // For sanity with the residual, we really solve Ax = -Ax_0,
            // where x has a zero initial guess, x_0 is a random vector.
            complex<double>* rand_guess = new complex<double>[mgstruct.curr_fine_size];
            complex<double>* Arand_guess = new complex<double>[mgstruct.curr_fine_size];

            for (i = 0; i < mgstruct.n_vector/null_partitions; i++)
            {
                // Create a gaussian random source. 
                gaussian<double>(rand_guess, mgstruct.curr_fine_size, generator);

                // Make orthogonal to previous solutions.
                if (i > 0) // If there are vectors to orthogonalize against...
                {
                    for (j = 0; j < i; j++) // Iterate over all of them...
                    {
                        for (k = 0; k < (do_ortho_eo ? null_partitions : 1); k++) // And then iterate over even/odd or corners!
                        {
                            orthogonal<double>(rand_guess, mgstruct.null_vectors[mgstruct.curr_level][j*null_partitions+k], mgstruct.curr_fine_size);
                            if (do_global_ortho_conj)
                            {
                                conj<double>(mgstruct.null_vectors[mgstruct.curr_level][j*null_partitions+k], mgstruct.curr_fine_size);
                                orthogonal<double>(rand_guess, mgstruct.null_vectors[mgstruct.curr_level][j*null_partitions+k], mgstruct.curr_fine_size);
                                conj<double>(mgstruct.null_vectors[mgstruct.curr_level][j*null_partitions+k], mgstruct.curr_fine_size);
                            }
                        }
                    }
                }

                // Solve the residual equation. 
                zero<double>(Arand_guess, mgstruct.curr_fine_size);

                fine_square_staggered(Arand_guess, rand_guess, (void*)&mgstruct);

                //square_staggered_u1(Arand_guess, rand_guess, (void*)&stagif);
                for (j = 0; j < mgstruct.curr_fine_size; j++)
                {
                   Arand_guess[j] = -Arand_guess[j]; 
                }
                zero<double>(mgstruct.null_vectors[mgstruct.curr_level][null_partitions*i], mgstruct.curr_fine_size);

                
                switch (null_gen)
                {
                    case NULL_GCR:
                        minv_vector_gcr(mgstruct.null_vectors[mgstruct.curr_level][null_partitions*i], Arand_guess, mgstruct.curr_fine_size, null_max_iter, null_precision, fine_square_staggered, &mgstruct, &verb);
                        break;
                    case NULL_BICGSTAB:
                        minv_vector_bicgstab(mgstruct.null_vectors[mgstruct.curr_level][null_partitions*i], Arand_guess, mgstruct.curr_fine_size, null_max_iter, null_precision, fine_square_staggered, &mgstruct, &verb);
                        break;
                    case NULL_CG:
                        minv_vector_cg(mgstruct.null_vectors[mgstruct.curr_level][null_partitions*i], Arand_guess, mgstruct.curr_fine_size, null_max_iter, null_precision, fine_square_staggered, &mgstruct, &verb);
                        break;
                    case NULL_MINRES:
                        minv_vector_minres(mgstruct.null_vectors[mgstruct.curr_level][null_partitions*i], Arand_guess, mgstruct.curr_fine_size, null_max_iter, null_precision, fine_square_staggered, &mgstruct, &verb);
                        break;
                }


                for (j = 0; j < mgstruct.curr_fine_size; j++)
                {
                    mgstruct.null_vectors[mgstruct.curr_level][null_partitions*i][j] += rand_guess[j];
                }

                // Split into eo now if need be, otherwise we do it later.
                if (do_ortho_eo)
                {
                    // Aggregate in chirality (or corners) as needed.  
                    // This is handled differently if we're on the top level or further down. 
                    if (mgstruct.curr_level == 0) // top level routines
                    {
                        null_partition_staggered(&mgstruct, i, bstrat, &Lat);
                    }
                    else // not on the top level
                    {
                        null_partition_coarse(&mgstruct, i, bstrat);
                    }
                }
                

                // Normalize new vectors.
                for (k = 0; k < (do_ortho_eo ? null_partitions : 1); k++)
                {
                    normalize(mgstruct.null_vectors[mgstruct.curr_level][null_partitions*i+k], mgstruct.curr_fine_size);
                }

                // Orthogonalize against previous vectors. 
                if (i > 0)
                {
                    for (j = 0; j < i; j++)
                    {
                        cout << "[L" << mgstruct.curr_level+1 << "_NULLVEC]: Pre-orthog cosines of " << j << "," << i << " are: "; 
                        for (k = 0; k < (do_ortho_eo ? null_partitions : 1); k++)
                        {
                            // Check dot product before normalization.
                            cout << abs(dot<double>(mgstruct.null_vectors[mgstruct.curr_level][null_partitions*i+k], mgstruct.null_vectors[mgstruct.curr_level][null_partitions*j+k], mgstruct.curr_fine_size)/sqrt(norm2sq<double>(mgstruct.null_vectors[mgstruct.curr_level][null_partitions*i+k],mgstruct.curr_fine_size)*norm2sq<double>(mgstruct.null_vectors[mgstruct.curr_level][null_partitions*j+k],mgstruct.curr_fine_size))) << " ";

                            orthogonal<double>(mgstruct.null_vectors[mgstruct.curr_level][null_partitions*i+k], mgstruct.null_vectors[0][null_partitions*j+k], mgstruct.curr_fine_size); 
                            if (do_global_ortho_conj)
                            {
                                conj<double>(mgstruct.null_vectors[mgstruct.curr_level][null_partitions*j+k], mgstruct.curr_fine_size);
                                orthogonal<double>(mgstruct.null_vectors[mgstruct.curr_level][null_partitions*i+k], mgstruct.null_vectors[0][null_partitions*j+k], mgstruct.curr_fine_size); 
                                conj<double>(mgstruct.null_vectors[mgstruct.curr_level][null_partitions*j+k], mgstruct.curr_fine_size);
                            }
                        }
                        cout << "\n";
                    }
                }

                // Normalize again.
                for (k = 0; k < (do_ortho_eo ? null_partitions : 1); k++)
                {
                    normalize(mgstruct.null_vectors[mgstruct.curr_level][null_partitions*i+k], mgstruct.curr_fine_size);
                }
                
            }
        
            // If we didn't split null vectors before, we do it now.
            if (!do_ortho_eo)
            {
                for (i = 0; i < mgstruct.n_vector/null_partitions; i++)
                {
                    // Aggregate in chirality (or corners) as needed.  
                    // This is handled differently if we're on the top level or further down. 
                    if (mgstruct.curr_level == 0) // top level routines
                    {
                        null_partition_staggered(&mgstruct, i, bstrat, &Lat);
                    }
                    else // not on the top level
                    {
                        null_partition_coarse(&mgstruct, i, bstrat);
                    }
                    
                    for (k = 0; k < null_partitions; k++)
                    {
                        normalize(mgstruct.null_vectors[mgstruct.curr_level][null_partitions*i+k], mgstruct.curr_fine_size);
                    }
                }
            }

            delete[] rand_guess; 
            delete[] Arand_guess; 



    #ifdef PRINT_NULL_VECTOR
            // Print vector.
            /*
            cout << "\n\nPrinting null vectors:\n"; 
            for (int n = 0; n < mgstruct.n_vector; n++)
            {
                cout << "\nVector " << n << "\n";
                for (int y = 0; y < mgstruct.curr_y_fine; y++)
                {
                    for (int x = 0; x < mgstruct.curr_x_fine; x++)
                    {
                        cout << "(";
                        for (int c = 0; c < mgstruct.n_vector; c++)
                        {
                            cout << mgstruct.null_vectors[mgstruct.curr_level][n][y*mgstruct.curr_x_fine*mgstruct.curr_dof_fine+x*mgstruct.curr_dof_fine+c] << ",";
                        }
                        cout << ") ";
                    }
                    cout << "\n";
                }
            }*/
    #endif // PRINT_NULL_VECTOR

            cout << "[MG]: Performing block orthonormalize of null vectors...\n";
            block_orthonormalize(&mgstruct); 
            
        }
        
        // Reset the tricks we've done. 
        
        mgstruct.matrix_vector = op; // Reset op to solved op.
        stagif.mass = MASS; 
        verb.verb_prefix = verb_string; 
        
        // Un-pop to the finest level.
        for (int n = 1; n < mgstruct.n_refine; n++)
        {
            level_up(&mgstruct);
        }
        
        
    } // end skipping generation if we're only doing a top level or smoother test. 
    
#endif // generate null vector. 
    
    
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
    
    if (do_eigentest)
    {
        for (int lev = 0; lev < mgstruct.n_refine; lev++)
        {
            int n_eigen = 0;
            int n_cv = 0;
            complex<double>* evals = 0;
            complex<double>** evecs = 0;
            
            if (set_eigen == -1 && set_cv == -1) // generate all eigenvalues, eigenvectors. 
            {
                // Allocate space for all eigenvalues, eigenvectors. 
                n_eigen = mgstruct.curr_fine_size;
                n_cv = mgstruct.curr_fine_size; 
                cout << "[L" << lev+1 << "_ARPACK]: Number of eigenvalues: " << n_eigen << " Number of cv: " << n_cv << "\n";
                
                evals = new complex<double>[mgstruct.curr_fine_size];
                evecs = new complex<double>*[mgstruct.curr_fine_size];
                for (i = 0; i < n_eigen; i++)
                {
                    evecs[i] = new complex<double>[mgstruct.curr_fine_size];
                }
                
                // Get low mag half
                arpack_dcn_t* ar_strc = arpack_dcn_init(mgstruct.curr_fine_size, n_eigen/2, n_cv); // max eigenvectors, internal vecs
                char eigtype[3]; strcpy(eigtype, "SM"); // Smallest magnitude eigenvalues.
                arpack_solve_t info_solve = arpack_dcn_getev(ar_strc, evals, evecs, mgstruct.curr_fine_size, n_eigen/2, n_cv, 4000, eigtype, 1e-7, 0.0, fine_square_staggered, (void*)&mgstruct); 
                //arpack_dcn_free(&ar_strc);

                // Print info about the eigensolve.
                cout << "[L" << lev+1 << "_ARPACK]: Number of converged eigenvalues: " << info_solve.nconv << "\n";
                cout << "[L" << lev+1 << "_ARPACK]: Number of iteration steps: " << info_solve.niter << "\n";
                cout << "[L" << lev+1 << "_ARPACK]: Number of matrix multiplies: " << info_solve.nops << "\n";

                // Get high mag half
                //arpack_dcn_t* ar_strc = arpack_dcn_init(mgstruct.curr_fine_size, n_eigen, n_cv); // max eigenvectors, internal vecs
                strcpy(eigtype, "LM"); // Smallest magnitude eigenvalues.
                info_solve = arpack_dcn_getev(ar_strc, evals+(mgstruct.curr_fine_size/2), evecs+(mgstruct.curr_fine_size/2), mgstruct.curr_fine_size, n_eigen/2, n_cv, 4000, eigtype, 1e-7, 0.0, fine_square_staggered, (void*)&mgstruct); 
                arpack_dcn_free(&ar_strc);

                // Print info about the eigensolve.
                cout << "[L" << lev+1 << "_ARPACK]: Number of converged eigenvalues: " << info_solve.nconv << "\n";
                cout << "[L" << lev+1 << "_ARPACK]: Number of iteration steps: " << info_solve.niter << "\n";
                cout << "[L" << lev+1 << "_ARPACK]: Number of matrix multiplies: " << info_solve.nops << "\n";

                // End of arpack bindings!   
            }
            else if (set_eigen != -1 && set_cv == -1) // generate n_eigen eigenvalues, min(mgstruct.curr_fine_size, 2.5 n_eigen) cv.
            {
                n_eigen = set_eigen;
                n_cv = min(mgstruct.curr_fine_size, 2*n_eigen + n_eigen/2);
                cout << "[L" << lev+1 << "_ARPACK]: Number of eigenvalues: " << n_eigen << " Number of cv: " << n_cv << "\n";
                
                evals = new complex<double>[n_eigen];
                evecs = new complex<double>*[n_eigen];
                for (i = 0; i < n_eigen; i++)
                {
                    evecs[i] = new complex<double>[mgstruct.curr_fine_size];
                }
                
                // Get low mag half
                arpack_dcn_t* ar_strc = arpack_dcn_init(mgstruct.curr_fine_size, n_eigen, n_cv); // max eigenvectors, internal vecs
                char eigtype[3]; strcpy(eigtype, "SM"); // Smallest magnitude eigenvalues.
                arpack_solve_t info_solve = arpack_dcn_getev(ar_strc, evals, evecs, mgstruct.curr_fine_size, n_eigen, n_cv, 4000, eigtype, 1e-7, 0.0, fine_square_staggered, (void*)&mgstruct); 
                arpack_dcn_free(&ar_strc);
                
                // Print info about the eigensolve.
                cout << "[L" << lev+1 << "_ARPACK]: Number of converged eigenvalues: " << info_solve.nconv << "\n";
                cout << "[L" << lev+1 << "_ARPACK]: Number of iteration steps: " << info_solve.niter << "\n";
                cout << "[L" << lev+1 << "_ARPACK]: Number of matrix multiplies: " << info_solve.nops << "\n";
                
            }
            else // generate n_eigen eigenvalues, min(mgstruct.curr_fine_size, n_cv) cv.
            {
                n_eigen = set_eigen;
                n_cv = set_cv;
                
                cout << "[L" << lev+1 << "_ARPACK]: Number of eigenvalues: " << n_eigen << " Number of cv: " << n_cv << "\n";
                
                evals = new complex<double>[n_eigen];
                evecs = new complex<double>*[n_eigen];
                for (i = 0; i < n_eigen; i++)
                {
                    evecs[i] = new complex<double>[mgstruct.curr_fine_size];
                }
                
                // Get low mag half
                arpack_dcn_t* ar_strc = arpack_dcn_init(mgstruct.curr_fine_size, n_eigen, n_cv); // max eigenvectors, internal vecs
                char eigtype[3]; strcpy(eigtype, "SM"); // Smallest magnitude eigenvalues.
                arpack_solve_t info_solve = arpack_dcn_getev(ar_strc, evals, evecs, mgstruct.curr_fine_size, n_eigen, n_cv, 4000, eigtype, 1e-7, 0.0, fine_square_staggered, (void*)&mgstruct); 
                arpack_dcn_free(&ar_strc);
                
                // Print info about the eigensolve.
                cout << "[L" << lev+1 << "_ARPACK]: Number of converged eigenvalues: " << info_solve.nconv << "\n";
                cout << "[L" << lev+1 << "_ARPACK]: Number of iteration steps: " << info_solve.niter << "\n";
                cout << "[L" << lev+1 << "_ARPACK]: Number of matrix multiplies: " << info_solve.nops << "\n";
            }

            // Sort eigenvalues (done differently depending on the operator).
            for (i = 0; i < n_eigen; i++)
            {
                for (j = 0; j < n_eigen-1; j++)
                {
                    switch (opt)
                    {
                        case STAGGERED:
                            if (abs(imag(evals[j])) > abs(imag(evals[j+1])))
                            {
                                complex<double> teval = evals[j]; evals[j] = evals[j+1]; evals[j+1] = teval;
                                complex<double>* tevec = evecs[j]; evecs[j] = evecs[j+1]; evecs[j+1] = tevec;
                            }
                            break;
                        case LAPLACE:
                        case LAPLACE_NC2:
                        case G5_STAGGERED:
                        case STAGGERED_NORMAL:
                        case STAGGERED_INDEX:
                            if (abs(real(evals[j])) > abs(real(evals[j+1])))
                            {
                                complex<double> teval = evals[j]; evals[j] = evals[j+1]; evals[j+1] = teval;
                                complex<double>* tevec = evecs[j]; evecs[j] = evecs[j+1]; evecs[j+1] = tevec;
                            }
                            break; 
                    }
                }
            }
            
            
            cout << "\n\nAll eigenvalues:\n";
            for (i = 0; i < n_eigen; i++)
            {
                cout << "[L" << lev+1 << "_FINEVAL]: Mass " << MASS << " Num " << i << " Eval " << evals[i] << "\n";
                normalize<double>(evecs[i], mgstruct.curr_fine_size);
            }

            complex<double>* evec_Pdag = new complex<double>[mgstruct.curr_coarse_size];
            complex<double>* evec_Pdag2 = new complex<double>[mgstruct.curr_coarse_size];
            complex<double>* evec_PPdag = new complex<double>[mgstruct.curr_fine_size];

            // Test overlap of null vectors with eigenvectors.
            // Formally, this is looking at the magnitude of (1 - P P^\dag) eigenvector.

            for (i = 0; i < n_eigen; i++)
            {
                // Zero out.
                zero<double>(evec_Pdag, mgstruct.curr_coarse_size);
                zero<double>(evec_PPdag, mgstruct.curr_fine_size);

                // Restrict eigenvector.
                restrict(evec_Pdag, evecs[i], &mgstruct);

                // Prolong.
                prolong(evec_PPdag, evec_Pdag, &mgstruct);

                // Subtract off eigenvector, take norm.
                for (j = 0; j < mgstruct.curr_fine_size; j++)
                {
                    evec_PPdag[j] -= evecs[i][j];
                }

                cout << "[L" << lev+1 << "_1mPPDAG]: Num " << i << " Overlap " << sqrt(norm2sq<double>(evec_PPdag, mgstruct.curr_fine_size)) << "\n";
            }

            // Test how good of a preconditioner the coarse operator is.
            // Formally, this is looking at the magnitude of (1 - P ( P^\dag A P )^(-1) P^\dag A) eigenvector.
            for (i = 0; i < n_eigen; i++)
            {
                // Zero out.
                zero<double>(evec_Pdag, mgstruct.curr_coarse_size);
                zero<double>(evec_Pdag2, mgstruct.curr_coarse_size);
                zero<double>(evec_PPdag, mgstruct.curr_fine_size);

                // Apply A.
                fine_square_staggered(evec_PPdag, evecs[i], (void*)&mgstruct);

                // Restrict.
                restrict(evec_Pdag, evec_PPdag, &mgstruct);

                // Invert A_coarse against it.
                invif = minv_vector_gcr_restart(evec_Pdag2, evec_Pdag, mgstruct.curr_coarse_size, 10000, 1e-7, 64, coarse_square_staggered, (void*)&mgstruct);

                // Prolong.
                zero<double>(evec_PPdag, mgstruct.curr_coarse_size);
                prolong(evec_PPdag, evec_Pdag2, &mgstruct);

                // Subtract off eigenvector, take norm.
                for (j = 0; j < mgstruct.curr_fine_size; j++)
                {
                    evec_PPdag[j] -= evecs[i][j];
                }

                cout << "[L" << lev+1 << "_1mP_Ac_PDAG_A]: Num " << i << " Overlap " << sqrt(norm2sq<double>(evec_PPdag, mgstruct.curr_fine_size)) << "\n";
            }


            delete[] evec_Pdag;
            delete[] evec_PPdag;
            delete[] evec_Pdag2;

            for (i = 0; i < n_eigen; i++)
            {
                delete[] evecs[i];
            }
            delete[] evecs;
            delete[] evals;

            if (lev < mgstruct.n_refine-1)
            {
                level_down(&mgstruct);
            }
        }

        for (int lev = mgstruct.n_refine-2; lev >= 0; lev--)
        {
            level_up(&mgstruct);
        }
    
    } // do_eigentest
#endif // EIGEN_TEST
    
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
    
#endif // COARSE_ONLY
    
    if (my_test == TOP_LEVEL_ONLY)
    {
        // Try a direct solve.
        cout << "\n[ORIG]: Solve fine system.\n";

        invif = minv_vector_gcr_restart(lhs, rhs, Lat.get_lattice_size(), 100000, outer_precision, outer_restart, op, (void*)&stagif, &verb);
        //invif = minv_vector_gcr_restart(lhs, rhs, Lat.get_lattice_size(), 100000, outer_precision, outer_restart, square_staggered_u1, (void*)&stagif);

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
        (*op)(check, lhs, (void*)&stagif); 
        //square_staggered_u1(check, lhs, (void*)&stagif);

        explicit_resid = 0.0;
        for (i = 0; i < Lat.get_lattice_size(); i++)
        {
          explicit_resid += real(conj(rhs[i] - check[i])*(rhs[i] - check[i]));
        }
        explicit_resid = sqrt(explicit_resid)/bnorm;

        printf("[ORIG]: [check] should equal [rhs]. The relative residual is %15.20e.\n", explicit_resid);

    } // TOP_LEVEL_ONLY
    
#ifdef COARSE_ONLY
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

    if (my_test == SMOOTHER_ONLY || my_test == TWO_LEVEL || my_test == THREE_LEVEL)
    {
        // Let's actually test a multigrid solve!
        cout << "\n[MG]: Test MG solve.\n";

        // Block normalize the null vectors.
        block_normalize(&mgstruct); 

        

        // Well, maybe this will work?
        zero<double>(lhs, Lat.get_lattice_size());
        //invif = minv_vector_cg_flex_precond_restart(lhs, rhs, Lat.get_lattice_size(), 10000, outer_precision, outer_restart, op, (void*)&stagif, mg_preconditioner, (void*)&mgprecond, &verb); 
        invif = minv_vector_gcr_var_precond_restart(lhs, rhs, Lat.get_lattice_size(), 10000, outer_precision, outer_restart, op, (void*)&stagif, mg_preconditioner, (void*)&mgprecond, &verb); 
        //invif = minv_vector_gcr_var_precond_restart(lhs, rhs, Lat.get_lattice_size(), 10000, outer_precision, outer_restart, square_staggered_u1, (void*)&stagif, mg_preconditioner, (void*)&mgprecond); 

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
    delete[] mgstruct.null_vectors; 
    
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



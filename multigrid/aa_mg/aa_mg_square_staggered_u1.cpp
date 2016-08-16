// This is a sample code that constructs general kinetic terms!

#include <iostream>
#include <iomanip> // to set output precision.
#include <cmath>
#include <string>
#include <sstream>
#include <complex>
#include <random>
#include <cstring> // should be replaced by using sstream


#include "generic_inverters.h"
#include "generic_inverters_precond.h"
#include "generic_vector.h"
#include "verbosity.h"
#include "mg.h"
#include "mg_real.h"
#include "mg_complex.h"
#include "u1_utils.h"
#include "lattice.h"

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

// How should we generate null vectors?
enum mg_null_gen_type
{
    NULL_GCR = 0,                    // Generate null vectors with GCR
    NULL_BICGSTAB = 1,               // Generate null vectors with BiCGStab
    NULL_CG = 2,                    // Generate null vectors with CG
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

// Square laplace 2d operator w/out u1 function.
void square_laplace(complex<double>* lhs, complex<double>* rhs, void* extra_data);

// Square staggered 2d operator w/out u1 function.
void square_staggered(complex<double>* lhs, complex<double>* rhs, void* extra_data);

// Square staggered 2d operator w/ u1 function.
void square_staggered_u1(complex<double>* lhs, complex<double>* rhs, void* extra_data);

// \gamma_5
void gamma_5(complex<double>* lhs, complex<double>* rhs, void* extra_data);

// Square \gamma_5 staggered 2d operator w/ u1 function.
void square_staggered_gamma5_u1(complex<double>* lhs, complex<double>* rhs, void* extra_data);  

// Staggered normal equations.
void square_staggered_normal_u1(complex<double>* lhs, complex<double>* rhs, void* extra_data);  

enum op_type
{
    STAGGERED = 0,
    LAPLACE = 1,
    LAPLACE_NC2 = 2,
    G5_STAGGERED = 3,
    STAGGERED_NORMAL = 4
};

enum src_type
{
    POINT = 0,
    RANDOM_GAUSSIAN = 1,
    ORIGIN_POINT = 2 // for correlator test. 
};

struct staggered_u1_op
{
    complex<double> *lattice;
    double mass;
    int x_fine;
    int y_fine; 
    Lattice* Lat; 
    int Nc; // only relevant for square laplace. 
};

int main(int argc, char** argv)
{  
    // Declare some variables.
    cout << setiosflags(ios::scientific) << setprecision(6);
    int i, j, k, x, y;
    complex<double> *lattice; // Holds the gauge field.
    complex<double> *lhs, *rhs, *check; // For some Kinetic terms.
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
    mg_null_gen_type null_gen = NULL_BICGSTAB; // NULL_BICGSTAB, NULL_GCR, NULL_CG
    
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
    int eo = 1; // 0 for no even/odd aggregation, 1 for even/odd aggregation.
    
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
    
    
    // Inner solver.
    inner_solver in_solve = GCR; //CR; //GCR; 
    double inner_precision = 1e-3;
    int inner_restart = 64;
    int inner_max = 1000;
    if (my_test == SMOOTHER_ONLY)
    {
        in_solve = NONE; 
    }
    
    // Smoother
    inner_solver in_smooth = GCR; //NONE; //GCR; BICGSTAB
    double omega_smooth = 0.67; // for MR only. 
    int pre_smooth = 6; // Can set on command line.
    int post_smooth = 6; // Can set on command line.
    
    // Gauge field information.
    double BETA = 6.0; // For random gauge field, phase angles have std.dev. 1/sqrt(beta).
                       // For heatbath gauge field, corresponds to non-compact beta.
                       // Can be set on command line with --beta.
    
    
    /////////////////////////////////////////////
    // Get a few parameters from command line. //
    /////////////////////////////////////////////
    for (i = 1; i < argc; i++)
    {
        if (strcmp(argv[i], "--help") == 0)
        {
            cout << "--lattice-size [32, 64, 128] {##}      (default 32x32)\n";
            cout << "--operator [laplace, laplace2, staggered\n";
            cout << "       g5_staggered, normal_staggered] (default staggered)\n";
            cout << "--null-operator [laplace, laplace2, staggered\n";
            cout << "       g5_staggered, normal_staggered] (default staggered)\n";
            cout << "--null-solver [gcr, bicgstab, cg]      (default bicgstab)\n";
            cout << "--null-precision [null prec]           (default 5e-5)\n";
            cout << "--null-eo [corner, yes, no]            (default yes)\n";
            cout << "--mass [mass]                          (default 1e-2)\n";
            cout << "--blocksize [blocksize]                (default 4)\n";
            cout << "--nvec [nvec]                          (default 4)\n";
            cout << "--nrefine [number coarse]              (default 1)\n";
            cout << "--gauge [unit, load, random]           (default load)\n";
            cout << "--gauge-transform [yes, no]            (default no)\n";
            cout << "--beta [3.0, 6.0, 10.0, 10000.0]       (default 6.0)\n";
            cout << "--npre-smooth [presmooth steps]        (default 6)\n";
            cout << "--npost-smooth [postsmooth steps]      (default 6)\n";
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
                i++;
            }
            else if (strcmp(argv[i], "--null-eo") == 0)
            {
                if (strcmp(argv[i+1], "yes") == 0)
                {
                    eo = 1;
                }
                else if (strcmp(argv[i+1], "corner") == 0)
                {
                    eo = 3; //Weird but hacky, I know.
                }
                else // none.
                {
                    eo = 0;
                }
                i++;
            }
            else if (strcmp(argv[i], "--nrefine") == 0)
            {
                n_refine = atoi(argv[i+1]);
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
                    lattice_size_y = lattice_size_y; // At the end, don't try to grab the next element!
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
            else
            {
                cout << argv[i] << " is not a valid flag.\n";
                cout << "--lattice-size [32, 64, 128] {##}      (default 32x32)\n";
                cout << "--operator [laplace, laplace2, staggered\n";
                cout << "       g5_staggered, normal_staggered] (default staggered)\n";
                cout << "--null-operator [laplace, laplace2, staggered\n";
                cout << "       g5_staggered, normal_staggered] (default staggered)\n";
                cout << "--null-solver [gcr, bicgstab, cg]      (default bicgstab)\n";
                cout << "--null-precision [null prec]           (default 5e-5)\n";
                cout << "--null-eo [yes, no]                    (default yes)\n";
                cout << "--mass [mass]                          (default 1e-2)\n";
                cout << "--blocksize [blocksize]                (default 4)\n";
                cout << "--nvec [nvec]                          (default 4)\n";
                cout << "--nrefine [number coarse]              (default 1)\n";
                cout << "--gauge [unit, load, random]           (default load)\n";
                cout << "--gauge-transform [yes, no]            (default no)\n";
                cout << "--beta [3.0, 6.0, 10.0, 10000.0]       (default 6.0)\n";
                cout << "--npre-smooth [presmooth steps]        (default 6)\n";
                cout << "--npost-smooth [postsmooth steps]      (default 6)\n";
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
        eo = 0;
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
    
    // In case we need a gauge transform.
    complex<double>* gauge_trans = new complex<double>[Lat.get_lattice_size()];
    
    // Zero it out.
    zero<double>(lattice, 2*Lat.get_lattice_size());
    zero<double>(rhs, Lat.get_lattice_size());
    zero<double>(lhs, Lat.get_lattice_size());
    zero<double>(check, Lat.get_lattice_size());
    zero<double>(gauge_trans, Lat.get_lattice_size());
    //
    
    // Fill stagif.
    stagif.lattice = lattice;
    stagif.mass = MASS; 
    stagif.x_fine = Lat.get_lattice_dimension(0); // DEPRECIATE
    stagif.y_fine = Lat.get_lattice_dimension(1); // DEPRECIATE
    stagif.Lat = &Lat; 
    stagif.Nc = Nc; // Only relevant for laplace test only.
    
    // Create the verbosity structure.
    inversion_verbose_struct verb;
    verb.verbosity = VERB_DETAIL;
    verb.verb_prefix = "[L1]: ";
    verb.precond_verbosity = VERB_NONE;
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
            internal_load_gauge_u1(lattice, Lat.get_lattice_dimension(0), Lat.get_lattice_dimension(1), BETA); // defined at end of file.
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
    for (i = 0; i < n_refine; i++)
    {
        mgstruct.blocksize_x[i] = X_BLOCKSIZE;
        mgstruct.blocksize_y[i] = Y_BLOCKSIZE;
    }


#if defined GEN_NULL_VECTOR
    mgstruct.eo = eo;
    mgstruct.n_vector = (eo+1)*n_null_vector; //eo+1 = 1 for no e/o, 2 for even/odd, 4 for corner. 
#else
    mgstruct.eo = eo;
    mgstruct.n_vector = (eo+1)*n_null_vector;
#endif
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

    mgprecond.in_smooth_type = in_smooth; // What inner smoother? MR or GCR.
    mgprecond.omega_smooth = omega_smooth; // What relaxation parameter should we use (MR only!)
    mgprecond.n_pre_smooth = pre_smooth; // 6 MR smoother steps before coarsening.
    mgprecond.n_post_smooth = post_smooth; // 6 MR smoother steps after refining.
    mgprecond.in_solve_type = in_solve; // What inner solver? NONE, MR, CG, GCR, BICGSTAB
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
        
        printf("About to generate null vector.\n"); fflush(stdout);
        
        // We generate null vectors by solving Ax = 0, with a
        // gaussian initial guess.
        // For sanity with the residual, we really solve Ax = -Ax_0,
        // where x has a zero initial guess, x_0 is a random vector.
        complex<double>* rand_guess = new complex<double>[Lat.get_lattice_size()];
        complex<double>* Arand_guess = new complex<double>[Lat.get_lattice_size()];
    
        // Temporarily set the mass to zero for the null vector generation. 
        stagif.mass = stagif.mass*1e-2;

        for (i = 0; i < mgstruct.n_vector/(mgstruct.eo+1); i++)
        {
            // Create a gaussian random source. 
            gaussian<double>(rand_guess, Lat.get_lattice_size(), generator);
            
            // Make orthogonal to previous solutions.
            if (i > 0) // If there are vectors to orthogonalize against...
            {
                for (j = 0; j < i; j++) // Iterate over all of them...
                {
                    for (k = 0; k <= mgstruct.eo; k++) // And then iterate over even/odd or color!
                    {
                        orthogonal<double>(rand_guess, mgstruct.null_vectors[0][j*(mgstruct.eo+1)+k], Lat.get_lattice_size());
                    }
                    /*if (mgstruct.eo == 1)
                    {
                        orthogonal<double>(rand_guess, mgstruct.null_vectors[0][j], Lat.get_lattice_size()); 
                        orthogonal<double>(rand_guess, mgstruct.null_vectors[0][j+mgstruct.n_vector/2], Lat.get_lattice_size()); 
                    }
                    else if (mgstruct.eo == 3) // corner.
                    {
                        orthogonal<double>(rand_guess, mgstruct.null_vectors[0][j], Lat.get_lattice_size()); 
                        orthogonal<double>(rand_guess, mgstruct.null_vectors[0][j+mgstruct.n_vector/4], Lat.get_lattice_size()); 
                        orthogonal<double>(rand_guess, mgstruct.null_vectors[0][j+2*mgstruct.n_vector/4], Lat.get_lattice_size()); 
                        orthogonal<double>(rand_guess, mgstruct.null_vectors[0][j+3*mgstruct.n_vector/4], Lat.get_lattice_size()); 
                    }
                    else // no eo.
                    {
                        orthogonal<double>(rand_guess, mgstruct.null_vectors[0][j], Lat.get_lattice_size()); 
                    }*/
                }
            }

            // Solve the residual equation. 
            zero<double>(Arand_guess, Lat.get_lattice_size());
            
            (*op_null)(Arand_guess, rand_guess, (void*)&stagif);
            
            //square_staggered_u1(Arand_guess, rand_guess, (void*)&stagif);
            for (j = 0; j < Lat.get_lattice_size(); j++)
            {
               Arand_guess[j] = -Arand_guess[j]; 
            }
            zero<double>(mgstruct.null_vectors[0][(mgstruct.eo+1)*i], Lat.get_lattice_size());

            switch (null_gen)
            {
                case NULL_GCR:
                    minv_vector_gcr(mgstruct.null_vectors[0][(mgstruct.eo+1)*i], Arand_guess, Lat.get_lattice_size(), null_max_iter, null_precision, op_null, (void*)&stagif, &verb); 
                    break;
                case NULL_BICGSTAB:
                    minv_vector_bicgstab(mgstruct.null_vectors[0][(mgstruct.eo+1)*i], Arand_guess, Lat.get_lattice_size(), null_max_iter, null_precision, op_null, (void*)&stagif, &verb); 
                    break;
                case NULL_CG:
                    minv_vector_cg(mgstruct.null_vectors[0][(mgstruct.eo+1)*i], Arand_guess, Lat.get_lattice_size(), null_max_iter, null_precision, op_null, (void*)&stagif, &verb); 
                    break;
            }
            

            for (j = 0; j < Lat.get_lattice_size(); j++)
            {
                mgstruct.null_vectors[0][(mgstruct.eo+1)*i][j] += rand_guess[j];
            }
            
            // Aggregate in chirality (or corners) as needed, orthogonalize against previous vectors.  
            if (mgstruct.eo == 1)
            {
                for (j = 0; j < Lat.get_lattice_size(); j++)
                {
                    if (Lat.index_is_even(j))
                    {
                        mgstruct.null_vectors[0][2*i+1][j] = mgstruct.null_vectors[0][2*i][j];
                        mgstruct.null_vectors[0][2*i][j] = 0.0;
                    }
                }
                normalize(mgstruct.null_vectors[0][2*i], Lat.get_lattice_size());
                normalize(mgstruct.null_vectors[0][2*i+1], Lat.get_lattice_size());
                
                // Orthogonalize against previous vectors. 
                if (i > 0)
                {
                    for (j = 0; j < i; j++)
                    {
                        // Check dot product before normalization.
                    cout << "[NULLVEC]: Pre-orthog cosines of " << j << "," << i << " are: " <<
                        abs(dot<double>(mgstruct.null_vectors[0][2*i], mgstruct.null_vectors[0][2*j], Lat.get_lattice_size())/sqrt(norm2sq<double>(mgstruct.null_vectors[0][2*i],Lat.get_lattice_size())*norm2sq<double>(mgstruct.null_vectors[0][2*j],Lat.get_lattice_size()))) << " " << 
                        abs(dot<double>(mgstruct.null_vectors[0][2*i+1], mgstruct.null_vectors[0][2*j+1], Lat.get_lattice_size())/sqrt(norm2sq<double>(mgstruct.null_vectors[0][2*i+1],Lat.get_lattice_size())*norm2sq<double>(mgstruct.null_vectors[0][2*j+1],Lat.get_lattice_size()))) << "\n"; 
                        
                        orthogonal<double>(mgstruct.null_vectors[0][2*i], mgstruct.null_vectors[0][2*j], Lat.get_lattice_size()); 
                        orthogonal<double>(mgstruct.null_vectors[0][2*i+1], mgstruct.null_vectors[0][2*j+1], Lat.get_lattice_size()); 
                    }
                }
                
                normalize(mgstruct.null_vectors[0][2*i], Lat.get_lattice_size());
                normalize(mgstruct.null_vectors[0][2*i+1], Lat.get_lattice_size());

            }
            else if (mgstruct.eo == 3) // corner
            {
                for (j = 0; j < Lat.get_lattice_size(); j++)
                {
                    // Find x and y component. 
                    Lat.index_to_coord(j, coord, nd);
                    if (coord[0]%2 == 1 && coord[1]%2 == 0)
                    {
                        mgstruct.null_vectors[0][4*i+1][j] = mgstruct.null_vectors[0][4*i][j];
                        mgstruct.null_vectors[0][4*i][j] = 0.0;
                    }
                    else if (coord[0]%2 == 0 && coord[1]%2 == 1)
                    {
                        mgstruct.null_vectors[0][4*i+2][j] = mgstruct.null_vectors[0][4*i][j];
                        mgstruct.null_vectors[0][4*i][j] = 0.0;
                    }
                    else if (coord[0]%2 == 1 && coord[1]%2 == 1)
                    {
                        mgstruct.null_vectors[0][4*i+3][j] = mgstruct.null_vectors[0][4*i][j];
                        mgstruct.null_vectors[0][4*i][j] = 0.0;
                    }
                }
                normalize(mgstruct.null_vectors[0][4*i], Lat.get_lattice_size());
                normalize(mgstruct.null_vectors[0][4*i+1], Lat.get_lattice_size());
                normalize(mgstruct.null_vectors[0][4*i+2], Lat.get_lattice_size());
                normalize(mgstruct.null_vectors[0][4*i+3], Lat.get_lattice_size());

                if (i > 0)
                {
                    for (j = 0; j < i; j++)
                    {
                        // Check dot product before normalization.
                        cout << "[NULLVEC]: Pre-orthog cosines of " << j << "," << i << " are: " <<
                            abs(dot<double>(mgstruct.null_vectors[0][4*i], mgstruct.null_vectors[0][4*j], Lat.get_lattice_size())/sqrt(norm2sq<double>(mgstruct.null_vectors[0][4*i],Lat.get_lattice_size())*norm2sq<double>(mgstruct.null_vectors[0][4*j],Lat.get_lattice_size()))) << " " << 
                            abs(dot<double>(mgstruct.null_vectors[0][4*i+1], mgstruct.null_vectors[0][4*j+1], Lat.get_lattice_size())/sqrt(norm2sq<double>(mgstruct.null_vectors[0][4*i+1],Lat.get_lattice_size())*norm2sq<double>(mgstruct.null_vectors[0][4*j+1],Lat.get_lattice_size()))) << " " << 
                            abs(dot<double>(mgstruct.null_vectors[0][4*i+2], mgstruct.null_vectors[0][4*j+2], Lat.get_lattice_size())/sqrt(norm2sq<double>(mgstruct.null_vectors[0][4*i+2],Lat.get_lattice_size())*norm2sq<double>(mgstruct.null_vectors[0][4*j+2],Lat.get_lattice_size()))) << " " << 
                            abs(dot<double>(mgstruct.null_vectors[0][4*i+3], mgstruct.null_vectors[0][4*j+3], Lat.get_lattice_size())/sqrt(norm2sq<double>(mgstruct.null_vectors[0][4*i+3],Lat.get_lattice_size())*norm2sq<double>(mgstruct.null_vectors[0][4*j+3],Lat.get_lattice_size()))) << " " << "\n"; 

                        orthogonal<double>(mgstruct.null_vectors[0][4*i], mgstruct.null_vectors[0][4*j], Lat.get_lattice_size()); 
                        orthogonal<double>(mgstruct.null_vectors[0][4*i+1], mgstruct.null_vectors[0][4*j+1], Lat.get_lattice_size()); 
                        orthogonal<double>(mgstruct.null_vectors[0][4*i+2], mgstruct.null_vectors[0][4*j+2], Lat.get_lattice_size()); 
                        orthogonal<double>(mgstruct.null_vectors[0][4*i+3], mgstruct.null_vectors[0][4*j+3], Lat.get_lattice_size()); 
                    }
                }

                normalize(mgstruct.null_vectors[0][4*i], Lat.get_lattice_size());
                normalize(mgstruct.null_vectors[0][4*i+1], Lat.get_lattice_size());
                normalize(mgstruct.null_vectors[0][4*i+2], Lat.get_lattice_size());
                normalize(mgstruct.null_vectors[0][4*i+3], Lat.get_lattice_size());

            }
            else // none
            {
                normalize(mgstruct.null_vectors[0][i], Lat.get_lattice_size());
                if (i > 0)
                {
                    for (j = 0; j < i; j++)
                    {
                        // Check dot product before normalization.
                    cout << "[NULLVEC]: Pre-orthog cosines of " << j << "," << i << " are: " <<
                        abs(dot<double>(mgstruct.null_vectors[0][i], mgstruct.null_vectors[0][j], Lat.get_lattice_size())/sqrt(norm2sq<double>(mgstruct.null_vectors[0][i],Lat.get_lattice_size())*norm2sq<double>(mgstruct.null_vectors[0][j],Lat.get_lattice_size()))) << "\n"; 
                        
                        orthogonal<double>(mgstruct.null_vectors[0][i], mgstruct.null_vectors[0][j], Lat.get_lattice_size()); 
                    }
                }
                normalize(mgstruct.null_vectors[0][i], Lat.get_lattice_size());
            }

        }

        // This causes a segfault related to the RNG when
        // the vector is initialized.
        delete[] rand_guess; 
        delete[] Arand_guess; 


    
#ifdef PRINT_NULL_VECTOR
        cout << "Check projector:\n"; 
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
        
        cout << "[MG]: Performing block orthonormalize of null vectors...\n";
        block_orthonormalize(&mgstruct); 
        
        
        // Do we need to generate more levels?
        if (mgstruct.n_refine > 1)
        {
            mgstruct.matrix_vector = op_null; // Trick op to null gen op.
            for (int n = 1; n < mgstruct.n_refine; n++)
            {
                level_down(&mgstruct);
                cout << "curr_fine_size: " << mgstruct.curr_fine_size << "\n";
                
                // Let's give it a whirl?!
                // We generate null vectors by solving Ax = 0, with a
                // gaussian initial guess.
                // For sanity with the residual, we really solve Ax = -Ax_0,
                // where x has a zero initial guess, x_0 is a random vector.
                complex<double>* c_rand_guess = new complex<double>[mgstruct.curr_fine_size];
                complex<double>* c_Arand_guess = new complex<double>[mgstruct.curr_fine_size];

                // Generate null vectors with the current level. 
                for (i = 0; i < mgstruct.n_vector/(mgstruct.eo+1); i++)
                {
                    gaussian<double>(c_rand_guess, mgstruct.curr_fine_size, generator);
                    
                    // Make orthogonal to previous solutions.
                    if (i > 0)
                    {
                        for (j = 0; j < i; j++)
                        {
                            if (mgstruct.eo == 1)
                            {
                                orthogonal<double>(c_rand_guess, mgstruct.null_vectors[mgstruct.curr_level][j], mgstruct.curr_fine_size); 
                                orthogonal<double>(c_rand_guess, mgstruct.null_vectors[mgstruct.curr_level][j+mgstruct.n_vector/2], mgstruct.curr_fine_size);
                            }
                            else // no eo.
                            {
                                orthogonal<double>(c_rand_guess, mgstruct.null_vectors[mgstruct.curr_level][j], mgstruct.curr_fine_size); 
                            }
                        }
                    }

                    zero<double>(c_Arand_guess, mgstruct.curr_fine_size); 
                    fine_square_staggered(c_Arand_guess, c_rand_guess, (void*)&mgstruct);
                    
                    
                    for (j = 0; j < mgstruct.curr_fine_size; j++)
                    {
                       c_Arand_guess[j] = -c_Arand_guess[j]; 
                    }
                    
                    zero<double>(mgstruct.null_vectors[mgstruct.curr_level][i], mgstruct.curr_fine_size);
                    
                    switch (null_gen)
                    {
                        case NULL_GCR:
                            minv_vector_gcr(mgstruct.null_vectors[mgstruct.curr_level][i], c_Arand_guess, mgstruct.curr_fine_size, null_max_iter, null_precision, fine_square_staggered, &mgstruct, &verb);
                            break;
                        case NULL_BICGSTAB:
                            minv_vector_bicgstab(mgstruct.null_vectors[mgstruct.curr_level][i], c_Arand_guess, mgstruct.curr_fine_size, null_max_iter, null_precision, fine_square_staggered, &mgstruct, &verb);
                            break;
                        case NULL_CG:
                            minv_vector_cg(mgstruct.null_vectors[mgstruct.curr_level][i], c_Arand_guess, mgstruct.curr_fine_size, null_max_iter, null_precision, fine_square_staggered, &mgstruct, &verb);
                            break;
                    }
                    

                    for (j = 0; j < mgstruct.curr_fine_size; j++)
                    {
                        mgstruct.null_vectors[mgstruct.curr_level][i][j] += c_rand_guess[j];
                    }


                    // Aggregate in chirality, orthogonalize against previous vectors.
                    if (mgstruct.eo == 1)
                    {
                        
                        for (j = 0; j < mgstruct.curr_fine_size; j++)
                        {
                            int c = j % mgstruct.n_vector; // What color index do we have?
                                                           // 0 to mgstruct.n_vector/2-1 is even, else is odd.
                            //int x_coord = (i - c)/mgstruct.n_vector % mgstruct.curr_x_fine;
                            //int y_coord = ((i - c)/mgstruct.n_vector - x_coord)/mgstruct.curr_x_fine;
                            
                            if (c >= mgstruct.n_vector/2)
                            {
                                mgstruct.null_vectors[mgstruct.curr_level][i+mgstruct.n_vector/2][j] = mgstruct.null_vectors[mgstruct.curr_level][i][j];
                                mgstruct.null_vectors[mgstruct.curr_level][i][j] = 0.0;
                            }
                        }
                        normalize(mgstruct.null_vectors[mgstruct.curr_level][i], mgstruct.curr_fine_size);
                        normalize(mgstruct.null_vectors[mgstruct.curr_level][i+mgstruct.n_vector/2], mgstruct.curr_fine_size);
                        
                        if (i > 0)
                        {
                            for (j = 0; j < i; j++)
                            {
                                // Check dot product before normalization.
                                cout << "[NULLVEC]: Pre-orthog cosines of " << j << "," << i << " are: " <<
                                    abs(dot<double>(mgstruct.null_vectors[mgstruct.curr_level][i], mgstruct.null_vectors[mgstruct.curr_level][j], mgstruct.curr_fine_size)/sqrt(norm2sq<double>(mgstruct.null_vectors[mgstruct.curr_level][i],mgstruct.curr_fine_size)*norm2sq<double>(mgstruct.null_vectors[mgstruct.curr_level][j],mgstruct.curr_fine_size))) << " " << 
                                    abs(dot<double>(mgstruct.null_vectors[mgstruct.curr_level][i+mgstruct.n_vector/2], mgstruct.null_vectors[mgstruct.curr_level][j+mgstruct.n_vector/2], mgstruct.curr_fine_size)/sqrt(norm2sq<double>(mgstruct.null_vectors[mgstruct.curr_level][i+mgstruct.n_vector/2],mgstruct.curr_fine_size)*norm2sq<double>(mgstruct.null_vectors[mgstruct.curr_level][j+mgstruct.n_vector/2],mgstruct.curr_fine_size))) << "\n"; 
                                
                                    orthogonal<double>(mgstruct.null_vectors[mgstruct.curr_level][i], mgstruct.null_vectors[mgstruct.curr_level][j], mgstruct.curr_fine_size); 
                                    orthogonal<double>(mgstruct.null_vectors[mgstruct.curr_level][i+mgstruct.n_vector/2], mgstruct.null_vectors[mgstruct.curr_level][j+mgstruct.n_vector/2], mgstruct.curr_fine_size); 
                            }
                            
                            
                        }
                        
                        normalize(mgstruct.null_vectors[mgstruct.curr_level][i], mgstruct.curr_fine_size);
                        normalize(mgstruct.null_vectors[mgstruct.curr_level][i+mgstruct.n_vector/2], mgstruct.curr_fine_size);

                    } // Need to preserve 4 corners...
                    else
                    {
                        normalize(mgstruct.null_vectors[mgstruct.curr_level][i], mgstruct.curr_fine_size);
                        
                        if (i > 0)
                        {
                            for (j = 0; j < i; j++)
                            {
                                // Check dot product before normalization.
                                cout << "[NULLVEC]: Pre-orthog cosines of " << j << "," << i << " are: " <<
                                    abs(dot<double>(mgstruct.null_vectors[mgstruct.curr_level][i], mgstruct.null_vectors[mgstruct.curr_level][j], mgstruct.curr_fine_size)/sqrt(norm2sq<double>(mgstruct.null_vectors[mgstruct.curr_level][i],mgstruct.curr_fine_size)*norm2sq<double>(mgstruct.null_vectors[mgstruct.curr_level][j],mgstruct.curr_fine_size))) << "\n"; 
                                
                                    orthogonal<double>(mgstruct.null_vectors[mgstruct.curr_level][i], mgstruct.null_vectors[mgstruct.curr_level][j], mgstruct.curr_fine_size); 
                            }
                        }
                        
                        normalize(mgstruct.null_vectors[mgstruct.curr_level][i], mgstruct.curr_fine_size);
                    }

                }

                delete[] c_rand_guess; 
                delete[] c_Arand_guess; 
         
                block_orthonormalize(&mgstruct); 
                
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
            }
            
            // Un-pop to the finest level.
            for (i = 1; i < mgstruct.n_refine; i++)
            {
                level_up(&mgstruct);
            }
            
            mgstruct.matrix_vector = op; // Reset op to solved op.
        }
        
        // Reset the mass.
        stagif.mass = MASS; 
        
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

// Square lattice.
// Kinetic term for a 2D laplace w/ period bc. Applies lhs = A*rhs.
// The unit vectors are e_1 = xhat, e_2 = yhat.
void square_laplace(complex<double>* lhs, complex<double>* rhs, void* extra_data)
{
    // Declare variables.
    int i;
    int x,y,c;
    int tmp; 
    staggered_u1_op* stagif = (staggered_u1_op*)extra_data;
    double mass = stagif->mass; 
    int x_fine = stagif->x_fine;
    int y_fine = stagif->y_fine;
    int Nc = stagif->Nc; // 1 is the trivial laplace. 
    
    // Apply the stencil.
    for (i = 0; i < x_fine*y_fine*Nc; i++)
    {
        lhs[i] = 0.0;
        c = i%Nc;      // get color.
        tmp = (i-c)/Nc;
        x = tmp%x_fine; // integer mod.
        y = tmp/x_fine; // integer divide.

        // + e1.
        lhs[i] = lhs[i]-rhs[y*x_fine*Nc+((x+1)%x_fine)*Nc+c];
        
        // - e1.
        lhs[i] = lhs[i]- rhs[y*x_fine*Nc+((x+x_fine-1)%x_fine)*Nc+c]; // The extra +N is because of the % sign convention.

        // + e2.
        lhs[i] = lhs[i]- rhs[((y+1)%y_fine)*x_fine*Nc+x*Nc+c];

        // - e2.
        lhs[i] = lhs[i]- rhs[((y+y_fine-1)%y_fine)*x_fine*Nc+x*Nc+c];

        // 0
        // Added mass term here.
        lhs[i] = lhs[i]+ (4+mass)*rhs[i];
    }

}

// Square lattice.
// Kinetic term for a 2D staggered w/ period bc. Applies lhs = A*rhs.
// The unit vectors are e_1 = xhat, e_2 = yhat.
// The "extra_data" doesn't include anything.
// Apsi[x][y] = m psi[x,y] - U[y][x,x+1] 
void square_staggered(complex<double>* lhs, complex<double>* rhs, void* extra_data)
{
    // Declare variables.
    int i;
    int x,y;
    double eta1;
    staggered_u1_op* stagif = (staggered_u1_op*)extra_data;
    double mass = stagif->mass; 
    int x_fine = stagif->x_fine;
    int y_fine = stagif->y_fine;

    // For a 2D square lattice, the stencil is:
    //   1 |  0 -eta1  0 |
    //   - | +1    0  -1 |  , where eta1 = (-1)^x
    //   2 |  0 +eta1  0 |
    //
    // e2 = yhat
    // ^
    // | 
    // |-> e1 = xhat

    // Apply the stencil.
    for (i = 0; i < x_fine*y_fine; i++)
    {
        lhs[i] = 0.0;
        x = i%x_fine; // integer mod.
        y = i/x_fine; // integer divide.
        eta1 = 1 - 2*(x%2);

        // + e1.
        lhs[i] = lhs[i]-rhs[y*x_fine+((x+1)%x_fine)];

        // - e1.
        lhs[i] = lhs[i]+ rhs[y*x_fine+((x+x_fine-1)%x_fine)]; // The extra +N is because of the % sign convention.

        // + e2.
        lhs[i] = lhs[i]- eta1*rhs[((y+1)%y_fine)*x_fine+x];

        // - e2.
        lhs[i] = lhs[i]+ eta1*rhs[((y+y_fine-1)%y_fine)*x_fine+x];

        // Normalization.
        lhs[i] = 0.5*lhs[i];

        // 0
        // Added mass term here.
        lhs[i] = lhs[i]+ mass*rhs[i];
    }

}

// Square lattice.
// Kinetic term for a 2D staggered w/ period bc. Applies lhs = A*rhs.
// The unit vectors are e_1 = xhat, e_2 = yhat.
// The "extra_data" is a cast to a complex gauge_field[N*N*2], 
//    loaded by the function read_lattice_u1. 
// Apsi[x][y] = m psi[x,y] - U[y][x,x+1] 
void square_staggered_u1(complex<double>* lhs, complex<double>* rhs, void* extra_data)
{
    // Declare variables.
    int i;
    int x,y;
    double eta1; 
    staggered_u1_op* stagif = (staggered_u1_op*)extra_data;
    complex<double>* lattice = stagif->lattice;
    double mass = stagif->mass; 
    int x_fine = stagif->x_fine;
    int y_fine = stagif->y_fine;

    // For a 2D square lattice, the stencil is:
    //   1 |  0 -eta1  0 |
    //   - | +1    0  -1 |  , where eta1 = (-1)^x
    //   2 |  0 +eta1  0 |
    //
    // e2 = yhat
    // ^
    // | 
    // |-> e1 = xhat

    // Apply the stencil.
    for (i = 0; i < x_fine*y_fine; i++)
    {
        lhs[i] = 0.0;
        x = i%x_fine; // integer mod.
        y = i/x_fine; // integer divide.
        eta1 = 1 - 2*(x%2);

        // + e1.
        lhs[i] = lhs[i]-lattice[y*x_fine*2+x*2]*rhs[y*x_fine+((x+1)%x_fine)];

        // - e1.
        lhs[i] = lhs[i]+ conj(lattice[y*x_fine*2+((x+x_fine-1)%x_fine)*2])*rhs[y*x_fine+((x+x_fine-1)%x_fine)]; // The extra +N is because of the % sign convention.

        // + e2.
        lhs[i] = lhs[i]- eta1*lattice[y*x_fine*2+x*2+1]*rhs[((y+1)%y_fine)*x_fine+x];

        // - e2.
        lhs[i] = lhs[i]+ eta1*conj(lattice[((y+y_fine-1)%y_fine)*x_fine*2+x*2+1])*rhs[((y+y_fine-1)%y_fine)*x_fine+x];

        // Normalization.
        lhs[i] = 0.5*lhs[i];

        // 0
        // Added mass term here.
        lhs[i] = lhs[i]+ mass*rhs[i];

        // Apply a gamma5.
        /*if ((x+y)%2 == 1)
        {
            lhs[i] = -lhs[i];
        }*/
    }
}

// \gamma_5
void gamma_5(complex<double>* lhs, complex<double>* rhs, void* extra_data)
{
    // Declare variables.
    int i;
    int x,y;
    staggered_u1_op* stagif = (staggered_u1_op*)extra_data;
    int x_fine = stagif->x_fine;
    int y_fine = stagif->y_fine;

    // Apply the stencil.
    for (i = 0; i < x_fine*y_fine; i++)
    {
        x = i%x_fine; // integer mod.
        y = i/x_fine; // integer divide.
        
        lhs[i] = ((double)(1 - 2*((x+y)%2)))*rhs[i];
    }
}

// Square \gamma_5 staggered 2d operator w/ u1 function.
void square_staggered_gamma5_u1(complex<double>* lhs, complex<double>* rhs, void* extra_data)
{
    staggered_u1_op* stagif = (staggered_u1_op*)extra_data;
    complex<double>* tmp = new complex<double>[stagif->x_fine*stagif->y_fine];
    
    square_staggered_u1(tmp, rhs, extra_data);
    gamma_5(lhs, tmp, extra_data);
    
    delete[] tmp;
    
}

// Staggered normal equations.
void square_staggered_normal_u1(complex<double>* lhs, complex<double>* rhs, void* extra_data)
{
    staggered_u1_op* stagif = (staggered_u1_op*)extra_data;
    complex<double>* tmp = new complex<double>[stagif->x_fine*stagif->y_fine];
    
    square_staggered_u1(tmp, rhs, extra_data);
    gamma_5(lhs, tmp, extra_data);
    square_staggered_u1(tmp, lhs, extra_data);
    gamma_5(lhs, tmp, extra_data);
    
    delete[] tmp;
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



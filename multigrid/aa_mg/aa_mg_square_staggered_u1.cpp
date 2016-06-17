// This is a sample code that constructs general kinetic terms!

#include <iostream>
#include <iomanip> // to set output precision.
#include <cmath>
#include <string>
#include <sstream>
#include <complex>
#include <random>

#include "generic_inverters.h"
#include "generic_inverters_precond.h"
#include "generic_vector.h"
#include "mg.h"
#include "mg_real.h"
#include "mg_complex.h"
#include "u1_utils.h"

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

// How many GCR iterations do we use?
//#define GEN_NULL_VECTOR_STEP 300
//#define GEN_NULL_VECTOR_REL_RESID 1e-4

//#define AGGREGATE_FOUR

//#define AGGREGATE_EOCONJ

// Are we testing a random gauge rotation?
//#define TEST_RANDOM_GAUGE

// Are we testing a random field?
//#define TEST_RANDOM_FIELD

// Are we loading a gauge field?
#define LOAD_GAUGE_FIELD
// Is it a heatbath field?
#define HEATBATH

// What type of test should we do?
enum mg_test_types
{
    TOP_LEVEL_ONLY = 0,       // only test the top level solver.
    SMOOTHER_ONLY = 1,        // Top level + smoother
    TWO_LEVEL = 2,            // Two level MG 
    THREE_LEVEL = 3           // Three level MG
};

// Square laplace 2d operator w/out u1 function.
void square_laplace(complex<double>* lhs, complex<double>* rhs, void* extra_data);

// Square staggered 2d operator w/out u1 function.
void square_staggered(complex<double>* lhs, complex<double>* rhs, void* extra_data);

// Square staggered 2d operator w/ u1 function.
void square_staggered_u1(complex<double>* lhs, complex<double>* rhs, void* extra_data);  

enum op_type
{
    STAGGERED = 0,
    LAPLACE = 1,
    LAPLACE_NC2 = 2
};

struct staggered_u1_op
{
    complex<double> *lattice;
    double mass;
    int x_fine;
    int y_fine; 
    int Nc; // only relevant for square laplace. 
};

int main(int argc, char** argv)
{  
    // Declare some variables.
    cout << setiosflags(ios::scientific) << setprecision(6);
    int i, j, x, y;
    complex<double> *lattice; // Holds the gauge field.
    complex<double> *lhs, *rhs, *check; // For some Kinetic terms.
    double explicit_resid = 0.0;
    double bnorm = 0.0;
    std::mt19937 generator (1337u); // RNG, 1337u is the seed. 
    inversion_info invif;
    staggered_u1_op stagif;
    
    // Set parameters. 
    
    // What operator are we using? (Laplace is free only.)
    op_type opt = STAGGERED; // STAGGERED, LAPLACE, LAPLACE_NC2

    // What test are we performing?
    mg_test_types my_test = TWO_LEVEL; //THREE_LEVEL; // TWO_LEVEL is the default which won't override anything.
    
    // L_x = L_y = Dimension for a square lattice.
    int square_size = 32; 
    
    // Describe the staggered fermions.
    double MASS = 1e-12;
    
    // Outer Inverter information.
    double outer_precision = 1e-6; 
    int outer_restart = 32; 
    
    // Multigrid information. 
    int n_refine = 1; // 1 = two level V cycle, 2 = three level V cycle, etc. 
    if (my_test == THREE_LEVEL) // FOR TEST ONLY
    {
        n_refine = 2;
    }
    int X_BLOCKSIZE = 8; 
    int Y_BLOCKSIZE = 8;
    int eo = 1; // 0 for no even/odd aggregation, 1 for even/odd aggregation.
    if (opt == LAPLACE || opt == LAPLACE_NC2) // FOR TEST ONLY
    {
        eo = 0;
    }
    
    // Null vector generation
    
// If GEN_NULL_VECTOR isn't defined:
//   1 for just const vector, 2 for const + even/odd vector, 4 for each corner
//    of the hypercube.
// If GEN_NULL_VECTOR is defined and eo = 0:
//   Generate "n_null_vector" null vectors which are block orthogonalized.
// IF GEN_NULL_VECTOR is defined and eo = 1:
//   Generate "n_null_vector" null vectors, partition into even and odd.
//    Total number of null vectors is 2*VECTOR_COUNT. 
    int n_null_vector = 4; // Note: Gets multiplied by 2 for LAPLACE_NC2 test.
    int null_max_iter = 100;
    double null_precision = 1e-4;
    
    // Advanced:
    // IF GEN_NULL_VECTOR is defined and AGGREGATE_EOCONJ is defined:
    //   Generate VECTOR_COUNT null vectors, partition into even and odd, duplicate complex conj.
    //    Total number of null vectors is 4*VECTOR_COUNT. 
    // IF GEN_NULL_VECTOR is defined and AGGREGATE_EO is defined:
    //   Generate VECTOR_COUNT null vectors, partition into corners of hypercube.
    //    Total number of null vectors is 4*VECTOR_COUNT. 
    
    
    // Inner solver.
    inner_solver in_solve = GCR; //GCR; 
    double inner_precision = 1e-3;
    int inner_restart = 10000;
    if (my_test == SMOOTHER_ONLY)
    {
        in_solve = NONE; 
    }
    
    // Smoother
    inner_solver in_smooth = GCR; //NONE; //GCR; 
    double omega_smooth = 0.67; // for minres only. 
    int pre_smooth = 6;
    int post_smooth = 6;
    
    // Gauge field information.
    double BETA = 6.0; // For random gauge field, phase angles have std.dev. 1/sqrt(beta).
                       // For heatbath gauge field, corresponds to non-compact beta.
    
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
    }
    cout << "[OP]: Operator " << op_name << " Mass " << MASS << "\n";
    
    // Only relevant for free laplace test.
    int Nc = 1;  // Only value that matters for staggered
    if (opt == LAPLACE_NC2)
    {
        Nc = 2;
    }
    if (opt == LAPLACE || opt == LAPLACE_NC2) // FOR TEST ONLY
    {
        n_null_vector *= Nc;
    }
    
    // Describe the fine lattice. 
    int x_fine = square_size;
    int y_fine = square_size;
    int fine_size = x_fine*y_fine*Nc;
    
    cout << "[VOL]: X " << x_fine << " Y " << y_fine << " Volume " << x_fine*y_fine;
    if (opt == LAPLACE || opt == LAPLACE_NC2) // FOR TEST ONLY
    {
        cout << " Nc " << Nc;
    }
    cout << "\n";
    
    // Do some allocation.
    // Initialize the lattice. Indexing: index = y*N + x.
    lattice = new complex<double>[2*fine_size];
    lhs = new complex<double>[fine_size];
    rhs = new complex<double>[fine_size];   
    check = new complex<double>[fine_size];   
    // Zero it out.
    zero<double>(lattice, 2*fine_size);
    zero<double>(rhs, fine_size);
    zero<double>(lhs, fine_size);
    zero<double>(check, fine_size);
    //
    
    // Fill stagif.
    stagif.lattice = lattice;
    stagif.mass = MASS; 
    stagif.x_fine = x_fine;
    stagif.y_fine = y_fine; 
    stagif.Nc = Nc; // Only relevant for laplace test only.
    
    
    
    // Describe the gauge field. 
    cout << "[GAUGE]: Creating a gauge field.\n";
    unit_gauge_u1(lattice, x_fine, y_fine);
    
#ifdef TEST_RANDOM_FIELD
    gauss_gauge_u1(lattice, x_fine, y_fine, generator, BETA);
    cout << "[GAUGE]: Created a U(1) gauge field with angle standard deviation " << 1.0/sqrt(BETA) << "\n";
#endif
#ifdef LOAD_GAUGE_FIELD
    if (x_fine == 32 && y_fine == 32)
    {
        if (abs(BETA - 3.0) < 1e-8)
        {
#ifdef HEATBATH
            read_gauge_u1(lattice, x_fine, y_fine, "./cfg/l32t32b30_heatbath.dat");
#else
            read_gauge_u1(lattice, x_fine, y_fine, "./cfg/l32t32b30.dat");
#endif
            cout << "[GAUGE]: Loaded a U(1) gauge field with non-compact beta = " << 1.0/sqrt(BETA) << "\n";
        }
        else if (abs(BETA - 6.0) < 1e-8)
        {
#ifdef HEATBATH
            read_gauge_u1(lattice, x_fine, y_fine, "./cfg/l32t32b60_heatbath.dat");
#else
            read_gauge_u1(lattice, x_fine, y_fine, "./cfg/l32t32b60.dat");
#endif
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
#ifdef HEATBATH
            read_gauge_u1(lattice, x_fine, y_fine, "./cfg/l32t32b60_heatbath.dat");
#else
            read_gauge_u1(lattice, x_fine, y_fine, "./cfg/l32t32b60.dat");
#endif
            cout << "[GAUGE]: Loaded a U(1) gauge field with non-compact beta = " << 1.0/sqrt(BETA) << "\n";
        }
        else if (abs(BETA - 6.0) < 1e-8)
        {
#ifdef HEATBATH
            read_gauge_u1(lattice, x_fine, y_fine, "./cfg/l64t64b60_heatbath.dat");
#else
            read_gauge_u1(lattice, x_fine, y_fine, "./cfg/l64t64b60.dat");
#endif
            cout << "[GAUGE]: Loaded a U(1) gauge field with non-compact beta = " << 1.0/sqrt(BETA) << "\n";
        }
        else
        {
            cout << "[GAUGE]: Saved U(1) gauge field with correct beta, volume does not exist.\n";
        }
    }
    else
    {
        cout << "[GAUGE]: Saved U(1) gauge field with correct beta, volume does not exist.\n";
    }
#endif // LOAD_GAUGE_FIELD
    
#ifdef TEST_RANDOM_GAUGE
    // Generate and perform a random gauge transformation.
    complex<double>* gauge_trans = new complex<double>[fine_size];
    rand_trans_u1(gauge_trans, x_fine, y_fine, generator);
    apply_gauge_trans_u1(lattice, gauge_trans, x_fine, y_fine);
    cout << "[GAUGE]: Performed a random gauge rotation.\n";
#endif
    
    cout << "[GAUGE]: The average plaquette is " << get_plaquette_u1(lattice, x_fine, y_fine) << ".\n";
    
    // Build an mg_struct.
    mg_operator_struct_complex mgstruct;
    mgstruct.x_fine = x_fine;
    mgstruct.y_fine = y_fine; 
    mgstruct.Nc = Nc; // only matters for square laplace.
    mgstruct.n_refine = n_refine; 
    mgstruct.blocksize_x = new int[n_refine];
    mgstruct.blocksize_y = new int[n_refine];
    for (i = 0; i < n_refine; i++)
    {
        mgstruct.blocksize_x[i] = X_BLOCKSIZE;
        mgstruct.blocksize_y[i] = Y_BLOCKSIZE;
    }

#if defined GEN_NULL_VECTOR && (defined AGGREGATE_FOUR || defined AGGREGATE_EOCONJ)
    mgstruct.n_vector = 4*n_null_vector;
    mgstruct.eo = 1;
#elif defined GEN_NULL_VECTOR
    mgstruct.eo = eo;
    mgstruct.n_vector = (eo+1)*n_null_vector;
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

    mgprecond.in_smooth_type = in_smooth; // What inner smoother? MINRES or GCR.
    mgprecond.omega_smooth = omega_smooth; // What relaxation parameter should we use (MINRES only!)
    mgprecond.n_pre_smooth = pre_smooth; // 6 MinRes smoother steps before coarsening.
    mgprecond.n_post_smooth = post_smooth; // 6 MinRes smoother steps after refining.
    mgprecond.in_solve_type = in_solve; // What inner solver? NONE, MINRES, CG, or GCR.
    mgprecond.n_step = inner_restart; // max number of steps to use for inner solver.
    mgprecond.rel_res = inner_precision; // Maximum relative residual for inner solver.
    mgprecond.mgstruct = &mgstruct; // Contains null vectors, fine operator. (Since we don't construct the fine op.)
    mgprecond.coarse_matrix_vector = coarse_square_staggered; // Function which applies the coarse operator. 
    mgprecond.fine_matrix_vector = fine_square_staggered; // Function which applies the fine operator. 
    mgprecond.matrix_extra_data = (void*)&mgstruct; // What extra_data the coarse operator expects. 

    // Set a point on the rhs.
    for (i = 0; i < Nc; i++)
    {
        rhs[(x_fine/2+(y_fine/2)*x_fine)*Nc+i] = 1.0;
    }
    //gaussian<double>(rhs, fine_size, generator);

    // Get norm for rhs.
    bnorm = sqrt(norm2sq<double>(rhs, fine_size));

    // Set a point on the lhs.
    //lhs[x_fine/2+(y_fine/2)*x_fine+1] = 1.0;
    
    // Create a projector.
    mgstruct.null_vectors = new complex<double>**[mgstruct.n_refine];
    // The top level is special since there are no color indices.
    mgstruct.null_vectors[0] = new complex<double>*[mgstruct.n_vector];
    for (j = 0; j < mgstruct.n_vector; j++)
    {
        mgstruct.null_vectors[0][j] = new complex<double>[fine_size];
        zero<double>(mgstruct.null_vectors[0][j], fine_size);
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
        for (i = 0; i < fine_size; i++)
        {
            mgstruct.null_vectors[0][0][i] = 1;
#ifdef TEST_RANDOM_GAUGE
            mgstruct.null_vectors[0][0][i] *= gauge_trans[i];
#endif
        }
    }
    else if (mgstruct.n_vector == 2) // constant, even/odd phase. 
    {
        cout << "[MG]: Null vector 1 is a constant.\n";
        cout << "[MG]: Null vector 2 is an even/odd phase.\n";
        for (i = 0; i < fine_size; i++)
        {
            mgstruct.null_vectors[0][0][i] = 1;
            x = i % x_fine;
            y = i / x_fine;
            mgstruct.null_vectors[0][1][i] = ((x+y)%2 == 0) ? complex<double>(0.0,1.0) : complex<double>(0.0,-1.0);
#ifdef TEST_RANDOM_GAUGE
            mgstruct.null_vectors[0][0][i] *= (gauge_trans[i]);
            mgstruct.null_vectors[0][1][i] *= (gauge_trans[i]);
#endif
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
        for (i = 0; i < fine_size; i++)
        {
            x = i % x_fine;
            y = i / x_fine;
            mgstruct.null_vectors[0][2*(y%2)+(x%2)][i] = 1.0;
            mgstruct.null_vectors[0][2*(y%2)+(x%2)][i] = dist(generator);
            
#ifdef TEST_RANDOM_GAUGE
            mgstruct.null_vectors[0][2*(y%2)+(x%2)][i] *= (gauge_trans[i]);
#endif
        }
    }
    else // invalid.
    {
        cout << "Unless you are generating null vectors, you can only use 1, 2, or 4 null vectors.\n";
        return 0;
    }

    
#ifdef TEST_RANDOM_GAUGE
    delete[] gauge_trans;
#endif
    
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
        complex<double>* rand_guess = new complex<double>[fine_size];
        complex<double>* Arand_guess = new complex<double>[fine_size];
    
        // Temporarily set the mass to zero for the null vector generation. 
        stagif.mass = 0.0;
#if defined GEN_NULL_VECTOR && (defined AGGREGATE_FOUR || defined AGGREGATE_EOCONJ)
        for (i = 0; i < mgstruct.n_vector/4; i++) // Because we partition fourfold afterwards.
#else // GEN_NULL_VECTOR is defined.
        for (i = 0; i < mgstruct.n_vector/(mgstruct.eo+1); i++)
#endif
        {
            gaussian<double>(rand_guess, fine_size, generator);

            zero<double>(Arand_guess, fine_size);
            
            (*op)(Arand_guess, rand_guess, (void*)&stagif);
            
            //square_staggered_u1(Arand_guess, rand_guess, (void*)&stagif);
            for (j = 0; j < fine_size; j++)
            {
               Arand_guess[j] = -Arand_guess[j]; 
            }
            zero<double>(mgstruct.null_vectors[0][i], fine_size);

            minv_vector_gcr(mgstruct.null_vectors[0][i], Arand_guess, fine_size, null_max_iter, null_precision, op, (void*)&stagif); 
            //minv_vector_gcr(mgstruct.null_vectors[0][i], Arand_guess, fine_size, GEN_NULL_VECTOR_STEP, GEN_NULL_VECTOR_REL_RESID, square_staggered_u1, (void*)&stagif); 

            for (j = 0; j < fine_size; j++)
            {
                mgstruct.null_vectors[0][i][j] += rand_guess[j];
            }

            
            normalize(mgstruct.null_vectors[0][i], fine_size); 
        }

        // This causes a segfault related to the RNG when
        // the vector is initialized.
        delete[] rand_guess; 
        delete[] Arand_guess; 

#if defined GEN_NULL_VECTOR && defined AGGREGATE_FOUR
        for (int n = 0; n < mgstruct.n_vector/4; n++)
        {
            for (i = 0; i < fine_size; i++)
            {
                x = i % x_fine;
                y = i / y_fine;
                if (x%2 == 1 && y%2 == 0)
                {
                    mgstruct.null_vectors[0][n+mgstruct.n_vector/4][i] = mgstruct.null_vectors[0][n][i];
                    mgstruct.null_vectors[0][n][i] = 0.0;
                }
                else if (x%2 == 0 && y%2 == 1)
                {
                    mgstruct.null_vectors[0][n+2*mgstruct.n_vector/4][i] = mgstruct.null_vectors[0][n][i];
                    mgstruct.null_vectors[0][n][i] = 0.0;
                }
                else if (x%2 == 1 && y%2 == 1)
                {
                    mgstruct.null_vectors[0][n+3*mgstruct.n_vector/4][i] = mgstruct.null_vectors[0][n][i];
                    mgstruct.null_vectors[0][n][i] = 0.0;
                }
            }
            normalize(mgstruct.null_vectors[0][n], fine_size);
            normalize(mgstruct.null_vectors[0][n+mgstruct.n_vector/2], fine_size);

        }
#elif defined GEN_NULL_VECTOR && defined AGGREGATE_EOCONJ
        for (int n = 0; n < mgstruct.n_vector/4; n++)
        {
            for (i = 0; i < fine_size; i++)
            {
                x = i % x_fine;
                y = i / y_fine;
                if ((x+y)%2 == 1)
                {
                    mgstruct.null_vectors[0][n+mgstruct.n_vector/4][i] = mgstruct.null_vectors[0][n][i];
                    mgstruct.null_vectors[0][n][i] = 0.0;
                }
            }
            normalize(mgstruct.null_vectors[0][n], fine_size);
            normalize(mgstruct.null_vectors[0][n+mgstruct.n_vector/4], fine_size);
            copy<double>(mgstruct.null_vectors[0][n+2*mgstruct.n_vector/4], mgstruct.null_vectors[0][n], fine_size);
            copy<double>(mgstruct.null_vectors[0][n+3*mgstruct.n_vector/4], mgstruct.null_vectors[0][n+mgstruct.n_vector/4], fine_size);
            conj(mgstruct.null_vectors[0][n+2*mgstruct.n_vector/4], fine_size);
            conj(mgstruct.null_vectors[0][n+3*mgstruct.n_vector/4], fine_size);

        }
#else // even/odd
        if (mgstruct.eo == 1)
        {
            for (int n = 0; n < mgstruct.n_vector/2; n++)
            {
                for (i = 0; i < fine_size; i++)
                {
                    x = i % x_fine;
                    y = i / y_fine;
                    if ((x+y)%2 == 1)
                    {
                        mgstruct.null_vectors[0][n+mgstruct.n_vector/2][i] = mgstruct.null_vectors[0][n][i];
                        mgstruct.null_vectors[0][n][i] = 0.0;
                    }
                }
                normalize(mgstruct.null_vectors[0][n], fine_size);
                normalize(mgstruct.null_vectors[0][n+mgstruct.n_vector/2], fine_size);

            }
        }
#endif // defined GEN_NULL_VECTOR && defined AGGREGATE_FOUR
    
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

                    zero<double>(c_Arand_guess, mgstruct.curr_fine_size); 
                    fine_square_staggered(c_Arand_guess, c_rand_guess, (void*)&mgstruct);
                    
                    
                    for (j = 0; j < mgstruct.curr_fine_size; j++)
                    {
                       c_Arand_guess[j] = -c_Arand_guess[j]; 
                    }
                    
                    zero<double>(mgstruct.null_vectors[mgstruct.curr_level][i], mgstruct.curr_fine_size);
                    
                    // Invert!
                    minv_vector_gcr(mgstruct.null_vectors[mgstruct.curr_level][i], c_Arand_guess, mgstruct.curr_fine_size, null_max_iter, null_precision, fine_square_staggered, &mgstruct);
                    

                    for (j = 0; j < mgstruct.curr_fine_size; j++)
                    {
                        mgstruct.null_vectors[mgstruct.curr_level][i][j] += c_rand_guess[j];
                    }


                    normalize(mgstruct.null_vectors[mgstruct.curr_level][i], mgstruct.curr_fine_size); 
                }

                delete[] c_rand_guess; 
                delete[] c_Arand_guess; 
                
                // Aggregate even/odd if we need to!
                if (mgstruct.eo == 1)
                {
                    for (int n = 0; n < mgstruct.n_vector/2; n++)
                    {
                        for (i = 0; i < mgstruct.curr_fine_size; i++)
                        {
                            int c = i % mgstruct.n_vector; // What color index do we have?
                                                           // 0 to mgstruct.n_vector/2-1 is even, else is odd.
                            //int x_coord = (i - c)/mgstruct.n_vector % mgstruct.curr_x_fine;
                            //int y_coord = ((i - c)/mgstruct.n_vector - x_coord)/mgstruct.curr_x_fine;
                            
                            if (c >= mgstruct.n_vector/2)
                            {
                                mgstruct.null_vectors[mgstruct.curr_level][n+mgstruct.n_vector/2][i] = mgstruct.null_vectors[mgstruct.curr_level][n][i];
                                mgstruct.null_vectors[mgstruct.curr_level][n][i] = 0.0;
                            }
                        }
                        normalize(mgstruct.null_vectors[mgstruct.curr_level][n], mgstruct.curr_fine_size);
                        normalize(mgstruct.null_vectors[mgstruct.curr_level][n+mgstruct.n_vector/2], mgstruct.curr_fine_size);

                    }
                }
                
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
        complex<double> *rhs_PdagP = new complex<double>[fine_size];
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
        complex<double>* rhs_PAP_fine = new complex<double>[fine_size];
        zero<double>(rhs_PAP_fine, fine_size);
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
        mgstruct.null_vectors[1] = new complex<double>[fine_size];

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
        complex<double>* rhs_PdagP = new complex<double>[fine_size];
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

        invif = minv_vector_gcr_restart(lhs, rhs, fine_size, 100000, outer_precision, outer_restart, op, (void*)&stagif);
        //invif = minv_vector_gcr_restart(lhs, rhs, fine_size, 100000, outer_precision, outer_restart, square_staggered_u1, (void*)&stagif);

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
        zero<double>(check, fine_size);
        (*op)(check, lhs, (void*)&stagif); 
        //square_staggered_u1(check, lhs, (void*)&stagif);

        explicit_resid = 0.0;
        for (i = 0; i < fine_size; i++)
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
    for (i = 0; i < fine_size; i++)
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
        zero<double>(lhs, fine_size);
        invif = minv_vector_gcr_var_precond_restart(lhs, rhs, fine_size, 10000, outer_precision, outer_restart, op, (void*)&stagif, mg_preconditioner, (void*)&mgprecond, NULL); 
        //invif = minv_vector_gcr_var_precond_restart(lhs, rhs, fine_size, 10000, outer_precision, outer_restart, square_staggered_u1, (void*)&stagif, mg_preconditioner, (void*)&mgprecond); 

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

        explicit_resid = 0.0;
        for (i = 0; i < fine_size; i++)
        {
            explicit_resid += real(conj(rhs[i] - check[i])*(rhs[i] - check[i]));
        }
        explicit_resid = sqrt(explicit_resid)/bnorm;

        printf("[MG]: [check] should equal [rhs]. The relative residual is %15.20e.\n", explicit_resid);
    } // SMOOTHER_ONLY or TWO_LEVEL
    
    // Free the lattice.
    delete[] lattice;
    delete[] lhs;
    delete[] rhs;
    delete[] check;
    
    // Clean up!
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


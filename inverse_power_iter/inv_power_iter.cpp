// Perform inverse power iterations to get eigenvalues. 
// Can also be used to get multiple eigenvalues. A real Lanczos algorithm would be better, but eh.
// Only really works for hermitian positive definite operators. 

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
#include "u1_utils.h"
#include "operators.h"

using namespace std; 

// Define pi.
#define PI 3.141592653589793

// Custom routine to load gauge field.
void internal_load_gauge_u1(complex<double>* lattice, int x_fine, int y_fine, double BETA);

int main(int argc, char** argv)
{  
    // Declare some variables.
    cout << setiosflags(ios::scientific) << setprecision(6);
    int i, j, x, y;
    complex<double> *lattice; // Holds the gauge field.
    complex<double> *lhs, *rhs, *check; // For some Kinetic terms.
    complex<double> *evals, **evecs; // Hold eigenvalues, eigenvectors. 
    complex<double> tmp; 
    double explicit_resid = 0.0;
    double bnorm = 0.0;
    std::mt19937 generator (1337u); // RNG, 1337u is the seed. 
    inversion_info invif;
    staggered_u1_op stagif;
    
    
    // Set parameters. 
    
    // L_x = L_y = Dimension for a square lattice.
    int square_size = 32; // Can be set on command line with --square_size. 
    
    // Describe the staggered fermions.
    double MASS = 1e-2; // for Rayleigh Quotient. 
    
    // Outer Inverter information.
    double outer_precision = 5e-13; 
    int outer_restart = 256; 
    
    // Gauge field information.
    double BETA = 6.0; // For random gauge field, phase angles have std.dev. 1/sqrt(beta).
                       // For heatbath gauge field, corresponds to non-compact beta.
                       // Can be set on command line with --beta.
    
    // Number of eigenvalues to get.
    int n_evals = 1;
    
    // Are we loading a config from a user-defined file?
    char* load_cfg = NULL;
    bool do_load = false;
    
    
    
    /////////////////////////////////////////////
    // Get a few parameters from command line. //
    /////////////////////////////////////////////
    for (i = 1; i < argc; i++)
    {
        if (strcmp(argv[i], "--help") == 0)
        {
            cout << "--beta [3.0, 6.0, 10.0, 10000.0]       (default 6.0)\n";
            cout << "--mass                                 (default 1e-2)\n";
            cout << "--square_size [32, 64, 128]            (default 32)\n";
            cout << "--n-evals [#]                          (default 1)\n";
            cout << "--load-cfg [path]                      (default do not load, overrides beta)\n";
            return 0;
        }
        if (i+1 != argc)
        {
            if (strcmp(argv[i], "--beta") == 0)
            {
                BETA = atof(argv[i+1]);
                i++;
            }
            else if (strcmp(argv[i], "--mass") == 0)
            {
                MASS = atof(argv[i+1]);
                i++;
            }
            else if (strcmp(argv[i], "--square_size") == 0)
            {
                square_size = atoi(argv[i+1]);
                i++;
            }
            else if (strcmp(argv[i], "--n-evals") == 0)
            {
                n_evals = atoi(argv[i+1]);
                i++;
            }
            else if (strcmp(argv[i], "--load-cfg") == 0)
            {
                load_cfg = argv[i+1];
                do_load = true;
                i++;
            }
            else
            {
                cout << "--beta [3.0, 6.0, 10.0, 10000.0]       (default 6.0)\n";
                cout << "--mass                                 (default 1e-2)\n";
                cout << "--square_size [32, 64, 128]            (default 32)\n";
                cout << "--n-evals [#]                          (default 1)\n";
                cout << "--load-cfg [path]                      (default do not load, overrides beta)\n";
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
    op = square_staggered_normal_u1; //square_staggered; //square_staggered_u1; //op = normal_staggered_u1;
    op_name = "Normal Staggered U(1)";
    cout << "[OP]: Operator " << op_name << " Mass " << MASS << "\n";
    
    // Only relevant for free laplace test.
    int Nc = 1;  // Only value that matters for staggered
    
    // Describe the fine lattice. 
    int x_fine = square_size;
    int y_fine = square_size;
    int fine_size = x_fine*y_fine*Nc;
    
    cout << "[VOL]: X " << x_fine << " Y " << y_fine << " Volume " << x_fine*y_fine;
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
    
    // Allocate space for eigenvalues, eigenvectors.
    if (n_evals > 0)
    {
        evals = new complex<double>[n_evals];
        evecs = new complex<double>*[n_evals];
        for (i = 0; i < n_evals; i++)
        {
            evecs[i] = new complex<double>[fine_size];
            zero<double>(evecs[i], fine_size); // Generate some initial vectors. 
        }
    }
    else
    {
        cout << "ERROR! Negative number of eigenvalues requested.\n";
        return 0;
    }
    cout << "[EVAL]: NEval " << n_evals << "\n";
    
    // Fill stagif.
    stagif.lattice = lattice;
    stagif.mass = MASS; 
    stagif.x_fine = x_fine;
    stagif.y_fine = y_fine; 
    stagif.Nc = Nc; // Only relevant for laplace test only.
    
    // Create the verbosity structure.
    inversion_verbose_struct verb;
    verb.verbosity = VERB_SUMMARY;
    verb.verb_prefix = "[L1]: ";
    verb.precond_verbosity = VERB_DETAIL;
    verb.precond_verb_prefix = "Prec ";
    
    // Describe the gauge field. 
    cout << "[GAUGE]: Creating a gauge field.\n";
    unit_gauge_u1(lattice, x_fine, y_fine);
    
    // Load the gauge field.
    if (do_load)
    {
        read_gauge_u1(lattice, x_fine, y_fine, load_cfg);
        cout << "[GAUGE]: Loaded a U(1) gauge field from " << load_cfg << "\n";
    }
    else // various predefined cfgs. 
    {
        internal_load_gauge_u1(lattice, x_fine, y_fine, BETA);
    }
    
    cout << "[GAUGE]: The average plaquette is " << get_plaquette_u1(lattice, x_fine, y_fine) << ", the topology is " << get_topo_u1(lattice, x_fine, y_fine) << ".\n";
    
    // Begin iterating. 
    for (i = 0; i < n_evals; i++)
    {
        // Create a new vector.
        gaussian<double>(rhs, fine_size, generator);
        
        // Orthogonalize it against previous vectors. 
        if (i > 0)
        {
            for (j = 0; j < i; j++)
            {
                orthogonal<double>(rhs, evecs[j], fine_size);
            }
        }
        
        // Print norm.
        cout << "[NORM]: " << norm2sq<double>(rhs, fine_size) << "\n";
        
        // Normalize. 
        normalize<double>(rhs, fine_size);
        
        // Zero out lhs.
        zero<double>(lhs, fine_size); 

        // Seed some initial guesses for the eigenvectors. 
        complex<double> curr, prev;
        curr = 1.0;
        prev = 0.5;

        int counter = 0;
        while (abs(1.0/curr - 1.0/prev)/abs(1.0/prev) > 1e-6)
        {

            prev = curr;

            // Perform inverse iteration.
            invif = minv_vector_cg(lhs, rhs, fine_size, 500, outer_precision, op, (void*)&stagif, &verb);
            //invif = minv_vector_gcr(lhs, rhs, fine_size, 500, outer_precision, op, (void*)&stagif, &verb);

            // Update Lambda.
            curr = dot<double>(rhs, lhs, fine_size)/norm2sq<double>(rhs, fine_size);

            // Reorthogonalize against previous vectors to avoid spurious re-introduction of eigenvalues.
            if (i > 0)
            {
                for (j = 0; j < i; j++)
                {
                    orthogonal<double>(lhs, evecs[j], fine_size);
                }
            }
            normalize<double>(lhs, fine_size);
            copy<double>(rhs, lhs, fine_size);
            zero<double>(lhs, fine_size);

            // Print
            cout << "[EVAL]: Val " << i << " Iteration " << counter++ << " Guess " << 1.0/curr << " RelResid " << abs(1.0/curr - 1.0/prev)/abs(1.0/prev) << "\n";
        }
        
        // Copy vector back into eigenvectors array, save eigenvalue.
        copy<double>(evecs[i], rhs, fine_size);
        evals[i] = 1.0/curr; 
    }
        
    // Sort and print eigenvalues. 
    for (i = 0; i < n_evals; i++)
    {
        for (j = 0; j < n_evals; j++)
        {
            if (real(evals[i]) < real(evals[j])) { tmp = evals[i]; evals[i] = evals[j]; evals[j] = tmp; } 
        }
    }
    
    cout << "\n\nAll eigenvalues:\n";
    for (i = 0; i < n_evals; i++)
    {
        cout << "[EVAL]: Val " << i << " Guess " << evals[i] << "\n";
    }
    
    // Free the lattice.
    delete[] lattice;
    delete[] lhs;
    delete[] rhs;
    delete[] check;
    
    for (i = 0; i < n_evals; i++)
    {
        delete[] evecs[i];
    }
    delete[] evecs;
    delete[] evals; 
    
    
    return 0; 
}

// Custom routine to load gauge field.
void internal_load_gauge_u1(complex<double>* lattice, int x_fine, int y_fine, double BETA)
{
    if (x_fine == 32 && y_fine == 32)
    {
        if (abs(BETA - 3.0) < 1e-8)
        {
            read_gauge_u1(lattice, x_fine, y_fine, "../multigrid/aa_mg/cfg/l32t32b30_heatbath.dat");
            cout << "[GAUGE]: Loaded a U(1) gauge field with non-compact beta = " << 1.0/sqrt(BETA) << "\n";
        }
        else if (abs(BETA - 6.0) < 1e-8)
        {
            read_gauge_u1(lattice, x_fine, y_fine, "../multigrid/aa_mg/cfg/l32t32b60_heatbath.dat");
            cout << "[GAUGE]: Loaded a U(1) gauge field with non-compact beta = " << 1.0/sqrt(BETA) << "\n";
        }
        else if (abs(BETA - 10.0) < 1e-8)
        {
            read_gauge_u1(lattice, x_fine, y_fine, "../multigrid/aa_mg/cfg/l32t32b100_heatbath.dat");
            cout << "[GAUGE]: Loaded a U(1) gauge field with non-compact beta = " << 1.0/sqrt(BETA) << "\n";
        }
        else if (abs(BETA - 10000.0) < 1e-8)
        {
            read_gauge_u1(lattice, x_fine, y_fine, "../multigrid/aa_mg/cfg/l32t32bperturb_heatbath.dat");
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
            read_gauge_u1(lattice, x_fine, y_fine, "../multigrid/aa_mg/cfg/l64t64b30_heatbath.dat");
            cout << "[GAUGE]: Loaded a U(1) gauge field with non-compact beta = " << 1.0/sqrt(BETA) << "\n";
        }
        else if (abs(BETA - 6.0) < 1e-8)
        {
            read_gauge_u1(lattice, x_fine, y_fine, "../multigrid/aa_mg/cfg/l64t64b60_heatbath.dat");
            cout << "[GAUGE]: Loaded a U(1) gauge field with non-compact beta = " << 1.0/sqrt(BETA) << "\n";
        }
        else if (abs(BETA - 10.0) < 1e-8)
        {
            read_gauge_u1(lattice, x_fine, y_fine, "../multigrid/aa_mg/cfg/l64t64b100_heatbath.dat");
            cout << "[GAUGE]: Loaded a U(1) gauge field with non-compact beta = " << 1.0/sqrt(BETA) << "\n";
        }
        else if (abs(BETA - 10000.0) < 1e-8)
        {
            read_gauge_u1(lattice, x_fine, y_fine, "../multigrid/aa_mg/cfg/l64t64bperturb_heatbath.dat");
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
            read_gauge_u1(lattice, x_fine, y_fine, "../multigrid/aa_mg/cfg/l128t128b60_heatbath.dat");
            cout << "[GAUGE]: Loaded a U(1) gauge field with non-compact beta = " << 1.0/sqrt(BETA) << "\n";
        }
        else if (abs(BETA - 10.0) < 1e-8)
        {
            read_gauge_u1(lattice, x_fine, y_fine, "../multigrid/aa_mg/cfg/l128t128b100_heatbath.dat");
            cout << "[GAUGE]: Loaded a U(1) gauge field with non-compact beta = " << 1.0/sqrt(BETA) << "\n";
        }
        else if (abs(BETA - 10000.0) < 1e-8)
        {
            read_gauge_u1(lattice, x_fine, y_fine, "../multigrid/aa_mg/cfg/l128t128bperturb_heatbath.dat");
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

}
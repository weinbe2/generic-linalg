// Calculate level crossings for the hermitian index staggered operator.
// Uses inverse power iterations to get the eigenvalues of the operator squared,
// then figure out if the proper eigenvalue is +/- sqrt(eigen).

#include <iostream>
#include <iomanip> // to set output precision.
#include <cmath>
#include <string>
#include <sstream>
#include <complex>
#include <random>
#include <cstring> // should be replaced by using sstream

#include "arpack_interface.h"
#include "generic_vector.h"
#include "verbosity.h"
#include "u1_utils.h"
#include "operators.h"

using namespace std; 

// Define pi.
#define PI 3.141592653589793

enum op_type
{
    STAGGERED = 1,
    LAPLACE = 0,
    G5_STAGGERED = 2,
    STAGGERED_NORMAL = 3,
    STAGGERED_INDEX = 4
};

// Custom routine to load gauge field.
void internal_load_gauge_u1(complex<double>* lattice, int x_fine, int y_fine, double BETA);

int main(int argc, char** argv)
{  
    // Declare some variables.
    cout << setiosflags(ios::scientific) << setprecision(6);
    int i,j;
    complex<double> *lattice; // Holds the gauge field.
    complex<double> *lhs, *rhs, *check; // For some Kinetic terms.
    complex<double> *evals, **evecs; // Hold eigenvalues, eigenvectors. 
    complex<double> tmp; 
    complex<double>* tmp2; 
    std::mt19937 generator (1337u); // RNG, 1337u is the seed. 
    inversion_info invif;
    staggered_u1_op stagif;
    
    
    // Set parameters. 
    
    // L_x = L_y = Dimension for a square lattice.
    int square_size = 32; // Can be set on command line with --square_size. 
    
    // Describe the staggered fermions.
    double MASS = 1e-2; // for Rayleigh Quotient. 
    
    // Eigensolver precision.
    double precision = 1e-7; 
    
    // What operator should we use?
    op_type opt = STAGGERED; // STAGGERED, LAPLACE, G5_STAGGERED, STAGGERED_NORMAL
    
    // Gauge field information.
    double BETA = 6.0; // For random gauge field, phase angles have std.dev. 1/sqrt(beta).
                       // For heatbath gauge field, corresponds to non-compact beta.
                       // Can be set on command line with --beta.
    
    // Number of eigenvalues to get.
    int n_evals = 1;
    
    // Number of internal values to get. By default n_evals*2.5.
    int n_cv = -1;
    
    // Load an external cfg?
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
            cout << "--operator [laplace, staggered, g5_staggered, \n";
            cout << "       normal_staggered, index]        (default staggered)\n";
            cout << "--mass                                 (default 1e-2)\n";
            cout << "--lattice-size [32, 64, 128]           (default 32)\n";
            cout << "--n-evals [#]                          (default 1)\n";
            cout << "--n-internal-vals [#]                  (default 2.5*n-evals)\n";
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
            else if (strcmp(argv[i], "--operator") == 0)
            {
                if (strcmp(argv[i+1], "laplace") == 0)
                {
                    opt = LAPLACE; 
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
            else if (strcmp(argv[i], "--n-internal-vals") == 0)
            {
                n_cv = atoi(argv[i+1]);
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
                cout << "--operator [laplace, staggered,\n";
                cout << "       g5_staggered, normal_staggered] (default staggered)\n";
                cout << "--mass                                 (default 1e-2)\n";
                cout << "--lattice-size [32, 64, 128]           (default 32)\n";
                cout << "--n-evals [#]                          (default 1)\n";
                cout << "--n-internal-vals [#]                  (default 2.5*n-evals)\n";
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
    
    // Check that number of internal vals is sane.
    if (n_cv == -1 || n_cv < n_evals)
    {
        n_cv = (2*n_evals + n_evals/2);
    }
    cout << "[EVAL]: NEval " << n_evals << " Internal vals " << n_cv << "\n";
    
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
    
    cout << "[GAUGE]: The average plaquette is " << get_plaquette_u1(lattice, x_fine, y_fine) << ".\n";
    
    // Okay! Time to pull in the arpack bindings. 
    
    // Make an arpack struct.
    arpack_dcn_t* ar_strc = arpack_dcn_init(fine_size, n_evals, n_cv); // max eigenvectors, max internal vectors
    char eigtype[3]; strcpy(eigtype, "SM"); // This could be SM (smallest magnitude).
                 // This could also be LM (largest magnitude)
                 // SR (smallest real), SI (smallest imaginary),
                 // and similar for largest.
    arpack_solve_t info_solve = arpack_dcn_getev(ar_strc, evals, evecs, fine_size, n_evals, n_cv, 4000, eigtype, precision, 0.0, op, (void*)&stagif); 
    
    // Print info about the eigensolve.
    cout << "[ARPACK]: Number of converged eigenvalues: " << info_solve.nconv << "\n";
    cout << "[ARPACK]: Number of iteration steps: " << info_solve.niter << "\n";
    cout << "[ARPACK]: Number of matrix multiplies: " << info_solve.nops << "\n";
    
    // End of arpack bindings!
    
    // Sort eigenvalues (done differently depending on the operator).
    for (i = 0; i < n_evals; i++)
    {
        for (j = 0; j < n_evals-1; j++)
        {
            switch (opt)
            {
                case STAGGERED:
                    if (imag(evals[j]) > imag(evals[j+1]))
                    {
                        tmp = evals[j]; evals[j] = evals[j+1]; evals[j+1] = tmp;
                        tmp2 = evecs[j]; evecs[j] = evecs[j+1]; evecs[j+1] = tmp2;
                    }
                    break;
                case LAPLACE:
                case G5_STAGGERED:
                case STAGGERED_NORMAL:
                case STAGGERED_INDEX:
                    if (real(evals[j]) > real(evals[j+1]))
                    {
                        tmp = evals[j]; evals[j] = evals[j+1]; evals[j+1] = tmp;
                        tmp2 = evecs[j]; evecs[j] = evecs[j+1]; evecs[j+1] = tmp2;
                    }
                    break; 
            }
        }
    }
                    
    
    cout << "\n\nAll eigenvalues:\n";
    for (i = 0; i < n_evals; i++)
    {
        cout << "[FINEVAL]: Mass " << MASS << " Num " << i << " Eval " << evals[i] << "\n";
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

// Hermitian index operator squared. 
void staggered_index_operator_sq(complex<double>* lhs, complex<double>* rhs, void* extra_data)
{
    staggered_u1_op* stagif = (staggered_u1_op*)extra_data;
    complex<double>* tmp = new complex<double>[stagif->x_fine*stagif->y_fine];
    
    staggered_index_operator(tmp, rhs, extra_data);
    staggered_index_operator(lhs, tmp, extra_data);
    
    delete[] tmp;
    
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
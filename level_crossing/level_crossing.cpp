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

// Form the index operator squared.
void staggered_index_operator_sq(complex<double>* lhs, complex<double>* rhs, void* extra_data);

// A shifted versions of the index op for rayleigh quotients.
void shift_op(complex<double>* lhs, complex<double>* rhs, void* extra_data);
struct shift_struct
{
    void* extra_data;
    void (*op)(complex<double>* a, complex<double>* b, void* extra_data);
    complex<double> shift;
    int length; 
};

int main(int argc, char** argv)
{  
    // Declare some variables.
    cout << setiosflags(ios::scientific) << setprecision(6);
    int i, j, x, y;
    complex<double> *lattice; // Holds the gauge field.
    complex<double> *lhs, *rhs, *check; // For some Kinetic terms.
    complex<double> *evals, *index_evals, **evecs; // Hold eigenvalues, eigenvectors. 
    complex<double> tmp; 
    double explicit_resid = 0.0;
    double p_bnorm = 0.0;
    double m_bnorm = 0.0;
    double relnorm = 0.0;
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
            else
            {
                cout << "--beta [3.0, 6.0, 10.0, 10000.0]       (default 6.0)\n";
                cout << "--mass                                 (default 1e-2)\n";
                cout << "--square_size [32, 64, 128]            (default 32)\n";
                cout << "--n-evals [#]                          (default 1)\n";
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
    op_name = "Squared staggered index operator.";
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
        index_evals = new complex<double>[n_evals];
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
    internal_load_gauge_u1(lattice, x_fine, y_fine, BETA);
    
    cout << "[GAUGE]: The average plaquette is " << get_plaquette_u1(lattice, x_fine, y_fine) << ".\n";
    
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
        while ((relnorm = abs(1.0/curr - 1.0/prev)/abs(1.0/prev)) > 1e-10)
        {

            prev = curr;

            // Switch over to rayleigh quotients at some point.
            if (counter < 10 || relnorm > 1e-4)  // inverse iterations.
            {
                // Perform inverse iteration.
                invif = minv_vector_cg(lhs, rhs, fine_size, 500, outer_precision,  staggered_index_operator_sq, (void*)&stagif, &verb);
                //invif = minv_vector_gcr(lhs, rhs, fine_size, 500, outer_precision, op, (void*)&stagif, &verb);
                
                curr = dot<double>(rhs, lhs, fine_size)/norm2sq<double>(rhs, fine_size);
            }
            else // rayleigh
            {
                shift_struct sif;
                sif.extra_data = (void*)&stagif;
                sif.op = staggered_index_operator_sq;
                sif.length = fine_size;
                sif.shift = 1.0/curr;
                
                invif = minv_vector_cg(lhs, rhs, fine_size, 500, outer_precision, shift_op, (void*)&sif, &verb);
                
                normalize<double>(lhs, fine_size);
                staggered_index_operator_sq(check, lhs, (void*)&stagif);
                curr = 1.0/(dot<double>(lhs, check, fine_size)/norm2sq<double>(lhs, fine_size));
            }

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
        
        // Okay, moment of truth: find if + or - is the right eigenvalue of the original unsquared op. 
        
        // First, try plus.
        zero<double>(lhs, fine_size);
        staggered_index_operator(lhs, evecs[i], (void*)&stagif);
        p_bnorm = 0.0;
        for (j = 0; j < fine_size; j++)
        {
            p_bnorm += real(conj(lhs[j] - sqrt(evals[i])*evecs[i][j])*(lhs[j] - sqrt(evals[i])*evecs[i][j]));
        }
        cout << "[EVAL]: Plus " << p_bnorm << " Minus ";
        
        // Next, try minus.
        m_bnorm = 0.0;
        for (j = 0; j < fine_size; j++)
        {
            m_bnorm += real(conj(lhs[j] + sqrt(evals[i])*evecs[i][j])*(lhs[j] + sqrt(evals[i])*evecs[i][j]));
        }
        cout << m_bnorm << "\n";
        
        index_evals[i] = (p_bnorm < m_bnorm) ? sqrt(evals[i]) : (-sqrt(evals[i]));
        
    }
        
    // Sort and print eigenvalues. 
    for (i = 0; i < n_evals; i++)
    {
        for (j = 0; j < n_evals; j++)
        {
            if (real(index_evals[i]) < real(index_evals[j]))
            {
                tmp = evals[i]; evals[i] = evals[j]; evals[j] = tmp;
                tmp = index_evals[i]; index_evals[i] = index_evals[j]; index_evals[j] = tmp;
            }
        }
    }
    
    cout << "\n\nAll eigenvalues:\n";
    for (i = 0; i < n_evals; i++)
    {
        cout << "[FINEVAL]: Mass " << MASS << " Num " << i << " SquareEval " << real(evals[i]) << " IndexEval " << real(index_evals[i]) << "\n";
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
    delete[] index_evals; 
    
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
    
// A shifted versions of the index op for rayleigh quotients.
void shift_op(complex<double>* lhs, complex<double>* rhs, void* extra_data)
{
    shift_struct* sop = (shift_struct*)extra_data;
    (*sop->op)(lhs, rhs, sop->extra_data);
    for (int i = 0; i < sop->length; i++)
    {
        lhs[i] -= sop->shift*rhs[i];
    }
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
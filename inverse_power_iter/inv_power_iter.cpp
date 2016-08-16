// Measure pions!

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

using namespace std; 

// Define pi.
#define PI 3.141592653589793

// Square staggered 2d operator w/ u1 function.
void square_staggered_u1(complex<double>* lhs, complex<double>* rhs, void* extra_data);  

void normal_staggered_u1(complex<double>* lhs, complex<double>* rhs, void* extra_data);  

struct staggered_u1_op
{
    complex<double> *lattice;
    double mass; // for Rayleigh quotient iteration
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
        }
    }
    
    //printf("Mass %.8e Blocksize %d %d Null Vectors %d\n", MASS, X_BLOCKSIZE, Y_BLOCKSIZE, n_null_vector);
    //return 0;
    
    ///////////////////////////////////////
    // End of human-readable parameters! //
    ///////////////////////////////////////
    
    
    string op_name;
    void (*op)(complex<double>*, complex<double>*, void*);
    op = square_staggered_u1; //op = normal_staggered_u1;
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

    
    cout << "[GAUGE]: The average plaquette is " << get_plaquette_u1(lattice, x_fine, y_fine) << ".\n";
    
    // Guess an initial vector.
    for (i = 0; i < fine_size; i++)
    {
        rhs[i] = 1.0 + 1e-4*std::polar(1.0, 0.5*(i%4)*PI);
    }
    normalize<double>(rhs, fine_size);
    
    complex<double> curr, prev;
    curr = 1.0;
    prev = 0.5;
    
    int counter = 0;
    while (abs(1.0/curr - 1.0/prev)/abs(1.0/prev) > 1e-12)
    {
        
        prev = curr;
        
        // Perform inverse iteration.
        //invif = minv_vector_cg(lhs, rhs, fine_size, 500, outer_precision, op, (void*)&stagif, &verb);
        invif = minv_vector_bicgstab(lhs, rhs, fine_size, 500, outer_precision, op, (void*)&stagif, &verb);
        
        // Update Lambda.
        curr = dot<double>(rhs, lhs, fine_size)/norm2sq<double>(rhs, fine_size);
        
        // Copy over.
        normalize<double>(lhs, fine_size);
        copy<double>(rhs, lhs, fine_size);
        
        // Print
        cout << "[EVAL]: " << counter++ << " " << 1.0/curr << " " << abs(1.0/curr - 1.0/prev)/abs(1.0/prev) << "\n";
        
    }
    
    // Free the lattice.
    delete[] lattice;
    delete[] lhs;
    delete[] rhs;
    delete[] check;
    
    
    return 0; 
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
    double mass = stagif->mass;  // for Rayleigh. 
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

// Implement the normal operator, taking advantage of gamma 5 hermiticity.
void normal_staggered_u1(complex<double>* lhs, complex<double>* rhs, void* extra_data)
{
    // Declare variables.
    int i;
    int x,y;
    double phase;
    staggered_u1_op* stagif = (staggered_u1_op*)extra_data;
    
    int fine_size = stagif->x_fine*stagif->y_fine*stagif->Nc;
    complex<double>* temp = new complex<double>[fine_size];
    
    // Apply D
    zero<double>(temp, fine_size);
    square_staggered_u1(temp, rhs, extra_data);
    
    // Apply gamma 5.
    for (i = 0; i < fine_size; i++)
    {
        x = i%stagif->x_fine; // integer mod.
        y = i/stagif->x_fine; // integer divide.
        
        temp[i] = (2.0*((x+y)%2)-1.0)*temp[i];
    }
    
    // Apply D
    zero<double>(lhs, fine_size);
    square_staggered_u1(lhs, temp, extra_data);
    
    // Apply gamma 5.
    for (i = 0; i < fine_size; i++)
    {
        x = i%stagif->x_fine; // integer mod.
        y = i/stagif->x_fine; // integer divide.
        
        lhs[i] = (2.0*((x+y)%2)-1.0)*lhs[i];
    }
}


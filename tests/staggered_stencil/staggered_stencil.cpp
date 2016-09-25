// Test file to compare BiCGStab with BiCGStab-l, as well as restarted versions. 

#include <iostream>
#include <iomanip> // to set output precision.
#include <cmath>
#include <string>
#include <sstream>
#include <complex>
#include <random>
#include <cstring> // should be replaced by using sstream

// Things for timing. 
#include <stdlib.h>
#include <unistd.h>
#include <time.h>


#include "generic_vector.h"
#include "verbosity.h"
#include "u1_utils.h"
#include "operators.h"
#include "operators_stencil.h"
#include "lattice.h"
#include "coarse_stencil.h"


using namespace std; 

// Define pi.
#define PI 3.141592653589793

// Custom routine to load gauge field.
void internal_load_gauge_u1(complex<double>* lattice, int x_fine, int y_fine, double BETA);

// Timing routine.
timespec diff(timespec start, timespec end)
{
    timespec tmp;
    if ((end.tv_nsec-start.tv_nsec) < 0)
    {
        tmp.tv_sec = end.tv_sec-start.tv_sec-1;
        tmp.tv_nsec = 1000000000+end.tv_nsec-start.tv_nsec;
    }
    else
    {
        tmp.tv_sec = end.tv_sec-start.tv_sec;
        tmp.tv_nsec = end.tv_nsec-start.tv_nsec;
    }
    return tmp;
}

int main(int argc, char** argv)
{  
    // Declare some variables.
    cout << setiosflags(ios::scientific) << setprecision(6);
    int i;
    complex<double> *lattice; // Holds the gauge field.
    complex<double> *lhs, *rhs, *check; // For some Kinetic terms.
    complex<double> tmp; 
    complex<double>* tmp2; 
    std::mt19937 generator (1337u); // RNG, 1337u is the seed. 
    staggered_u1_op stagif;
    double bnorm; 
    stringstream sstream; 
    
    // Timing.
    timespec time1, time2, timediff;

    // Set parameters. 
    
    // L_x = L_y = Dimension for a square lattice.
    int square_size = 64; // Can be set on command line with --square_size. 
    
    // What masses should we use?
    double mass = 1e-2;
    
    
    // Gauge field information.
    double BETA = 6.0; // For random gauge field, phase angles have std.dev. 1/sqrt(beta).
                       // For heatbath gauge field, corresponds to non-compact beta.
                       // Can be set on command line with --beta.
    
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
            cout << "--lattice-size [32, 64, 128]           (default 64)\n";
            cout << "--mass [#]                             (default 1e-2)\n";
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
            else if (strcmp(argv[i], "--lattice-size") == 0)
            {
                square_size = atoi(argv[i+1]);
                i++;
            }
            else if (strcmp(argv[i], "--mass") == 0)
            {
                mass = atof(argv[i+1]);
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
                cout << "--lattice-size [32, 64, 128]           (default 32)\n";
                cout << "--mass [#]                             (default 1e-2)\n";
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
    
    // Only relevant for free laplace test.
    int Nc = 1;  // Only value that matters for staggered
    
    // Describe the fine lattice. 
    int x_fine = square_size;
    int y_fine = square_size;
    int fine_size = x_fine*y_fine*Nc;
    
    // Create a lattice objet.
    int lattice_size[2];
    lattice_size[0] = x_fine;
    lattice_size[1] = y_fine;
    
    Lattice Lat(2, lattice_size, Nc); 
    
    cout << "[VOL]: X " << x_fine << " Y " << y_fine << " Volume " << x_fine*y_fine;
    cout << "\n";
    
    // Do some allocation.
    // Initialize the lattice. Indexing: index = y*N + x.
    lattice = new complex<double>[2*fine_size];
    lhs = new complex<double>[fine_size];
    rhs = new complex<double>[fine_size];   
    check = new complex<double>[fine_size];   
    tmp2 = new complex<double>[fine_size];
    // Zero it out.
    zero<double>(lattice, 2*fine_size);
    zero<double>(rhs, fine_size);
    zero<double>(lhs, fine_size);
    zero<double>(check, fine_size);
    zero<double>(tmp2, fine_size);
    //
    
    // Fill stagif.
    stagif.lattice = lattice;
    stagif.mass = mass;  
    stagif.x_fine = x_fine;
    stagif.y_fine = y_fine; 
    stagif.Nc = Nc; // Only relevant for laplace test.
    
    // Create the verbosity structure.
    inversion_verbose_struct verb;
    verb.verbosity = VERB_DETAIL;
    verb.verb_prefix = "[BICGSTAB]: ";
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
    cout << "[GAUGE]: The topological charge is " << get_topo_u1(lattice, x_fine, y_fine) << ".\n";
    
    // Compare applying the staggered operator to getting the staggered stencil and applying it using the
    // generic stencil interface.
    
    // Set a random source. Should be the same every time b/c we hard code the seed above.
    //gaussian<double>(rhs, fine_size, generator); 
    
    zero<double>(rhs, fine_size);
    rhs[x_fine+1] = 1.0;
    
    // Apply the staggered operator.
    cout << "[L1]: Applying staggered operator function.\n";
    zero<double>(lhs, fine_size);
    square_staggered_u1(lhs, rhs, (void*)&stagif);
    cout << "clov\t+x\t+y\t-x\t-y\n";
    cout << lhs[x_fine+1] << "\t" << lhs[x_fine+2] << "\t" << lhs[2*x_fine+1] << "\t" << lhs[x_fine] << "\t" << lhs[1] << "\n";
    
    // Get and apply a stencil.
    cout << "[L1]: Applying staggered operator stencil.\n";
    stencil_2d* stenc = new stencil_2d(&Lat, get_stencil_size(STAGGERED));
    get_square_staggered_u1_stencil(stenc, &stagif);
    zero<double>(tmp2, fine_size);
    apply_stencil_2d(tmp2, rhs, (void*)stenc);
    cout << "clov\t+x\t+y\t-x\t-y\n";
    cout << tmp2[x_fine+1] << "\t" << tmp2[x_fine+2] << "\t" << tmp2[2*x_fine+1] << "\t" << tmp2[x_fine] << "\t" << tmp2[1] << "\n";
    
    // Compare.
    cout << "The difference between the two methods is " << diffnorm2sq<double>(lhs, tmp2, fine_size) << ".\n\n";
    
    // Apply the g5_staggered operator.
    cout << "[L1]: Applying g5_staggered operator function.\n";
    zero<double>(lhs, fine_size);
    square_staggered_gamma5_u1(lhs, rhs, (void*)&stagif);
    cout << "clov\t+x\t+y\t-x\t-y\n";
    cout << lhs[x_fine+1] << "\t" << lhs[x_fine+2] << "\t" << lhs[2*x_fine+1] << "\t" << lhs[x_fine] << "\t" << lhs[1] << "\n";
    
    // Get and apply a stencil.
    cout << "[L1]: Applying g5_staggered operator stencil.\n";
    stenc->clear_stencils();
    get_square_staggered_gamma5_u1_stencil(stenc, &stagif);
    zero<double>(tmp2, fine_size);
    apply_stencil_2d(tmp2, rhs, (void*)stenc);
    cout << "clov\t+x\t+y\t-x\t-y\n";
    cout << tmp2[x_fine+1] << "\t" << tmp2[x_fine+2] << "\t" << tmp2[2*x_fine+1] << "\t" << tmp2[x_fine] << "\t" << tmp2[1] << "\n";
    
    // Compare.
    cout << "The difference between the two methods is " << diffnorm2sq<double>(lhs, tmp2, fine_size) << ".\n";
    
    // Free the lattice.
    delete stenc;
    delete[] lattice;
    delete[] lhs;
    delete[] rhs;
    delete[] check;
    delete[] tmp2; 
    
    
    return 0; 
}


// Custom routine to load gauge field.
void internal_load_gauge_u1(complex<double>* lattice, int x_fine, int y_fine, double BETA)
{
    if (x_fine == 32 && y_fine == 32)
    {
        if (abs(BETA - 3.0) < 1e-8)
        {
            read_gauge_u1(lattice, x_fine, y_fine, "../../multigrid/aa_mg/cfg/l32t32b30_heatbath.dat");
            cout << "[GAUGE]: Loaded a U(1) gauge field with non-compact beta = " << 1.0/sqrt(BETA) << "\n";
        }
        else if (abs(BETA - 6.0) < 1e-8)
        {
            read_gauge_u1(lattice, x_fine, y_fine, "../../multigrid/aa_mg/cfg/l32t32b60_heatbath.dat");
            cout << "[GAUGE]: Loaded a U(1) gauge field with non-compact beta = " << 1.0/sqrt(BETA) << "\n";
        }
        else if (abs(BETA - 10.0) < 1e-8)
        {
            read_gauge_u1(lattice, x_fine, y_fine, "../../multigrid/aa_mg/cfg/l32t32b100_heatbath.dat");
            cout << "[GAUGE]: Loaded a U(1) gauge field with non-compact beta = " << 1.0/sqrt(BETA) << "\n";
        }
        else if (abs(BETA - 10000.0) < 1e-8)
        {
            read_gauge_u1(lattice, x_fine, y_fine, "../../multigrid/aa_mg/cfg/l32t32bperturb_heatbath.dat");
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
            read_gauge_u1(lattice, x_fine, y_fine, "../../multigrid/aa_mg/cfg/l64t64b30_heatbath.dat");
            cout << "[GAUGE]: Loaded a U(1) gauge field with non-compact beta = " << 1.0/sqrt(BETA) << "\n";
        }
        else if (abs(BETA - 6.0) < 1e-8)
        {
            read_gauge_u1(lattice, x_fine, y_fine, "../../multigrid/aa_mg/cfg/l64t64b60_heatbath.dat");
            cout << "[GAUGE]: Loaded a U(1) gauge field with non-compact beta = " << 1.0/sqrt(BETA) << "\n";
        }
        else if (abs(BETA - 10.0) < 1e-8)
        {
            read_gauge_u1(lattice, x_fine, y_fine, "../../multigrid/aa_mg/cfg/l64t64b100_heatbath.dat");
            cout << "[GAUGE]: Loaded a U(1) gauge field with non-compact beta = " << 1.0/sqrt(BETA) << "\n";
        }
        else if (abs(BETA - 10000.0) < 1e-8)
        {
            read_gauge_u1(lattice, x_fine, y_fine, "../../multigrid/aa_mg/cfg/l64t64bperturb_heatbath.dat");
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
            read_gauge_u1(lattice, x_fine, y_fine, "../../multigrid/aa_mg/cfg/l128t128b60_heatbath.dat");
            cout << "[GAUGE]: Loaded a U(1) gauge field with non-compact beta = " << 1.0/sqrt(BETA) << "\n";
        }
        else if (abs(BETA - 10.0) < 1e-8)
        {
            read_gauge_u1(lattice, x_fine, y_fine, "../../multigrid/aa_mg/cfg/l128t128b100_heatbath.dat");
            cout << "[GAUGE]: Loaded a U(1) gauge field with non-compact beta = " << 1.0/sqrt(BETA) << "\n";
        }
        else if (abs(BETA - 10000.0) < 1e-8)
        {
            read_gauge_u1(lattice, x_fine, y_fine, "../../multigrid/aa_mg/cfg/l128t128bperturb_heatbath.dat");
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


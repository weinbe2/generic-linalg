// This test verifies that the various staggered Dslash functions are working correctly.
// * Compares the \gamma5 D function with applying D then \gamma_5.
// * Compares the D^\dagger function with applying \gamma_5 then D then \gamma_5
// * Compares summing applying D_{eo}, D_{oe}, and the mass term to applying D.
// * Compares summing applying -D_{eo}D_{oe}, -D_{oe}D_{eo}, and the mass term squared to applying D^\dagger D.
// * Compares inverting D directly with even-odd reconstruction.

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


#include "generic_bicgstab_l.h"
#include "generic_cg.h"
#include "generic_vector.h"
#include "verbosity.h"
#include "u1_utils.h"
#include "operators.h"

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
    int i,j,x,y;
    complex<double> *lattice; // Holds the gauge field.
    complex<double> *lhs, *lhs2, *rhs, *rhs2, *check; // For some Kinetic terms.
    complex<double> tmp; 
    complex<double> *tmp2, *tmp3; 
    std::mt19937 generator (1337u); // RNG, 1337u is the seed. 
    inversion_info invif;
    staggered_u1_op stagif;
    double bnorm; 
    stringstream sstream; 
    
    // Timing.
    timespec time1, time2, timediff;

    // Set parameters. 
    
    // L_x = L_y = Dimension for a square lattice.
    int square_size = 64; // Can be set on command line with --square_size. 
    
    // What masses should we use?
    double mass = 1e-1;
    
    // Solver precision.
    double outer_precision = 1e-6; 
    
    // Restart iterations.
    int outer_restart = 64;
    
    
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
            cout << "--mass [#]                             (default 1e-1)\n";
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
                cout << "--operator [laplace, staggered, g5_staggered, \n";
                cout << "       normal_staggered, index]        (default staggered)\n";
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
    
    cout << "[VOL]: X " << x_fine << " Y " << y_fine << " Volume " << x_fine*y_fine;
    cout << "\n";
    
    // Do some allocation.
    // Initialize the lattice. Indexing: index = y*N + x.
    lattice = new complex<double>[2*fine_size];
    lhs = new complex<double>[fine_size];
    rhs = new complex<double>[fine_size];   
    rhs2 = new complex<double>[fine_size];   
    lhs2 = new complex<double>[fine_size];
    check = new complex<double>[fine_size];   
    tmp2 = new complex<double>[fine_size];
    tmp3 = new complex<double>[fine_size];
    // Zero it out.
    zero<double>(lattice, 2*fine_size);
    zero<double>(rhs, fine_size);
    zero<double>(rhs2, fine_size);
    zero<double>(lhs, fine_size);
    zero<double>(lhs2, fine_size);
    zero<double>(check, fine_size);
    zero<double>(tmp2, fine_size);
    zero<double>(tmp3, fine_size); 
    //
    
    // Fill stagif.
    stagif.lattice = lattice;
    stagif.mass = mass;  
    stagif.x_fine = x_fine;
    stagif.y_fine = y_fine; 
    stagif.Nc = Nc; // Only relevant for laplace test only.
    
    // Create the verbosity structure.
    inversion_verbose_struct verb;
    verb.verbosity = VERB_SUMMARY; //VERB_DETAIL;
    
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
    
    
    /**************
    * BEGIN TESTS *
    **************/
    
    // Set a random source. Should be the same every time b/c we hard code the seed above.
    gaussian<double>(rhs, fine_size, generator); 
    
    ////////////////////
    // * Compares the \gamma5 D function with applying D then \gamma_5.
    ////////////////////
    
    // Apply \gamma_5 D
    square_staggered_gamma5_u1(lhs, rhs, (void*)&stagif);
    
    // Apply D then \gamma_5
    square_staggered_u1(tmp2, rhs, (void*)&stagif);
    gamma_5(lhs2, tmp2, (void*)&stagif);
    
    cout << "[TEST1]: \\gamma_5 D test: relative residual " << diffnorm2sq(lhs, lhs2, fine_size) << ".\n";
    
    ////////////////////
    // * Compares the D^\dagger function with applying \gamma_5 then D then \gamma_5
    ////////////////////
    
    // Apply D^\dagger
    square_staggered_dagger_u1(lhs, rhs, (void*)&stagif);
    
    // Apply \gamma_5 then D then \gamma_5
    gamma_5(lhs2, rhs, (void*)&stagif);
    square_staggered_u1(tmp2, lhs2, (void*)&stagif);
    gamma_5(lhs2, tmp2, (void*)&stagif);
    
    cout << "[TEST2]: D^\\dagger test:  relative residual " << diffnorm2sq(lhs, lhs2, fine_size) << ".\n";
    
    ////////////////////
    // * Compares summing applying D_{eo}, D_{oe}, and the mass term to applying D.
    ////////////////////
    
    // Apply D
    square_staggered_u1(lhs, rhs, (void*)&stagif);
    
    // Apply D_{eo}, D_{oe}, m in pieces.
    square_staggered_deo_u1(lhs2, rhs, (void*)&stagif); // D_{eo} piece in lhs2.
    square_staggered_doe_u1(tmp2, rhs, (void*)&stagif); // D_{oe} piece in tmp2.
    for (i = 0; i < fine_size; i++)
    {
        lhs2[i] = lhs2[i] + tmp2[i] + stagif.mass*rhs[i]; // combine, add mass.
    }
    
    cout << "[TEST3]: D_{eo}+D_{oe}+m test:  relative residual " << diffnorm2sq(lhs, lhs2, fine_size) << ".\n";
    
    ////////////////////
    // * Compares summing applying -D_{eo}D_{oe}, -D_{oe}D_{eo}, and the mass term squared to applying D^\dagger D.
    ////////////////////
    
     // Apply D^\dagger D
    square_staggered_normal_u1(lhs, rhs, (void*)&stagif);
    
    // Apply D_{eo}D_{oe}.
    square_staggered_doe_u1(tmp2, rhs, (void*)&stagif);
    square_staggered_deo_u1(lhs2, tmp2, (void*)&stagif);
    
    // Apply D_{oe}D_{eo}.
    square_staggered_deo_u1(tmp2, rhs, (void*)&stagif);
    square_staggered_doe_u1(tmp3, tmp2, (void*)&stagif);
    
    // Combine into m^2 - D_{eo}D_{oe} - D_{oe}D_{eo}.
    for (i = 0; i < fine_size; i++)
    {
        lhs2[i] = stagif.mass*stagif.mass*rhs[i] - lhs2[i] - tmp3[i];
    }
    
    cout << "[TEST4]: D^\\dagger D test:  relative residual " << diffnorm2sq(lhs, lhs2, fine_size) << ".\n";
    
    
    ////////////////////
    // * Compares inverting D directly with even-odd reconstruction.
    ////////////////////
    
    // Invert D directly using BiCGstab-L.
    verb.verb_prefix = "[D]: ";
    invif = minv_vector_bicgstab_l(lhs, rhs, fine_size, 100000, outer_precision, 4, square_staggered_u1, (void*)&stagif, &verb);
    
    // Invert D using an even/odd decomposition, solving the even system with CG, reconstructing odd.
    
    verb.verb_prefix = "[D_EO_PRE]: ";
    // Prepare rhs: m rhs_e - D_{eo} rhs_o
    square_staggered_deo_u1(rhs2, rhs, (void*)&stagif);
    for (i = 0; i < fine_size; i++)
    {
        x = i%x_fine;
        y = i/x_fine;
        if ((x+y)%2 == 0) // even
        {
            rhs2[i] = stagif.mass*rhs[i] - rhs2[i];
        }
    }
    
    // Perform even/odd inversion
    invif = minv_vector_cg(lhs2, rhs2, fine_size, 1000000, outer_precision, square_staggered_m2mdeodoe_u1, (void*)&stagif, &verb);
    
    // Reconstruct odd: m^{-1}*(rhs_o - D_{oe} lhs2_e)
    square_staggered_doe_u1(tmp2, lhs2, (void*)&stagif);
    for (i = 0; i < fine_size; i++)
    {
        x = i%x_fine;
        y = i/x_fine;
        if ((x+y)%2 == 1) // odd
        {
            lhs2[i] = 1.0/stagif.mass*(rhs[i] - tmp2[i]);
        }
    }
    
    cout << "[TEST5]: Even/odd precond test:  relative residual " << diffnorm2sq(lhs, lhs2, fine_size) << ".\n";
    
    
    /************
    * END TESTS *
    ************/
    
    // Free the lattice.
    delete[] lattice;
    delete[] lhs;
    delete[] lhs2; 
    delete[] rhs;
    delete[] rhs2;
    delete[] check;
    delete[] tmp2;
    delete[] tmp3; 
    
    
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


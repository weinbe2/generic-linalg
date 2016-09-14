// Test file to compare GCR and CGNE of the staggered operator
// on an even-only source.  

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


#include "generic_cg.h"
#include "generic_gmres.h"
#include "generic_gcr.h"
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
    complex<double> *lhs, *rhs, *check; // For some Kinetic terms.
    complex<double> tmp; 
    complex<double>* tmp2; 
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
    double mass = 1e-2;
    
    // Solver precision.
    double outer_precision = 1e-10; 
    
    // Restart iterations.
    int outer_restart = 64;
    
    // What operator should we use?
    op_type opt = STAGGERED; 
    
    // Are we solving an even or an odd source only?
    bool is_even = true; 
    
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
            cout << "--source-type [even, odd]              (default even)\n";
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
            else if (strcmp(argv[i], "--source-type") == 0)
            {
                if (strcmp(argv[i+1], "even") == 0)
                {
                    is_even = true;
                }
                else if (strcmp(argv[i+1], "odd") == 0)
                {
                    is_even = false; 
                }
                i++;
            }
            else
            {
                cout << "--beta [3.0, 6.0, 10.0, 10000.0]       (default 6.0)\n";
                cout << "--lattice-size [32, 64, 128]           (default 32)\n";
                cout << "--mass [#]                             (default 1e-2)\n";
                cout << "--load-cfg [path]                      (default do not load, overrides beta)\n";
                cout << "--source-type [even, odd]              (default even)\n";
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
    op = square_staggered_u1;
    op_name = "Staggered U(1)";

    cout << "[OP]: Operator " << op_name << " Mass " << mass << "\n";
    
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
    stagif.Nc = Nc; // Only relevant for laplace test only.
    
    // Create the verbosity structure.
    inversion_verbose_struct verb;
    verb.verbosity = VERB_DETAIL;
    verb.verb_prefix = "[GCR]: ";
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
    
    // Set a random source. Should be the same every time b/c we hard code the seed above.
    gaussian<double>(rhs, fine_size, generator); 

    // Project onto even or odd only.
    gamma_5(tmp2, rhs, (void*)&stagif);
    double factor = (is_even) ? 1.0 : -1.0;
    for (i = 0; i < fine_size; i++)
    {
         rhs[i] = 0.5*(rhs[i]+factor*tmp2[i]);
    }
    
    // Test GCR(\infty)
    verb.verb_prefix = "[GCR]: ";
        
    zero<double>(lhs, fine_size);
    clock_gettime(CLOCK_PROCESS_CPUTIME_ID, &time1);
    invif = minv_vector_gcr(lhs, rhs, fine_size, 100000, outer_precision, op, (void*)&stagif, &verb); 
    clock_gettime(CLOCK_PROCESS_CPUTIME_ID, &time2); timediff = diff(time1, time2); cout << "Time " << ((double)(long int)(1000000000*timediff.tv_sec + timediff.tv_nsec))*1e-9 << "\n";

    // Test GMRES(\infty)
    zero<double>(lhs, fine_size);
    verb.verb_prefix = "[GMRES]: ";
        
    zero<double>(lhs, fine_size);
    clock_gettime(CLOCK_PROCESS_CPUTIME_ID, &time1);
    invif = minv_vector_gmres(lhs, rhs, fine_size, 100000, outer_precision, op, (void*)&stagif, &verb);
    clock_gettime(CLOCK_PROCESS_CPUTIME_ID, &time2); timediff = diff(time1, time2); cout << "Time " << ((double)(long int)(1000000000*timediff.tv_sec + timediff.tv_nsec))*1e-9 << "\n";
    
    // Test CGNE with applying D^\dag to rhs.
    verb.verb_prefix = "[CG_NORMAL]: ";
        
    zero<double>(lhs, fine_size);
    // Apply D^\dagger
    zero<double>(check, fine_size); gamma_5(check, rhs, (void*)&stagif);
    zero<double>(tmp2, fine_size); square_staggered_u1(tmp2, check, (void*)&stagif);
    zero<double>(check, fine_size); gamma_5(check, tmp2, (void*)&stagif);
        
    clock_gettime(CLOCK_PROCESS_CPUTIME_ID, &time1);
    invif = minv_vector_cg(lhs, check, fine_size, 100000, outer_precision, square_staggered_normal_u1, (void*)&stagif, &verb); // Remark: Dslash count should be doubled. 
    clock_gettime(CLOCK_PROCESS_CPUTIME_ID, &time2); timediff = diff(time1, time2); cout << "Time " << ((double)(long int)(1000000000*timediff.tv_sec + timediff.tv_nsec))*1e-9 << "\n";

    // Test CGNE withOUT applying D^\dag to rhs.
    verb.verb_prefix = "[CG_NORMAL_NODDAG]: ";
        
    zero<double>(lhs, fine_size);
        
    clock_gettime(CLOCK_PROCESS_CPUTIME_ID, &time1);
    invif = minv_vector_cg(lhs, rhs, fine_size, 100000, outer_precision, square_staggered_normal_u1, (void*)&stagif, &verb); // Remark: Dslash count should be doubled. 
    clock_gettime(CLOCK_PROCESS_CPUTIME_ID, &time2); timediff = diff(time1, time2); cout << "Time " << ((double)(long int)(1000000000*timediff.tv_sec + timediff.tv_nsec))*1e-9 << "\n";


    
    // Free the lattice.
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


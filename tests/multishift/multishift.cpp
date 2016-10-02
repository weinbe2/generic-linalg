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

#include "generic_bicgstab_m.h"
#include "generic_cg_m.h"
#include "generic_cr_m.h"
#include "generic_bicgstab.h"
#include "generic_cg.h"
#include "generic_cr.h"
#include "generic_vector.h"
#include "verbosity.h"
#include "u1_utils.h"
#include "operators.h"

using namespace std; 

// Define pi.
#define PI 3.141592653589793


// Custom routine to load gauge field.
void internal_load_gauge_u1(complex<double>* lattice, int x_fine, int y_fine, double BETA);

// Real free laplace.
void square_laplace(double* lhs, double* rhs, void* extra_data);

// Real free staggered.
void square_staggered(double* lhs, double* rhs, void* extra_data);

// Gamma 5
void gamma_5(double* lhs, double* rhs, void* extra_data);

// Real free g5 staggered.
void square_staggered_gamma5(double* lhs, double* rhs, void* extra_data);

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
    int i,n;
    complex<double> *lattice; // Holds the gauge field.
    complex<double> tmp; 
    std::mt19937 generator (1337u); // RNG, 1337u is the seed. 
    inversion_info invif, invif2; 
    staggered_u1_op stagif;
    double bnorm; 
    stringstream sstream; 
    int seq_inversion, multi_inversion; // keep track of inversions. 
    
    // Timing.
    timespec time1, time2, timediff;
    double seq_time, shift_time; 

    // Set parameters. 
    
    // L_x = L_y = Dimension for a square lattice.
    int square_size = 64; // Can be set on command line with --square_size. 
    
    // What masses should we use?
    double mass = 1e-2;
    
    // Solver precision.
    double outer_precision = 1e-10; 
    
    // Restart iterations.
    int outer_restart = 64;
    
    // Residual check frequency in CG-M.
    int resid_check_freq = 10;
    
    // What operator should we use?
    op_type opt = STAGGERED; // STAGGERED, LAPLACE, G5_STAGGERED, STAGGERED_NORMAL
    
    
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
            cout << "--operator [laplace, staggered, g5_staggered, \n";
            cout << "       normal_staggered, index]        (default staggered)\n";
            cout << "--lattice-size [32, 64, 128]           (default 64)\n";
            cout << "--mass [#]                             (default 1e-2)\n";
            cout << "--resid-check-freq [#]                 (default 10)\n";
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
            else if (strcmp(argv[i], "--resid-check-freq") == 0)
            {
                resid_check_freq = atoi(argv[i+1]);
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
                cout << "--resid-check-freq [#]                 (default 10)\n";
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
            op_name = "Laplace U(1)";
            op = square_laplace_u1;
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
    // Zero it out.
    zero<double>(lattice, 2*fine_size);
    //
    
    // Fill stagif.
    stagif.lattice = lattice;
    stagif.mass = mass;  
    stagif.x_fine = x_fine;
    stagif.y_fine = y_fine; 
    stagif.Nc = Nc; // Only relevant for laplace test only.
    
    // Create the verbosity structure.
    inversion_verbose_struct verb;
    verb.verbosity = VERB_SUMMARY;
    verb.verb_prefix = "[CG-M]: ";
    verb.precond_verbosity = VERB_SUMMARY;
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
    
    
    
    // Make some sources.
    int n_shift = 7;
    stagif.mass = 0.001;
    double* shifts = new double[n_shift];
    
    // Perform a real test of CG-M.
    cout << "\n[TEST]: Test CG-M on the free laplace for masses 0.001, 0.005, 0.01, 0.05, 0.1, 0.5, 1.\n";
    
    shifts[0] = 0.001;
    //shifts[1] = 0.005;
    //shifts[2] = 0.01;
    shifts[1] = 0.01;
    shifts[2] = 0.005;
    shifts[3] = 0.05;
    shifts[4] = 0.1;
    shifts[5] = 0.5;
    shifts[6] = 1.0;
    
    double** lhs_real = new double*[n_shift];
    for (n = 0; n < n_shift; n++)
    {
        lhs_real[n] = new double[fine_size];
        zero<double>(lhs_real[n], fine_size);
    }
    
    
    // Set a random source. Should be the same every time b/c we hard code the seed above.
    double* rhs_real = new double[fine_size];
    gaussian<double>(rhs_real, fine_size, generator); 
    bnorm = sqrt(norm2sq<double>(rhs_real, fine_size));
    
    // CG!
    
    // Compare against sequential inversions.
    seq_inversion = 0;
    
    // Heavy -> lightest.
    zero<double>(lhs_real[0], fine_size);
    clock_gettime(CLOCK_PROCESS_CPUTIME_ID, &time1);
    for (n = n_shift - 1; n >= 0; n--)
    {
        stagif.mass = shifts[n];
        sstream.str(string()); sstream << "[CG-REAL_mass=" << stagif.mass << "]: "; verb.verb_prefix = sstream.str(); // Update the verbosity string.
        invif = minv_vector_cg(lhs_real[0], rhs_real, fine_size, 10000, outer_precision, square_laplace, (void*)&stagif, &verb);
        seq_inversion += invif.ops_count;
    }
    clock_gettime(CLOCK_PROCESS_CPUTIME_ID, &time2); timediff = diff(time1, time2);
    seq_time = ((double)(long int)(1000000000*timediff.tv_sec + timediff.tv_nsec))*1e-9;
    
    
    // Now do multishift. Get ready!
    verb.verb_prefix = "[CG-M-REAL]: ";
    stagif.mass = shifts[0];
    for (n = n_shift - 1; n >= 0; n--)
    {
        shifts[n] -= shifts[0];
    }
    zero<double>(lhs_real[0], fine_size);
    
    clock_gettime(CLOCK_PROCESS_CPUTIME_ID, &time1);
    invif = minv_vector_cg_m(lhs_real, rhs_real, n_shift, fine_size, resid_check_freq, 10000, outer_precision, shifts, square_laplace, (void*)&stagif, true, &verb);
    multi_inversion = invif.ops_count; 
    // Check for failures, fix.
    for (n = 0; n < n_shift; n++)
    {
        if (sqrt(invif.resSqmrhs[n]) > outer_precision*bnorm)
        {
            stagif.mass += shifts[n];
            verb.verb_prefix = "[CG-M-REAL-FIX]: ";
            invif2 = minv_vector_cg(lhs_real[n], rhs_real, fine_size, 10000, outer_precision, square_laplace, (void*)&stagif, &verb);
            multi_inversion += invif2.ops_count; 
            stagif.mass -= shifts[n];
        }
    }
    clock_gettime(CLOCK_PROCESS_CPUTIME_ID, &time2); timediff = diff(time1, time2);
    shift_time = ((double)(long int)(1000000000*timediff.tv_sec + timediff.tv_nsec))*1e-9;
    
    cout << "[TEST]: CG-REAL Multirhs took " << multi_inversion << " matrix ops and " << shift_time << " seconds, sequential inversions took " << seq_inversion << " matrix ops and " << seq_time << " seconds.\n";
    
    
    // CR!
    cout << "\n[TEST]: Test CG-R on the free g5 staggered for masses 0.001, 0.005, 0.01, 0.05, 0.1, 0.5, 1.\n";
    
    // Compare against sequential inversions.
    seq_inversion = 0;
    
    // Reset masses. 
    shifts[0] = 0.001;
    //shifts[1] = 0.005;
    //shifts[2] = 0.01;
    shifts[1] = 0.01;
    shifts[2] = 0.005;
    shifts[3] = 0.05;
    shifts[4] = 0.1;
    shifts[5] = 0.5;
    shifts[6] = 1.0;
    
    // Heavy -> lightest.
    zero<double>(lhs_real[0], fine_size);
    clock_gettime(CLOCK_PROCESS_CPUTIME_ID, &time1);
    for (n = n_shift - 1; n >= 0; n--)
    {
        stagif.mass = shifts[n];
        sstream.str(string()); sstream << "[CR-REAL_mass=" << stagif.mass << "]: "; verb.verb_prefix = sstream.str(); // Update the verbosity string.
        invif = minv_vector_cr(lhs_real[0], rhs_real, fine_size, 10000, outer_precision, square_staggered_gamma5, (void*)&stagif, &verb);
        seq_inversion += invif.ops_count;
    }
    clock_gettime(CLOCK_PROCESS_CPUTIME_ID, &time2); timediff = diff(time1, time2);
    seq_time = ((double)(long int)(1000000000*timediff.tv_sec + timediff.tv_nsec))*1e-9;
    
    
    // Now do multishift. Get ready!
    verb.verb_prefix = "[CR-M-REAL]: ";
    stagif.mass = shifts[0];
    for (n = n_shift - 1; n >= 0; n--)
    {
        shifts[n] -= shifts[0];
    }
    zero<double>(lhs_real[0], fine_size);
    
    clock_gettime(CLOCK_PROCESS_CPUTIME_ID, &time1);
    invif = minv_vector_cr_m(lhs_real, rhs_real, n_shift, fine_size, resid_check_freq, 10000, outer_precision, shifts, square_staggered_gamma5, (void*)&stagif, true, &verb);
    multi_inversion = invif.ops_count; 
    // Check for failures, fix.
    for (n = 0; n < n_shift; n++)
    {
        if (sqrt(invif.resSqmrhs[n]) > outer_precision*bnorm)
        {
            stagif.mass += shifts[n];
            verb.verb_prefix = "[CR-M-REAL-FIX]: ";
            invif2 = minv_vector_cr(lhs_real[n], rhs_real, fine_size, 10000, outer_precision, square_staggered_gamma5, (void*)&stagif, &verb);
            multi_inversion += invif2.ops_count; 
            stagif.mass -= shifts[n];
        }
    }
    clock_gettime(CLOCK_PROCESS_CPUTIME_ID, &time2); timediff = diff(time1, time2);
    shift_time = ((double)(long int)(1000000000*timediff.tv_sec + timediff.tv_nsec))*1e-9;
    
    cout << "[TEST]: CR-REAL Multirhs took " << multi_inversion << " matrix ops and " << shift_time << " seconds, sequential inversions took " << seq_inversion << " matrix ops and " << seq_time << " seconds.\n";
    
    // BICGSTAB!
    cout << "\n[TEST]: Test BICGSTAB-M on the free staggered for masses 0.1, 0.2, 0.3, 0.4, 0.5, 0.8, 1.0\n";
    
    // Compare against sequential inversions.
    seq_inversion = 0;
    
    // Reset masses. 
    shifts[0] = 0.1;
    shifts[1] = 0.2;
    shifts[2] = 0.3;
    shifts[3] = 0.4;
    shifts[4] = 0.5;
    shifts[5] = 0.8;
    shifts[6] = 1.0;
    
    // Heavy -> lightest.
    zero<double>(lhs_real[0], fine_size);
    clock_gettime(CLOCK_PROCESS_CPUTIME_ID, &time1);
    for (n = n_shift - 1; n >= 0; n--)
    {
        stagif.mass = shifts[n];
        sstream.str(string()); sstream << "[BICGSTAB-REAL_mass=" << stagif.mass << "]: "; verb.verb_prefix = sstream.str(); // Update the verbosity string.
        invif = minv_vector_bicgstab(lhs_real[0], rhs_real, fine_size, 10000, outer_precision, square_staggered, (void*)&stagif, &verb);
        seq_inversion += invif.ops_count;
    }
    clock_gettime(CLOCK_PROCESS_CPUTIME_ID, &time2); timediff = diff(time1, time2);
    seq_time = ((double)(long int)(1000000000*timediff.tv_sec + timediff.tv_nsec))*1e-9;
    
    
    // Now do multishift. Get ready!
    verb.verb_prefix = "[BICGSTAB-M-REAL]: ";
    stagif.mass = shifts[0];
    for (n = n_shift - 1; n >= 0; n--)
    {
        shifts[n] -= shifts[0];
    }
    zero<double>(lhs_real[0], fine_size);
    
    clock_gettime(CLOCK_PROCESS_CPUTIME_ID, &time1);
    invif = minv_vector_bicgstab_m(lhs_real, rhs_real, n_shift, fine_size, resid_check_freq, 10000, outer_precision, shifts, square_staggered, (void*)&stagif, true, &verb);
    multi_inversion = invif.ops_count; 
    // Check for failures, fix.
    for (n = 0; n < n_shift; n++)
    {
        if (sqrt(invif.resSqmrhs[n]) > outer_precision*bnorm)
        {
            stagif.mass += shifts[n];
            verb.verb_prefix = "[BICGSTAB-M-REAL-FIX]: ";
            invif2 = minv_vector_bicgstab(lhs_real[n], rhs_real, fine_size, 10000, outer_precision, square_staggered, (void*)&stagif, &verb);
            multi_inversion += invif2.ops_count; 
            stagif.mass -= shifts[n];
        }
    }
    clock_gettime(CLOCK_PROCESS_CPUTIME_ID, &time2); timediff = diff(time1, time2);
    shift_time = ((double)(long int)(1000000000*timediff.tv_sec + timediff.tv_nsec))*1e-9;
    
    cout << "[TEST]: BICGSTAB-REAL Multirhs took " << multi_inversion << " matrix ops and " << shift_time << " seconds, sequential inversions took " << seq_inversion << " matrix ops and " << seq_time << " seconds.\n";
    
    for (n = 0; n < n_shift; n++)
    {
        delete[] lhs_real[n];
    }
    delete[] lhs_real;
    delete[] rhs_real; 
    stagif.mass = mass;  
    
    
    // Perform a complex test of CG-M.
    cout << "\n[TEST]: Test CG-M on the beta=6.0 gauge laplace for masses 0.001, 0.005, 0.01, 0.05, 0.1, 0.5, 1.\n";
    
    
    shifts[0] = 0.001;
    //shifts[1] = 0.005;
    //shifts[2] = 0.01;
    shifts[1] = 0.01;
    shifts[2] = 0.005;
    shifts[3] = 0.05;
    shifts[4] = 0.1;
    shifts[5] = 0.5;
    shifts[6] = 1.0;
    
    complex<double>** lhs_cplx = new complex<double>*[n_shift];
    for (n = 0; n < n_shift; n++)
    {
        lhs_cplx[n] = new complex<double>[fine_size];
        zero<double>(lhs_cplx[n], fine_size);
    }
    
    
    // Set a random source. Should be the same every time b/c we hard code the seed above.
    complex<double>* rhs_cplx = new complex<double>[fine_size];
    gaussian<double>(rhs_cplx, fine_size, generator); 
    bnorm = sqrt(norm2sq<double>(rhs_cplx, fine_size));
    
    // Compare against sequential inversions.
    seq_inversion = 0;
    
    // Heavy -> lightest.
    zero<double>(lhs_cplx[0], fine_size);
    clock_gettime(CLOCK_PROCESS_CPUTIME_ID, &time1);
    for (n = n_shift - 1; n >= 0; n--)
    {
        stagif.mass = shifts[n];
        sstream.str(string()); sstream << "[CG-CPLX_mass=" << stagif.mass << "]: "; verb.verb_prefix = sstream.str(); // Update the verbosity string.
        invif = minv_vector_cg(lhs_cplx[0], rhs_cplx, fine_size, 10000, outer_precision, square_laplace_u1, (void*)&stagif, &verb);
        seq_inversion += invif.ops_count;
    }
    clock_gettime(CLOCK_PROCESS_CPUTIME_ID, &time2); timediff = diff(time1, time2);
    seq_time = ((double)(long int)(1000000000*timediff.tv_sec + timediff.tv_nsec))*1e-9;
    
    
    // Now do multishift. Get ready!
    verb.verb_prefix = "[CG-CPLX-M]: ";
    stagif.mass = shifts[0];
    for (n = n_shift - 1; n >= 0; n--)
    {
        shifts[n] -= shifts[0];
    }
    zero<double>(lhs_cplx[0], fine_size);
    
    clock_gettime(CLOCK_PROCESS_CPUTIME_ID, &time1);
    invif = minv_vector_cg_m(lhs_cplx, rhs_cplx, n_shift, fine_size, resid_check_freq, 10000, outer_precision, shifts, square_laplace_u1, (void*)&stagif, true, &verb);
    multi_inversion = invif.ops_count; 
    // Check for failures, fix.
    for (n = 0; n < n_shift; n++)
    {
        if (sqrt(invif.resSqmrhs[n]) > outer_precision*bnorm)
        {
            stagif.mass += shifts[n];
            verb.verb_prefix = "[CG-M-CPLX-FIX]: ";
            invif2 = minv_vector_cg(lhs_cplx[n], rhs_cplx, fine_size, 10000, outer_precision, square_laplace_u1, (void*)&stagif, &verb);
            multi_inversion += invif2.ops_count; 
            stagif.mass -= shifts[n];
        }
    }
    clock_gettime(CLOCK_PROCESS_CPUTIME_ID, &time2); timediff = diff(time1, time2);
    shift_time = ((double)(long int)(1000000000*timediff.tv_sec + timediff.tv_nsec))*1e-9;
    
    cout << "[TEST]: CG-CPLX Multirhs took " << multi_inversion << " matrix ops and " << shift_time << " seconds, sequential inversions took " << seq_inversion << " matrix ops and " << seq_time << " seconds.\n";
    
    
    // CR-CPLX
    cout << "\n[TEST]: Test CR-M on the beta=6.0 g5 staggered for masses 0.001, 0.005, 0.01, 0.05, 0.1, 0.5, 1.\n";
    
    // Compare against sequential inversions.
    seq_inversion = 0;
    
    // Reset masses. 
    shifts[0] = 0.001;
    //shifts[1] = 0.005;
    //shifts[2] = 0.01;
    shifts[1] = 0.01;
    shifts[2] = 0.005;
    shifts[3] = 0.05;
    shifts[4] = 0.1;
    shifts[5] = 0.5;
    shifts[6] = 1.0;
    
    // Heavy -> lightest.
    zero<double>(lhs_cplx[0], fine_size);
    clock_gettime(CLOCK_PROCESS_CPUTIME_ID, &time1);
    for (n = n_shift - 1; n >= 0; n--)
    {
        stagif.mass = shifts[n];
        sstream.str(string()); sstream << "[CR-CPLX_mass=" << stagif.mass << "]: "; verb.verb_prefix = sstream.str(); // Update the verbosity string.
        invif = minv_vector_cr(lhs_cplx[0], rhs_cplx, fine_size, 10000, outer_precision, square_staggered_gamma5_u1, (void*)&stagif, &verb);
        seq_inversion += invif.ops_count;
    }
    clock_gettime(CLOCK_PROCESS_CPUTIME_ID, &time2); timediff = diff(time1, time2);
    seq_time = ((double)(long int)(1000000000*timediff.tv_sec + timediff.tv_nsec))*1e-9;
    
    
    // Now do multishift. Get ready!
    verb.verb_prefix = "[CR-M-CPLX]: ";
    stagif.mass = shifts[0];
    for (n = n_shift - 1; n >= 0; n--)
    {
        shifts[n] -= shifts[0];
    }
    zero<double>(lhs_cplx[0], fine_size);
    
    clock_gettime(CLOCK_PROCESS_CPUTIME_ID, &time1);
    invif = minv_vector_cr_m(lhs_cplx, rhs_cplx, n_shift, fine_size, resid_check_freq, 10000, outer_precision, shifts, square_staggered_gamma5_u1, (void*)&stagif, true, &verb);
    multi_inversion = invif.ops_count; 
    // Check for failures, fix.
    for (n = 0; n < n_shift; n++)
    {
        if (sqrt(invif.resSqmrhs[n]) > outer_precision*bnorm)
        {
            stagif.mass += shifts[n];
            verb.verb_prefix = "[CR-M-CPLX-FIX]: ";
            invif2 = minv_vector_cr(lhs_cplx[n], rhs_cplx, fine_size, 10000, outer_precision, square_staggered_gamma5_u1, (void*)&stagif, &verb);
            multi_inversion += invif2.ops_count; 
            stagif.mass -= shifts[n];
        }
    }
    clock_gettime(CLOCK_PROCESS_CPUTIME_ID, &time2); timediff = diff(time1, time2);
    shift_time = ((double)(long int)(1000000000*timediff.tv_sec + timediff.tv_nsec))*1e-9;
    
    cout << "[TEST]: CR-CPLX Multirhs took " << multi_inversion << " matrix ops and " << shift_time << " seconds, sequential inversions took " << seq_inversion << " matrix ops and " << seq_time << " seconds.\n";
    
    // BICGSTAB!
    cout << "\n[TEST]: Test BICGSTAB-M on the free staggered for masses 0.1, 0.2, 0.3, 0.4, 0.5, 0.8, 1.0\n";
    
    // Compare against sequential inversions.
    seq_inversion = 0;
    
    // Reset masses. 
    shifts[0] = 0.1;
    shifts[1] = 0.2;
    shifts[2] = 0.3;
    shifts[3] = 0.4;
    shifts[4] = 0.5;
    shifts[5] = 0.8;
    shifts[6] = 1.0;
    
    // Heavy -> lightest.
    zero<double>(lhs_cplx[0], fine_size);
    clock_gettime(CLOCK_PROCESS_CPUTIME_ID, &time1);
    for (n = n_shift - 1; n >= 0; n--)
    {
        stagif.mass = shifts[n];
        sstream.str(string()); sstream << "[BICGSTAB-CPLX_mass=" << stagif.mass << "]: "; verb.verb_prefix = sstream.str(); // Update the verbosity string.
        invif = minv_vector_bicgstab(lhs_cplx[0], rhs_cplx, fine_size, 10000, outer_precision, square_staggered_u1, (void*)&stagif, &verb);
        seq_inversion += invif.ops_count;
    }
    clock_gettime(CLOCK_PROCESS_CPUTIME_ID, &time2); timediff = diff(time1, time2);
    seq_time = ((double)(long int)(1000000000*timediff.tv_sec + timediff.tv_nsec))*1e-9;
    
    
    // Now do multishift. Get ready!
    verb.verb_prefix = "[BICGSTAB-M-CPLX]: ";
    stagif.mass = shifts[0];
    for (n = n_shift - 1; n >= 0; n--)
    {
        shifts[n] -= shifts[0];
    }
    zero<double>(lhs_cplx[0], fine_size);
    
    clock_gettime(CLOCK_PROCESS_CPUTIME_ID, &time1);
    invif = minv_vector_bicgstab_m(lhs_cplx, rhs_cplx, n_shift, fine_size, resid_check_freq, 10000, outer_precision, shifts, square_staggered, (void*)&stagif, true, &verb);
    multi_inversion = invif.ops_count; 
    // Check for failures, fix.
    for (n = 0; n < n_shift; n++)
    {
        if (sqrt(invif.resSqmrhs[n]) > outer_precision*bnorm)
        {
            stagif.mass += shifts[n];
            verb.verb_prefix = "[BICGSTAB-M-CPLX-FIX]: ";
            invif2 = minv_vector_bicgstab(lhs_cplx[n], rhs_cplx, fine_size, 10000, outer_precision, square_staggered_u1, (void*)&stagif, &verb);
            multi_inversion += invif2.ops_count; 
            stagif.mass -= shifts[n];
        }
    }
    clock_gettime(CLOCK_PROCESS_CPUTIME_ID, &time2); timediff = diff(time1, time2);
    shift_time = ((double)(long int)(1000000000*timediff.tv_sec + timediff.tv_nsec))*1e-9;
    
    cout << "[TEST]: BICGSTAB-CPLX Multirhs took " << multi_inversion << " matrix ops and " << shift_time << " seconds, sequential inversions took " << seq_inversion << " matrix ops and " << seq_time << " seconds.\n";
    
    
    for (n = 0; n < n_shift; n++)
    {
        delete[] lhs_cplx[n];
    }
    delete[] lhs_cplx;
    delete[] rhs_cplx; 
    stagif.mass = mass;  
    
    delete[] shifts;
    
    
    
    
    
    // Other tests!
    /*
    // NON-RESTARTED TESTS
    
    // Test an inversion.
    zero<double>(lhs, fine_size);
    clock_gettime(CLOCK_PROCESS_CPUTIME_ID, &time1);
    invif = minv_vector_bicgstab(lhs, rhs, fine_size, 100000, outer_precision, op, (void*)&stagif, &verb);
    clock_gettime(CLOCK_PROCESS_CPUTIME_ID, &time2); timediff = diff(time1, time2); cout << "Time " << ((double)(long int)(1000000000*timediff.tv_sec + timediff.tv_nsec))*1e-9 << "\n";
    
    // Test BiCGStab-l inversions in powers of 2.
    for (i = 1; i <= 16; i*=2)
    {
        sstream.str(string()); // Clear the string stream.
        sstream << "[BICGSTAB-" << i << "]: "; // Set the string stream.
        verb.verb_prefix = sstream.str(); // Update the verbosity string.
        
        zero<double>(lhs, fine_size);
        clock_gettime(CLOCK_PROCESS_CPUTIME_ID, &time1);
        invif = minv_vector_bicgstab_l(lhs, rhs, fine_size, 100000, outer_precision, i, op, (void*)&stagif, &verb);
        clock_gettime(CLOCK_PROCESS_CPUTIME_ID, &time2); timediff = diff(time1, time2); cout << "Time " << ((double)(long int)(1000000000*timediff.tv_sec + timediff.tv_nsec))*1e-9 << "\n";
    }
    
    // Test GCR(\infty)
    verb.verb_prefix = "[GCR]: ";
        
    zero<double>(lhs, fine_size);
    clock_gettime(CLOCK_PROCESS_CPUTIME_ID, &time1);
    invif = minv_vector_gcr(lhs, rhs, fine_size, 100000, outer_precision, op, (void*)&stagif, &verb); // Remark: Dslash count should be doubled. 
    clock_gettime(CLOCK_PROCESS_CPUTIME_ID, &time2); timediff = diff(time1, time2); cout << "Time " << ((double)(long int)(1000000000*timediff.tv_sec + timediff.tv_nsec))*1e-9 << "\n";
    

    // Suggestion from Kate: if the operator is staggered or g5_staggered, test the normal equations with CG.
    if (opt == STAGGERED || opt == G5_STAGGERED)
    {
        verb.verb_prefix = "[CG_NORMAL]: ";
        
        zero<double>(lhs, fine_size);
        
        // Apply D^\dagger
        zero<double>(check, fine_size); gamma_5(check, rhs, (void*)&stagif);
        zero<double>(tmp2, fine_size); square_staggered_u1(tmp2, check, (void*)&stagif);
        zero<double>(check, fine_size); gamma_5(check, tmp2, (void*)&stagif);
        
        clock_gettime(CLOCK_PROCESS_CPUTIME_ID, &time1);
        invif = minv_vector_cg(lhs, check, fine_size, 100000, outer_precision, square_staggered_normal_u1, (void*)&stagif, &verb); // Remark: Dslash count should be doubled. 
        clock_gettime(CLOCK_PROCESS_CPUTIME_ID, &time2); timediff = diff(time1, time2); cout << "Time " << ((double)(long int)(1000000000*timediff.tv_sec + timediff.tv_nsec))*1e-9 << "\n";
    }
    
    // Test restarted BiCGStab(64)
    /*
    int restart_freq = 64; 
    
    // Test an inversion.
    sstream.str(string()); // Clear the string stream.
    sstream << "[BICGSTAB(" << restart_freq << ")]: "; // Set the string stream.
    verb.verb_prefix = sstream.str(); // Update the verbosity string.
    
    zero<double>(lhs, fine_size);
    invif = minv_vector_bicgstab_restart(lhs, rhs, fine_size, 100000, outer_precision, restart_freq, op, (void*)&stagif, &verb);
    
    // Test BiCGStab-l inversions in powers of 2.
    for (i = 1; i <= 16; i*=2)
    {
        sstream.str(string()); // Clear the string stream.
        sstream << "[BICGSTAB-" << i << "(" << restart_freq << ")]: "; // Set the string stream.
        verb.verb_prefix = sstream.str(); // Update the verbosity string.
        
        zero<double>(lhs, fine_size);
        invif = minv_vector_bicgstab_l_restart(lhs, rhs, fine_size, 100000, outer_precision, restart_freq, i, op, (void*)&stagif, &verb);
    }*/
    
    // Free the lattice.
    
    return 0; 
}


// Square lattice.
// Kinetic term for a 2D laplace w/ period bc. Applies lhs = A*rhs.
// The unit vectors are e_1 = xhat, e_2 = yhat.
void square_laplace(double* lhs, double* rhs, void* extra_data)
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
void square_staggered(double* lhs, double* rhs, void* extra_data)
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

// \gamma_5
void gamma_5(double* lhs, double* rhs, void* extra_data)
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

// Square \gamma_5 staggered 2d operator w/out u1 function.
void square_staggered_gamma5(double* lhs, double* rhs, void* extra_data)
{
    staggered_u1_op* stagif = (staggered_u1_op*)extra_data;
    double* tmp = new double[stagif->x_fine*stagif->y_fine];
    
    square_staggered(tmp, rhs, extra_data);
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


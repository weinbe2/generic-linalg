// This test file sees if the staggered operator can be well preconditioned
// by a staggered operator with a laplace term.

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


// For Dslash pieces  
#include "generic_gcr.h"
#include "generic_gcr_var_precond.h"
#include "generic_vector.h"
#include "generic_precond.h"
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
  complex<double> *gfield; // Holds the gauge field.
  complex<double> *lhs, *lhs2, *rhs, *rhs2, *check; // For some Kinetic terms.
  complex<double> *rhs_internal, *lhs_internal, *lhs2_internal;
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
  int square_size = 32; // Can be set on command line with --square_size. 

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
          cout << "--lattice-size [32, 64, 128]           (default 32)\n";
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
  gfield = new complex<double>[2*fine_size];
  lhs = new complex<double>[fine_size];
  rhs = new complex<double>[fine_size];   
  rhs2 = new complex<double>[fine_size];   
  lhs2 = new complex<double>[fine_size];
  check = new complex<double>[fine_size];   
  tmp2 = new complex<double>[fine_size];
  tmp3 = new complex<double>[fine_size];
  lhs_internal = new complex<double>[fine_size];
  lhs2_internal = new complex<double>[fine_size];
  rhs_internal = new complex<double>[fine_size];
  // Zero it out.
  zero<double>(gfield, 2*fine_size);
  zero<double>(rhs, fine_size);
  zero<double>(rhs2, fine_size);
  zero<double>(lhs, fine_size);
  zero<double>(lhs2, fine_size);
  zero<double>(check, fine_size);
  zero<double>(tmp2, fine_size);
  zero<double>(tmp3, fine_size); 
  zero<double>(lhs_internal, fine_size);
  zero<double>(rhs_internal, fine_size); 
  zero<double>(lhs2_internal, fine_size); 
  //

  // Fill stagif.
  stagif.lattice = gfield;
  stagif.mass = mass;  
  stagif.x_fine = x_fine;
  stagif.y_fine = y_fine; 
  stagif.Nc = Nc; // Only relevant for laplace test only.

  // Create the verbosity structure.
  inversion_verbose_struct verb;
  verb.verbosity = VERB_SUMMARY; //VERB_DETAIL;
  verb.precond_verbosity = VERB_NONE; //VERB_SUMMARY; 

  // Describe the gauge field. 
  cout << "[GAUGE]: Creating a gauge field.\n";
  unit_gauge_u1(gfield, x_fine, y_fine);

  // Load the gauge field.
  if (do_load)
  {
      read_gauge_u1(gfield, x_fine, y_fine, load_cfg);
      cout << "[GAUGE]: Loaded a U(1) gauge field from " << load_cfg << "\n";
  }
  else // various predefined cfgs. 
  {
      internal_load_gauge_u1(gfield, x_fine, y_fine, BETA);
  }

  cout << "[GAUGE]: The average plaquette is " << get_plaquette_u1(gfield, x_fine, y_fine) << ".\n";
  cout << "[GAUGE]: The topological charge is " << get_topo_u1(gfield, x_fine, y_fine) << ".\n";


  /**************
  * BEGIN TESTS *
  **************/

  // Set a random source. Should be the same every time b/c we hard code the seed above.
  gaussian<double>(rhs, fine_size, generator); 

  //////
  // First test: invert the staggered Dslash with GCR(1024).
  //////
  
  zero<double>(lhs, fine_size);
  verb.verb_prefix = "[Dslash]: ";
  invif = minv_vector_gcr_restart(lhs, rhs, fine_size, 100000, outer_precision, 1024, square_staggered_u1, (void*)&stagif, &verb); 
  
  //////
  // Second test: invert the staggered Dslash with GCR(1024), preconditioned by GCR(8) of the staggered Dslash.
  //////
  
  // Create a stagif structure for the preconditioned solve.
  staggered_u1_op stagif_precond; 
  stagif_precond.lattice = gfield;
  stagif_precond.mass = mass;  
  stagif_precond.x_fine = x_fine;
  stagif_precond.y_fine = y_fine; 
  stagif_precond.Nc = 1; 
  stagif_precond.wilson_coeff = 0.0;
  
  // Create a GCR(8) preconditioner struct.
  gcr_precond_struct_complex gcrprec;
  gcrprec.n_step = 8;
  gcrprec.rel_res = 1e-10; // guarantee we always do 8 steps.
  gcrprec.matrix_vector = square_staggered_u1;
  gcrprec.matrix_extra_data = &stagif_precond;
  
  // Do the inversion
  zero<double>(lhs, fine_size);
  verb.verb_prefix = "[Dslash_W=0]: ";
  verb.precond_verb_prefix = "[Dslash_W=0_GCR(8)]: ";
  invif = minv_vector_gcr_var_precond_restart(lhs, rhs, fine_size, 100000, outer_precision, 1024, square_staggered_u1, (void*)&stagif, gcr_preconditioner, &gcrprec, &verb); 
  
  //////
  // Third test: invert the staggered Dslash with GCR(1024), preconditioned by GCR(8)
  // of the staggered + laplace Dslash, w/ coefficient 0.01.
  //////
  
  // Set the wilson coeff for this solve.
  stagif_precond.wilson_coeff = 0.01;
  stagif_precond.mass = mass;
  
  // Update GCR precond function with new fcn.
  gcrprec.matrix_vector = square_staggered_2linklaplace_u1;
  
  // Do the inversion
  zero<double>(lhs, fine_size);
  verb.verb_prefix = "[Dslash_W=0p01]: ";
  verb.precond_verb_prefix = "[Dslash_W=0p01_GCR(8)]: ";
  invif = minv_vector_gcr_var_precond_restart(lhs, rhs, fine_size, 100000, outer_precision, 1024, square_staggered_u1, (void*)&stagif, gcr_preconditioner, &gcrprec, &verb); 
  
  //////
  // Fourth test: invert the staggered Dslash with GCR(1024), preconditioned by GCR(8)
  // of the staggered + laplace Dslash, w/ coefficient 0.02.
  //////
  
  // Set the wilson coeff for this solve.
  stagif_precond.wilson_coeff = 0.05;
  stagif_precond.mass = mass;
  
  // Update GCR precond function with new fcn.
  gcrprec.matrix_vector = square_staggered_2linklaplace_u1;
  
  // Do the inversion
  zero<double>(lhs, fine_size);
  verb.verb_prefix = "[Dslash_W=0p05]: ";
  verb.precond_verb_prefix = "[Dslash_W=0p05_GCR(8)]: ";
  invif = minv_vector_gcr_var_precond_restart(lhs, rhs, fine_size, 100000, outer_precision, 1024, square_staggered_u1, (void*)&stagif, gcr_preconditioner, &gcrprec, &verb); 
  
  //////
  // Fifth test: invert the staggered Dslash with GCR(1024), preconditioned by GCR(8)
  // of the staggered + laplace Dslash, w/ coefficient 0.04.
  //////
  
  // Set the wilson coeff for this solve.
  stagif_precond.wilson_coeff = 0.12;
  stagif_precond.mass = mass;
  
  // Update GCR precond function with new fcn.
  gcrprec.matrix_vector = square_staggered_2linklaplace_u1;
  
  // Do the inversion
  zero<double>(lhs, fine_size);
  verb.verb_prefix = "[Dslash_W=0p12]: ";
  verb.precond_verb_prefix = "[Dslash_W=0p12_GCR(8)]: ";
  invif = minv_vector_gcr_var_precond_restart(lhs, rhs, fine_size, 100000, outer_precision, 1024, square_staggered_u1, (void*)&stagif, gcr_preconditioner, &gcrprec, &verb); 
  
  
  /************
  * END TESTS *
  ************/

  // Free the lattice.
  delete[] gfield;
  delete[] lhs;
  delete[] lhs2; 
  delete[] rhs;
  delete[] rhs2;
  delete[] check;
  delete[] tmp2;
  delete[] tmp3; 
  delete[] lhs_internal;
  delete[] rhs_internal; 
  delete[] lhs2_internal; 
  
  

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


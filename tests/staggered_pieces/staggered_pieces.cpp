// This test verifies that the various staggered Dslash functions are working correctly.
// * Compares the \gamma5 D function with applying D then \gamma_5.
// * Compares the D^\dagger function with applying \gamma_5 then D then \gamma_5
// * Compares summing applying D_{eo}, D_{oe}, and the mass term to applying D.
// * Compares summing applying -D_{eo}D_{oe}, -D_{oe}D_{eo}, and the mass term squared to applying D^\dagger D.
// * Compares inverting D directly with even-odd reconstruction.
// * Compares creating and applying a D stencil to the D function.
// * Compares creating and applying a \gamma5 D stencil to the \gamma5 D function.
// * Compares creating and applying a D^\dagger stencil to the D^\dagger function.
// * Compares applying a \gamma5 D stencil to applying a D stencil then a (-1)^coord function.
// * Compares applying a D^\dagger stencil to applying a (-1)^coord function then D stencil then (-1)^coord function.
// * Compares summing applying D_{eo}, D_{oe}, and the mass term via stencil to applying a D stencil.
// * Compares summing applying -D_{eo}D_{oe}, -D_{oe}D_{eo} via stencil, and the stencil shift term to applying D^\dagger D via stencils.
// * Compares inverting a D stencil directly with even-odd stencil reconstruction.
// * Compares applying a D stencil to applying a D stencil with the hypercube rotated into degrees of freedom.
// * Compares applying a \gamma5 D stencil to performing a dof rotation, applying D_internal, applying \sigma3, undoing dof rotation.
// * Compares applying a D^\dagger stencil to performing a dof rotation, applying \sigma3, D_internal, \sigma3, undoing dof rotation.
// * Compares applying a D stencil to performing a dof rotation, summing applying D_{tb}, D_{bt}, and the mass, then undoing dof rotation.
// * Compares summing applying an internal -D_{tb}D_{bt}, -D_{bt}D_{tb}, and the stencil shift term to applying D^\dagger D via stencils.
// * Compares inverting a D_internal stencil directly with top-bottom stencil reconstruction.




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


// For Dslash pieces and even/odd preconditioning tests. 
#include "generic_bicgstab_l.h"
#include "generic_cg.h"
#include "generic_vector.h"
#include "verbosity.h"
#include "u1_utils.h"
#include "operators.h"

// For stencil related tests.
#include "lattice.h"
#include "lattice_functions.h"
#include "operators_stencil.h"

// For internal dof related tests. 
#include "mg.h"
#include "mg_complex.h"
#include "null_gen.h"

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
  verb.verb_prefix = "[TEST5-D]: ";
  invif = minv_vector_bicgstab_l(lhs, rhs, fine_size, 100000, outer_precision, 4, square_staggered_u1, (void*)&stagif, &verb);

  // Invert D using an even/odd decomposition, solving the even system with CG, reconstructing odd.

  verb.verb_prefix = "[TEST5-D_EO_PRE]: ";
  
  // Prepare rhs: m rhs_e - D_{eo} rhs_o
  square_staggered_eoprec_prepare(rhs2, rhs, (void*)&stagif);
  
  // Perform even/odd inversion
  invif = minv_vector_cg(tmp2, rhs2, fine_size, 1000000, outer_precision, square_staggered_m2mdeodoe_u1, (void*)&stagif, &verb);

  // Reconstruct odd: m^{-1}*(rhs_o - D_{oe} lhs2_e)
  square_staggered_eoprec_reconstruct(lhs2, tmp2, rhs, (void*)&stagif);

  cout << "[TEST5]: Even/odd precond test:  relative residual " << diffnorm2sq(lhs, lhs2, fine_size) << ".\n";

  ////////////////////
  //* Compares creating and applying a D stencil to the D function.
  ////////////////////

  // Apply D
  square_staggered_u1(lhs, rhs, (void*)&stagif);

  // Prepare a lattice to hold the staggered stencil.
  const int nd = 2; // 2 dimensions
  int lattice_size[2]; lattice_size[0] = x_fine; lattice_size[1] = y_fine;
  Lattice latt_staggered(nd, lattice_size, 1); // 1 for one internal degree of freedom.

  // Prepare a stencil.
  stencil_2d stencil_staggered(&latt_staggered, get_stencil_size(STAGGERED)); // stencil size is 1.

  // Get the staggered stencil.
  get_square_staggered_u1_stencil(&stencil_staggered, &stagif);

  // Apply D via the staggered stencil.
  apply_stencil_2d(lhs2, rhs, &stencil_staggered);

  cout << "[TEST6]: D stencil test:  relative residual " << diffnorm2sq(lhs, lhs2, fine_size) << ".\n";

  ////////////////////
  //* Compares creating and applying a \gamma5 D stencil to the \gamma5 D function.
  ////////////////////

  // Apply \gamma5 D
  square_staggered_gamma5_u1(lhs, rhs, (void*)&stagif);

  // Prepare a lattice to hold the \gamma5 staggered stencil.
  Lattice latt_staggered_gamma5(nd, lattice_size, 1); // 1 for one internal degree of freedom.

  // Prepare a stencil.
  stencil_2d stencil_staggered_gamma5(&latt_staggered_gamma5, get_stencil_size(STAGGERED)); // stencil size is 1.

  // Get the staggered stencil.
  get_square_staggered_gamma5_u1_stencil(&stencil_staggered_gamma5, &stagif);

  // Apply \gamma5 D via the staggered stencil.
  apply_stencil_2d(lhs2, rhs, &stencil_staggered_gamma5);

  cout << "[TEST7]: \\gamma5 D stencil test:  relative residual " << diffnorm2sq(lhs, lhs2, fine_size) << ".\n";

  ////////////////////
  //* Compares creating and applying a D^\dagger stencil to the D^\dagger function.
  ////////////////////

  // Apply D^\dagger
  square_staggered_dagger_u1(lhs, rhs, (void*)&stagif);

  // Prepare a lattice to hold the staggered dagger stencil.
  Lattice latt_staggered_dagger(nd, lattice_size, 1); // 1 for one internal degree of freedom.

  // Prepare a stencil.
  stencil_2d stencil_staggered_dagger(&latt_staggered_dagger, get_stencil_size(STAGGERED)); // stencil size is 1.

  // Get the staggered stencil.
  get_square_staggered_dagger_u1_stencil(&stencil_staggered_dagger, &stagif);

  // Apply D^\dagger via the staggered stencil.
  apply_stencil_2d(lhs2, rhs, &stencil_staggered_dagger);

  cout << "[TEST8]: D^\\dagger stencil test:  relative residual " << diffnorm2sq(lhs, lhs2, fine_size) << ".\n";

  ////////////////////
  // * Compares applying a \gamma5 D stencil to applying a D stencil then a (-1)^coord function.
  ////////////////////
  
  // Apply \gamma5 D stencil.
  apply_stencil_2d(lhs, rhs, &stencil_staggered_gamma5);
  
  // Apply a D stencil.
  apply_stencil_2d(tmp2, rhs, &stencil_staggered);
  
  // Apply a lattice \epsilon function.
  lattice_epsilon(lhs2, tmp2, stencil_staggered.lat);
  
  cout << "[TEST9]: \\gamma5 D stencil vs epsilon D stencil test:  relative residual " << diffnorm2sq(lhs, lhs2, fine_size) << ".\n";
  
  ////////////////////
  // * Compares applying a D^\dagger stencil to applying a (-1)^coord function then D stencil then (-1)^coord function.
  ////////////////////
  
  // Apply D^\dagger stencil.
  apply_stencil_2d(lhs, rhs, &stencil_staggered_dagger);
  
  // Apply a lattice \epsilon function.
  lattice_epsilon(lhs2, rhs, stencil_staggered.lat);
  
  // Apply a D stencil.
  apply_stencil_2d(tmp2, lhs2, &stencil_staggered);
  
  // Apply a lattice \epsilon function.
  lattice_epsilon(lhs2, tmp2, stencil_staggered.lat);
  
  cout << "[TEST10]: D^\\dagger stencil vs epsilon D epsilon stencil test:  relative residual " << diffnorm2sq(lhs, lhs2, fine_size) << ".\n";
  
  ////////////////////
  // * Compares summing applying D_{eo}, D_{oe}, and the mass term via stencil to applying a D stencil.
  ////////////////////

  // Apply D via the staggered stencil.
  apply_stencil_2d(lhs, rhs, &stencil_staggered);

  // Apply D_{eo}, D_{oe}, m in pieces via stencil.
  apply_stencil_2d_eo(lhs2, rhs, &stencil_staggered); // D_{eo} piece in lhs2.
  apply_stencil_2d_oe(tmp2, rhs, &stencil_staggered); // D_{oe} piece in tmp2.
  for (i = 0; i < stencil_staggered.lat->get_lattice_size(); i++)
  {
      lhs2[i] = lhs2[i] + tmp2[i] + stencil_staggered.shift*rhs[i];
  }

  cout << "[TEST11]: D_{eo}+D_{oe}+m stencil test:  relative residual " << diffnorm2sq(lhs, lhs2, fine_size) << ".\n";

  ////////////////////
  // * Compares summing applying -D_{eo}D_{oe}, -D_{oe}D_{eo} via stencil, and the stencil shift term to applying D^\dagger D via stencils.
  ////////////////////
  
  // Apply D^\dagger D via stencils. (To do: get D^\dagger D stencil directly.)
  apply_stencil_2d(tmp2, rhs, &stencil_staggered);
  apply_stencil_2d(lhs, tmp2, &stencil_staggered_dagger);

  // Apply D_{eo}D_{oe} via stencils.
  apply_stencil_2d_eo(tmp2, rhs, &stencil_staggered);
  apply_stencil_2d_oe(lhs2, tmp2, &stencil_staggered);

  // Apply D_{oe}D_{eo} via stencils.
  apply_stencil_2d_oe(tmp2, rhs, &stencil_staggered);
  apply_stencil_2d_eo(tmp3, tmp2, &stencil_staggered);
  
  // Combine into m^2 - D_{eo}D_{oe} - D_{oe}D_{eo}.
  for (i = 0; i < stencil_staggered.lat->get_lattice_size(); i++)
  {
      lhs2[i] = stencil_staggered.shift*stencil_staggered.shift*rhs[i] - lhs2[i] - tmp3[i];
  }

  cout << "[TEST12]: D^\\dagger D stencil test:  relative residual " << diffnorm2sq(lhs, lhs2, fine_size) << ".\n";
  
  ////////////////////
  // * Compares inverting a D stencil directly with even-odd stencil reconstruction.
  ////////////////////
  
  // Invert D directly using BiCGstab-L.
  verb.verb_prefix = "[TEST13-D]: ";
  invif = minv_vector_bicgstab_l(lhs, rhs, fine_size, 100000, outer_precision, 4, apply_stencil_2d, &stencil_staggered, &verb);

  // Invert D using an even/odd decomposition, solving the even system with CG, reconstructing odd.

  verb.verb_prefix = "[TEST13-D_EO_PRE]: ";
  
  // Prepare rhs: m rhs_e - D_{eo} rhs_o
  apply_square_staggered_eoprec_prepare_stencil(rhs2, rhs, &stencil_staggered);
  
  // Perform even/odd inversion
  invif = minv_vector_cg(tmp2, rhs2, fine_size, 1000000, outer_precision, apply_square_staggered_m2mdeodoe_stencil, &stencil_staggered, &verb);

  // Reconstruct odd: m^{-1}*(rhs_o - D_{oe} lhs2_e)
  apply_square_staggered_eoprec_reconstruct_stencil(lhs2, tmp2, rhs, &stencil_staggered);

  cout << "[TEST13]: Even/odd precond stencil test:  relative residual " << diffnorm2sq(lhs, lhs2, fine_size) << ".\n";

  ////////////////////
  // * Compares applying a D stencil to applying a D stencil with the hypercube rotated into degrees of freedom.
  ////////////////////
  
  // Apply D via the staggered stencil.
  apply_stencil_2d(lhs, rhs, &stencil_staggered);
  
  // We're abusing the multigrid/null vector machinary we have to do this test.
  
  // Prepare a null vector structure. The null vector structure contains a lot of fields, many of which
  // are related to generating null vectors via relaxation. Since we only want to do a unitary transformation,
  // just rotating the hypercube into internal degrees of freedom, we only need to set a few pieces.
  null_vector_params nvec_params;
  nvec_params.opt_null = STAGGERED; // using a STAGGERED op.
  nvec_params.n_null_vector = 1; // only "one" null vector...
  nvec_params.null_partitions = 4; // partitioned into 4...
  nvec_params.bstrat = BLOCK_CORNER; // corners of the hypercube.
  
  // Prepare a lattice for the internal d.o.f. staggered operator.
  lattice_size[0] = x_fine/2; lattice_size[1] = y_fine/2; // We're rotating the 2^2 hypercube in.
  Lattice latt_staggered_internal(nd, lattice_size, 4); // 4 for the 4 internal degrees of freedom.
  
  // Prepare a multigrid operator struct. We're not doing anything related to multigrid explicitly,
  // but we're going to take advantage of the prolong and restrict functions to do the necessary
  // unitary rotations to go between hypercubes and internal degrees of freedom.
  // There's a lot of cleanup that could be done within the mg struct. For ex, there are a lot of state-related
  // variables that could just be queried on the fly from the current lattice.
  mg_operator_struct_complex mgstruct;
  
  mgstruct.x_fine = x_fine; // top level x dimension.
  mgstruct.y_fine = y_fine; // top level y dimension.
  mgstruct.n_refine = 1; // we have two "levels", so there's one "refinement": top level is the original lattice,
                             // bottom is the internal lattice.
  
  mgstruct.blocksize_x = new int[1]; mgstruct.blocksize_x[0] = 2; // block by 2.
  mgstruct.blocksize_y = new int[1]; mgstruct.blocksize_y[0] = 2;
  mgstruct.Nc = 1; // Nc = 1 on the top level, i.e., there are no internal degrees of freedom.
  
  mgstruct.latt = new Lattice*[2]; // array of pointers to the lattices.
  mgstruct.latt[0] = &latt_staggered;
  mgstruct.latt[1] = &latt_staggered_internal;
  
  mgstruct.stencils = new stencil_2d*[2]; // array of pointers to stencils.
  mgstruct.stencils[0] = &stencil_staggered; // the top level is the fine stencil.
  // We'll put the internal dof stencil in momentarily. 
  
  mgstruct.have_dagger_stencil = false; // we don't need the daggered stencil. (This is used for normal smoother solves.)
  mgstruct.dagger_stencils = 0;
  
  mgstruct.n_vector = 4; // We have 4 null vectors. This is currently hard coded to be constant for every level.
  mgstruct.null_vectors = new complex<double>**[1]; // This holds the null vectors. We have only one refinement.
  mgstruct.null_vectors[0] = new complex<double>*[mgstruct.n_vector]; // We'll have 4 "null vectors", one for each corner.
  for (i = 0; i < mgstruct.n_vector; i++)
  {
    mgstruct.null_vectors[0][i] = new complex<double>[latt_staggered.get_lattice_size()];
    zero<double>(mgstruct.null_vectors[0][i], mgstruct.latt[0]->get_lattice_size()); 
  }
  
  mgstruct.matrix_vector = square_staggered_u1; // What's the original fine function?
  mgstruct.matrix_vector_dagger = 0; // We don't need the daggered function.
  mgstruct.matrix_extra_data = &stagif; // This is the structure that goes along with the fine function.
  
  mgstruct.curr_level = 0; // state variable, the fine level is currently the top (0th) level.
  
  mgstruct.curr_dof_fine = mgstruct.latt[0]->get_nc(); // the fine level dof.
  mgstruct.curr_x_fine = mgstruct.latt[0]->get_lattice_dimension(0); // fine level x dimension
  mgstruct.curr_y_fine = mgstruct.latt[0]->get_lattice_dimension(1); // fine level y dimension.
  mgstruct.curr_fine_size = mgstruct.latt[0]->get_lattice_size(); // total size.
  
  mgstruct.curr_dof_coarse = mgstruct.latt[1]->get_nc(); // the coarse level dof.
  mgstruct.curr_x_coarse = mgstruct.latt[1]->get_lattice_dimension(0); // coarse level x dimension
  mgstruct.curr_y_coarse = mgstruct.latt[1]->get_lattice_dimension(1); // coarse level y dimension.
  mgstruct.curr_coarse_size = mgstruct.latt[1]->get_lattice_size(); // total size.
  
  // Whew, that was a lot. Now that we've got that all filled out, we have a function that'll automatically generate
  // the unitary transformation null vectors (which are the free field null vectors when you do a 2^2 block size).
  null_generate_free(&mgstruct, &nvec_params); 
  
  // We need to properly block normalize the vectors. This is what makes the transformation unitary.
  // By construction, the vectors are already orthogonal, so the "ortho" part doesn't do anything.
  block_orthonormalize(&mgstruct); 
  
  // Create a stencil for the internal dof operator.
  stencil_2d stencil_staggered_internal(&latt_staggered_internal, get_stencil_size(STAGGERED)); // stencil size is 1.
  mgstruct.stencils[1] = &stencil_staggered_internal;
  
  // Using the fine stencil and the transformation defined by the free null vectors, we can build the
  // internal dof stencil.
  generate_coarse_from_fine_stencil(&stencil_staggered_internal, &stencil_staggered, &mgstruct, true); // true -> don't build the mass explicitly into the stencil.
  stencil_staggered_internal.shift = stencil_staggered.shift; // the mass is unchanged by the unitary transformation.
  
  // Alright, everything's done! We're now ready to actually do a test.
  
  // Perform a unitary transformation to take the original staggered layout to the internal dof layout.
  // This function is why we went through all of the misery to set up a mg struct.
  restrict(rhs_internal, rhs, &mgstruct);
  
  // Apply the internal dof stencil.
  apply_stencil_2d(lhs_internal, rhs_internal, &stencil_staggered_internal);
  
  // Perform a unitary transformation to return to the original staggered layout.
  prolong(lhs2, lhs_internal, &mgstruct);
  
  cout << "[TEST14]: Internal dof stencil test:  relative residual " << diffnorm2sq(lhs, lhs2, fine_size) << ".\n";

  ////////////////////
  // * Compares applying a \gamma5 D stencil to performing a dof rotation, applying D_internal, applying \sigma3, undoing dof rotation.
  ////////////////////
  
  // Apply \gamma5 D via the staggered stencil.
  apply_stencil_2d(lhs, rhs, &stencil_staggered_gamma5);
  
  // Apply an internal unitary transformation.
  restrict(rhs_internal, rhs, &mgstruct);
  
  // Apply D_internal
  apply_stencil_2d(lhs_internal, rhs_internal, &stencil_staggered_internal);
  
  // Apply \sigma3
  lattice_sigma3(tmp2, lhs_internal, &latt_staggered_internal);
  
  // Apply a hypercube unitary transformation.
  prolong(lhs2, tmp2, &mgstruct);
  
  cout << "[TEST15]: \\gamma5 D stencil vs internal \\sigma3 D stencil test:  relative residual " << diffnorm2sq(lhs, lhs2, fine_size) << ".\n";
  
  ////////////////////
  // * Compares applying a D^\dagger stencil to performing a dof rotation, applying \sigma3, D_internal, \sigma3, undoing dof rotation.
  ////////////////////
  
  // Apply D^\dagger via the staggered stencil.
  apply_stencil_2d(lhs, rhs, &stencil_staggered_dagger);
  
  // Apply an internal unitary transformation.
  restrict(rhs_internal, rhs, &mgstruct);
  
  // Apply \sigma3
  lattice_sigma3(lhs_internal, rhs_internal, &latt_staggered_internal);
  
  // Apply D_internal
  apply_stencil_2d(tmp2, lhs_internal, &stencil_staggered_internal);
  
  // Apply \sigma3
  lattice_sigma3(lhs_internal, tmp2, &latt_staggered_internal);
  
  // Apply a hypercube unitary transformation.
  prolong(lhs2, lhs_internal, &mgstruct);
  
  cout << "[TEST16]: D^\\dagger stencil vs internal \\sigma3 D \\sigma3 stencil test:  relative residual " << diffnorm2sq(lhs, lhs2, fine_size) << ".\n";
  
  ////////////////////
  // * Compares applying a D stencil to performing a dof rotation, summing applying D_{tb}, D_{bt}, and the mass, then undoing dof rotation.
  ////////////////////
  
  // Apply D via the staggered stencil.
  apply_stencil_2d(lhs, rhs, &stencil_staggered);
  
  // Apply an internal unitary transformation.
  restrict(rhs_internal, rhs, &mgstruct);
  
  // Apply internal D_{tb}, D_{bt}, m in pieces.
  apply_stencil_2d_tb(lhs_internal, rhs_internal, &stencil_staggered_internal); // D_{tb} piece in lhs_internal.
  apply_stencil_2d_bt(tmp2, rhs_internal, &stencil_staggered_internal); // D_{bt} piece in tmp2.
  for (i = 0; i < latt_staggered_internal.get_lattice_size(); i++)
  {
      lhs_internal[i] = lhs_internal[i] + tmp2[i] + stencil_staggered_internal.shift*rhs_internal[i]; // combine, add mass.
  }
  
  // Apply a hypercube unitary transformation
  prolong(lhs2, lhs_internal, &mgstruct);
  
  cout << "[TEST17]: D stencil vs internal D_{tb}+D_{bt}+mass test:  relative residual " << diffnorm2sq(lhs, lhs2, fine_size) << ".\n";
  
  ////////////////////
  // * Compares summing applying an internal -D_{tb}D_{bt}, -D_{bt}D_{tb}, and the stencil shift term to applying D^\dagger D via stencils.
  ////////////////////
  
  // Apply D^\dagger D via stencils. (To do: get D^\dagger D stencil directly.)
  apply_stencil_2d(tmp2, rhs, &stencil_staggered);
  apply_stencil_2d(lhs, tmp2, &stencil_staggered_dagger);
  
  // Apply an internal unitary transformation.
  restrict(rhs_internal, rhs, &mgstruct);

  // Apply D_{tb}D_{bt} via stencils.
  apply_stencil_2d_bt(tmp2, rhs_internal, &stencil_staggered_internal);
  apply_stencil_2d_tb(lhs_internal, tmp2, &stencil_staggered_internal);

  // Apply D_{bt}D_{tb} via stencils.
  apply_stencil_2d_tb(tmp2, rhs_internal, &stencil_staggered_internal);
  apply_stencil_2d_bt(tmp3, tmp2, &stencil_staggered_internal);
  
  // Combine into m^2 - D_{tb}D_{bt} - D_{bt}D_{tb}.
  for (i = 0; i < stencil_staggered_internal.lat->get_lattice_size(); i++)
  {
      lhs_internal[i] = stencil_staggered_internal.shift*stencil_staggered_internal.shift*rhs_internal[i] - lhs_internal[i] - tmp3[i];
  }
  
  // Apply a hypercube unitary transformation
  prolong(lhs2, lhs_internal, &mgstruct);

  cout << "[TEST18]: D^\\dagger D internal stencil test:  relative residual " << diffnorm2sq(lhs, lhs2, fine_size) << ".\n";
  
  ////////////////////
  // * Compares inverting a D_internal stencil directly with top-bottom stencil reconstruction.
  ////////////////////
  
  // Apply an internal unitary transformation.
  restrict(rhs_internal, rhs, &mgstruct);
  
  // Invert D_internal directly using BiCGstab-L.
  verb.verb_prefix = "[TEST19-D_INTERNAL]: ";
  invif = minv_vector_bicgstab_l(lhs_internal, rhs_internal, fine_size, 100000, outer_precision, 4, apply_stencil_2d, &stencil_staggered_internal, &verb);

  // Invert D using an even/odd decomposition, solving the even system with CG, reconstructing odd.

  verb.verb_prefix = "[TEST19-D_INTERNAL_TB_PRE]: ";
  
  // Prepare rhs: m rhs_t - D_{tb} rhs_b
  apply_square_staggered_tbprec_prepare_stencil(rhs2, rhs_internal, &stencil_staggered_internal);
  
  // Perform even/odd inversion
  invif = minv_vector_cg(tmp2, rhs2, fine_size, 1000000, outer_precision, apply_square_staggered_m2mdtbdbt_stencil, &stencil_staggered_internal, &verb);

  // Reconstruct odd: m^{-1}*(rhs_o - D_{oe} lhs2_e)
  apply_square_staggered_tbprec_reconstruct_stencil(lhs2_internal, tmp2, rhs_internal, &stencil_staggered_internal);
  
  // Apply a hypercube unitary transformation
  prolong(lhs, lhs_internal, &mgstruct);
  prolong(lhs2, lhs2_internal, &mgstruct);

  cout << "[TEST19]: Internal top/bottom precond stencil test:  relative residual " << diffnorm2sq(lhs, lhs2, fine_size) << ".\n";
  
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
  
  // Free null vector and multigrid related things.
  delete[] mgstruct.blocksize_x;
  delete[] mgstruct.blocksize_y;
  delete[] mgstruct.latt;
  delete[] mgstruct.stencils; 
  for (i = 0; i < mgstruct.n_vector; i++)
  {
    delete[] mgstruct.null_vectors[0][i];
  }
  delete[] mgstruct.null_vectors[0];
  delete[] mgstruct.null_vectors;
  

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


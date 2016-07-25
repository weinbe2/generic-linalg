// Fri Jan 15 13:06:07 EST 2016
// Evan S Weinberg weinbe2@bu.edu
// C++ file for multigrid on square laplace using GMRES(n) for internal iteration.
// Also supports using SOR as a smoother. 
// Based heavily on "mg.c" by Rich Brower!

#include <cstdio>
#include <cmath>
#include <iostream>
#include <complex>
#include <sstream>

using namespace std;

// Where GMRES is defined!
#include "generic_inverters.h"
#include "generic_vector.h"

// Structure passed to generic square laplace function.
typedef struct{
  int N;  // Linear dimension of top level. Fixed at start.
  int Lmax; // Maximum number of levels. Fixed at start.
  int size[20]; // Size of each level. Fixed at start.
  double a[20]; // Lattice spacing at each level. Fixed at start.
  double m; // Mass entering laplace equation. Fixed at start.
  //double scale[20]; 
  int level; // what level are we iterating on? Changes. 
} param_t;

// Generic square laplace. Expects the void* to be of type param_t*. 
void square_laplacian(double* lhs, double* rhs, void* extra_data);

// Project the residual down coarse <- fine.
void project_residual(double* res_c, double* res_f, param_t p);

// Interpolate the coarse solution up and add it to the fine solution.
void interpolate_add(double* soln_f, double* soln_c, param_t p);

// Can optionally pass sizes in.
// Arg 1: size (take log2, round to nearest power of 2.)
// Arg 2: mass
// Arg 3: levels
// Arg 4: GMRES(n)
int main(int argc, char** argv)
{
  // Just a simple test! 
  
  // Declare some variables.
  int i, j, x, y;
  int nlev, lev;
  double *phi[20], *res[20], *check, *soln, *b; // For some Kinetic terms.
  double explicit_resid = 0.0;
  double bnorm, resid;
  int iter, rest; 
  int dslash[20];
  param_t p; 
  inversion_info invif;
  
  // For smoothing
  int smoother, smoother_steps;
  double omega; 
  
  // Set the lattice size, and other such things.
  if (argc == 5 || argc == 6)
  {
    stringstream ss;
    
    double size_tmp;
    ss << argv[1]; // size.
    ss >> size_tmp;
    p.Lmax = ((int)(log(size_tmp)/log(2)+0.5))-1;
    p.N = 2*(int)(pow(2,p.Lmax)+0.5);
    ss.clear();
    
    ss << argv[2]; // mass
    ss >> p.m;
    ss.clear();
    
    ss << argv[3]; // levels
    ss >> nlev;
    ss.clear();
    
    ss << argv[4]; // GMRES(n)
    ss >> rest;
    ss.clear();
    
    if (argc == 6)
    {
      ss << argv[5]; // number of smoother steps.
      ss >> smoother_steps;
      ss.clear();
    }
    else
    {
      smoother_steps = 0; // don't smooth.
    }
  }
  else
  {
    cout << "This program solves the 2D laplace equation \\ del^2+m^2 with periodic boundary conditions.\n";
    cout << "Call as: ./mg_square_laplace.cpp [length] [mass] [levels] [gmres(n)] {n_smooth}\n";
    cout << "   [length]:    The length in one dimension. We solve on an L by L space.\n";
    cout << "   [mass]:      The IR regulator mass.\n";
    cout << "   [levels]:    The number of levels to do down. The maximum is log_2(L)-1.\n";
    cout << "   [gmres(n)]:  The n for our GMRES(n) inner solver.\n";
    cout << "   {n_smooth}:  The number of smoothing steps to use. (Default: 0)\n";
    return 0;
  }
  
  if (smoother_steps <= 0) // if number of smoother steps is 0 or negative
  {
    smoother = 0; // don't smooth!
  }
  else // smooth
  {
    smoother = 1;
    
    // What is the relaxation parameter?
    // Needs to be set s.t. (1-omega*|largest eigenvalue|) < 1. 
    // For 2D laplace, |largest eigenvalue| is 4+m^2 -> use omega = (2/3)/(4+m^2)
  omega = 2.0/(3.0*(4.0+p.m*p.m));
  }
    
  
  
  
  
  printf("Start Length %d Mass %15.20e Levels %d GMRES %d\n", p.N, p.m, nlev, rest);
  if (smoother == 1)
  {
    printf("Using a %d step smoother with omega = %.8e.\n", smoother_steps, omega);
  }
  
  if(nlev  > p.Lmax)
  { 
    cout << "ERROR More levels than available in lattice!" << endl;
    return 0;
  }
  
  cout << endl << "V cycle for " << p.N << " by " << p.N << " lattice with nlev = ";
  cout << nlev << " out of max " << p.Lmax << "." << endl;
  
  // Initialize the size and lattice spacing for each level.
  p.size[0] = p.N;
  p.a[0] = 1.0;
  for (i = 1; i < p.Lmax+1; i++)
  {
    p.size[i] = p.size[i-1]/2;
    p.a[i] = 2.0*p.a[i-1];
  }
  
  // Initialize the lattice for one level. Indexing: index = y*N + x.
  for (i=0;i<=p.Lmax;i++)
  {
    phi[i] = new double[p.size[i]*p.size[i]];
    res[i] = new double[p.size[i]*p.size[i]];
    for (j=0;j<p.size[i]*p.size[i];j++)
    {
      phi[i][j] = 0.0;
      res[i][j] = 0.0;
    }
    dslash[i] = 0; // count of dslash operations. 
  }
  check = new double[p.N*p.N];   
  soln = new double[p.N*p.N];
  b = new double[p.N*p.N];
  
  // Zero it out.
  for (i=0; i < p.N*p.N; i++)
  {
    b[i] = 0.0;
    check[i] = 0.0;
    soln[i] = 0.0;
  }

  // Set a point on the rhs.
  b[p.N/2+(p.N/2)*p.N] = 1.0;

  // Get norm for rhs. We use this to get the relative residual.
  bnorm = 0.0;
  for (i=0;i<p.N*p.N;i++)
  {
    bnorm = bnorm + b[i]*b[i];
  }
  bnorm = sqrt(bnorm);

  // Set a point on the lhs.
  //phi[0][p.N/2+(p.N/2)*p.N+1] = 1.0;

  cout << "Solving A [lhs] = [rhs] for lhs, using a point source." << endl;

  // lhs = A^(-1) rhs
  // Arguments:
  // 1: lhs
  // 2: rhs
  // 3: size of vector
  // 4: maximum iterations
  // 5: residual
  // 5a for gmres_restart: how often to restart.
  // 6: function pointer
  // 7: "extra data": can set this to not-null to pass in gauge fields, etc.

  // Set up a zero level "cycle".
  resid = bnorm;
  for (i=0;i<p.N*p.N;i++)
  {
    res[0][i] = b[i];
  }
  iter = 0;
  while (resid/bnorm > 1e-6)
  {
    
    for (p.level=0;p.level<nlev;p.level++) // go down.
    {
      
      // Perform one iteration of GMRES(8).
      invif = minv_vector_gmres(phi[p.level], res[p.level], p.size[p.level]*p.size[p.level], rest, 1e-10, square_laplacian, &p);
      
      // Perform 8 iterations of BiCGStab
      //invif = minv_vector_bicgstab(phi[p.level], res[p.level], p.size[p.level]*p.size[p.level], rest, 1e-10, square_laplacian, &p);
      
      dslash[p.level] += rest; 
      
      // Project the residual coarse <- fine.
      project_residual(res[p.level+1], res[p.level], p);
      
      // Project residual down.
      /*for (x=0;x<p.size[p.level+1];x++)
      {
        for (y=0;y<p.size[p.level+1];y++)
        {
          // res[lev+1] = P^dag res[lev].
          res[p.level+1][x+y*p.size[p.level+1]] = 0.25*(
            res[p.level][2*x+2*y*p.size[p.level]] +
            res[p.level][2*x+1+2*y*p.size[p.level]] +
            res[p.level][2*x+(2*y+1)*p.size[p.level]] +
            res[p.level][2*x+1+(2*y+1)*p.size[p.level]]);
        }
      }*/

    } 
    
    // We're now at the lowest level. Come back up!
    for (p.level=nlev;p.level>=0;p.level--) // Come up
    {
      if (smoother == 1)
      {
        // Do "smoother_steps" iterations of SOR with parameter omega.
        minv_vector_sor(phi[p.level], res[p.level], p.size[p.level]*p.size[p.level], smoother_steps, 1e-10, omega, square_laplacian, &p);
      }
      
      
      // Perform one iteration of GMRES(8), perhaps using old
      // phi as an initial guess!
      invif = minv_vector_gmres(phi[p.level], res[p.level], p.size[p.level]*p.size[p.level], rest, 1e-10, square_laplacian, &p);
      
      // Perform 8 iterations of BiCGStab
      //invif = minv_vector_bicgstab(phi[p.level], res[p.level], p.size[p.level]*p.size[p.level], rest, 1e-10, square_laplacian, &p);
      
      dslash[p.level] += rest;
      
      // We've now done a solve at the lowest level!
      if (p.level > 0)
      {
        // Interpolate phi up to the next level!
        interpolate_add(phi[p.level-1], phi[p.level], p);
        
        // Clear out coarse phi for the next iteration.
        for(x=0;x<p.size[p.level];x++)
        {
          for(y=0;y<p.size[p.level];y++)
          {
            phi[p.level][x+y*p.size[p.level]] = 0.0;
          }
        }
      } 
    }
    
    p.level = 0;
    // Add phi[0] to soln. This gives us a new x_0.
    for (i=0;i<p.size[p.level]*p.size[p.level];i++)
    {
      soln[i] = soln[i] + phi[0][i];
    }
  
    // Compute new residual. Can technically pull this from top level minv_* call,
    // which is why it doesn't count toward dslash count. 
    square_laplacian(check, phi[0], &p);
    for (i=0;i<p.size[0]*p.size[0];i++)
    {
      res[0][i] = res[0][i] - check[i];
      phi[lev][i] = 0.0;
    }
    
    resid = sqrt(norm2sq<double>(res[0],p.size[p.level]*p.size[p.level])); 
    
    cout << "Iter " << ++iter << " Relative residual " << resid/bnorm << "." << endl;
  
  }

  //invif = minv_vector_cg(soln, b, p.N*p.N, 4000, 1e-8, square_laplacian, &p);
  //invif = minv_vector_gcr(lhs, rhs, p.N*p.N, 4000, 1e-8, square_laplacian, NULL);
  //invif = minv_vector_bicgstab(lhs, rhs, p.N*p.N, 4000, 1e-8, square_laplacian, NULL);
  //invif = minv_vector_gmres_norestart(phi[0], res[0], p.N*p.N, 10000, 1e-6, square_laplacian, &p);
  //invif = minv_vector_gmres_restart(soln, b, p.N*p.N, 50000, 1e-6, 8, square_laplacian, &p);
  //invif = minv_vector_gmres_restart(lhs, rhs, p.N*p.N, 4000, 1e-8, 20, square_laplacian, NULL);


  printf("Algorithm %s took %d iterations to reach a residual of %.8e.\n", invif.name.c_str(), invif.iter, sqrt(invif.resSq));

  printf("Computing [check] = A [lhs] as a confirmation.\n");

  // Check and make sure we get the right answer.
  square_laplacian(check, soln, &p);

  for (i = 0; i < p.N*p.N; i++)
  {
    explicit_resid += (b[i] - check[i])*(b[i] - check[i]);
  }
  explicit_resid = sqrt(explicit_resid);

  printf("[check] should equal [rhs]. The residual is %15.20e.\n", explicit_resid/bnorm);

  printf("-------------------------\n");
  
  // Find number of ops!
  // Specific to sq laplace, Each site on each iteration requires:
  // Floating point: 1 multiply, 5 adds, on top of other things!
  // Integer:        8 adds, 4 multiply, 4 mod.
  long fmul = 0;
  long fadd = 0;
  long imul = 0;
  long iadd = 0;
  long imod = 0;
  for (i=0;i<=p.Lmax;i++)
  {
    //      # ops             N*N         indiv arith
    fmul += dslash[i]*p.size[i]*p.size[i]*1;
    fadd += dslash[i]*p.size[i]*p.size[i]*5;
    imul += dslash[i]*p.size[i]*p.size[i]*4;
    iadd += dslash[i]*p.size[i]*p.size[i]*8;
    imod += dslash[i]*p.size[i]*p.size[i]*4;
  }
  long ops = fmul+fadd+imul+iadd+imod;
  
  printf("Fin Len %d Mass %1.5g Lvls %d GMRES %d Itr %d Ops %ld Fmul %ld Fadd %ld Imul %ld Iadd %ld Imod %ld Res %1.5e\n", p.N, p.m, nlev, rest, iter, ops, fmul, fadd, imul, iadd, imod, explicit_resid/bnorm);
  //printf("Fin \n", iter, explicit_resid/bnorm);
  for (i=0;i<p.Lmax;i++)
  {
    printf("Level %d Dslash %d\n", i, dslash[i]);
  }

  // Free the lattice.
  for (i=0;i<p.Lmax+1;i++)
  {
    delete[] phi[i];
    delete[] res[i];
  }
  delete[] check;
  delete[] soln;
  delete[] b;


  return 0;
}

// Project the residual coarse <- fine.
void project_residual(double* res_c, double* res_f, param_t p)
{
  int x, y;
  int L_fine, L_coarse;
  
  L_fine = p.size[p.level];
  L_coarse = p.size[p.level+1];
  
  for (x=0;x<L_coarse;x++)
  {
    for (y=0;y<L_coarse;y++)
    {
      // res[lev+1] = P^dag res[lev].
      res_c[x+y*L_coarse] = 0.25*(
        res_f[2*x+2*y*L_fine] +
        res_f[2*x+1+2*y*L_fine] +
        res_f[2*x+(2*y+1)*L_fine] +
        res_f[2*x+1+(2*y+1)*L_fine]);
    }
  }
}

// Interpolate the coarse solution up and add it to the fine solution.
void interpolate_add(double* soln_f, double* soln_c, param_t p)
{
  int x, y;
  int L_fine, L_coarse;
  
  L_fine = p.size[p.level-1];
  L_coarse = p.size[p.level];

  for(x=0;x<L_coarse;x++)
  {
    for(y=0;y<L_coarse;y++)
    {
      soln_f[2*x+2*y*L_fine] += soln_c[x+y*L_coarse];
      soln_f[2*x+1+2*y*L_fine] += soln_c[x+y*L_coarse];
      soln_f[2*x+(2*y+1)*L_fine] += soln_c[x+y*L_coarse];
      soln_f[2*x+1+(2*y+1)*L_fine] += soln_c[x+y*L_coarse];
    }
  }

}

// Square lattice.
// Kinetic term for a 2D laplacian w/ period bc. Applies lhs = A*rhs.
// The unit vectors are e_1 = xhat, e_2 = yhat.
// The "extra_data" is really a param_t type.
void square_laplacian(double* lhs, double* rhs, void* extra_data)
{
   // Declare variables.
   int i;
   int x,y;
   int L;
   param_t* p;
   
   // Unbox "extra_data".
   p = (param_t*)extra_data; 
   L = p->size[p->level];
   // For a 2D square lattice, the stencil is:
   //     |  0 -1  0 |
   //     | -1 +4 -1 |
   //     |  0 -1  0 |
   //
   // e2 = yhat
   // ^
   // | 
   // |-> e1 = xhat

   // Each iteration requires:
   // Floating point: 1 multiply, 5 adds, on top of other things!
   // Integer:        8 adds, 4 multiply, 4 mod.
   // Apply the stencil.
   for (i = 0; i < L*L; i++)
   {
      lhs[i] = 0.0;
      x = i%L; // integer mod.
      y = i/L; // integer divide.
      
      // + e1.
      lhs[i] = lhs[i]-rhs[y*L+((x+1)%L)];
     
      // - e1.
      lhs[i] = lhs[i]-rhs[y*L+((x+L-1)%L)]; // The extra +N is because of the % sign convention.
      
      // + e2.
      lhs[i] = lhs[i]-rhs[((y+1)%L)*L+x];
    
      // - e2.
      lhs[i] = lhs[i]-rhs[((y+L-1)%L)*L+x];

      // 0
      // Added mass term here.
      lhs[i] = lhs[i]+(4+p->m*p->m)*rhs[i]; // Could rescale 'm' on each level. Makes sense if each level is a blocking transform.
      //lhs[i] = lhs[i]+(4+p->m*p->m*p->a[p->level]*p->a[p->level])*rhs[i];
   }
       
}


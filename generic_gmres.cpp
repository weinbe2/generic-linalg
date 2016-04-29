// Fri Apr 29 11:54:22 EDT 2016
// Evan S Weinberg
// C++ file for GMRES inverter.
 
// To do:
// 1. Template to support float, double.
// 2. Improve GMRES to get residual from reduced matrix.
//    This saves one matrix op per iteration. (We currently
//    lazily use the matrix op again to get the residual.)

#include <iostream>
#include <cmath>
#include <string>
#include <sstream>
#include <complex>

#include "generic_inverters.h"
#include "generic_vector.h"

using namespace std;


// Solves lhs = A^(-1) rhs using gmres

// At one point during development, we solved the
// subspace normal equations within gmres with CG. 
// These functions and structs support this.
// At this point, we just use Gaussian elimination,
// but these are here for legacy.
/*
struct gmres_struct
{
  double** h;
  int size;
};
void gmres_hTh(double* phi, double* phi0, void* extra_data);
*/

inversion_info minv_vector_gmres_norestart(double  *phi, double  *phi0, int size, int max_iter, double eps, void (*matrix_vector)(double*,double*,void*), void* extra_info)
{
// GMRES solutions to Mphi = b 
//  see https://en.wikipedia.org/wiki/Generalized_minimal_residual_method
// A clearer article: 
// http://www.math.iit.edu/~fass/477577_Chapter_14.pdf
  
  // Solve a counting issue...
  max_iter++;
  
  // Initialize vectors.
  double **q, **h; // Q matrix, Hessenberg matrix.
  double **hTh; // stores normal eqn matrix for hTh. 
  double *res, *tmp, *tmp2;
  double *y, *bhTy; // Holds the reduced solution. 
  double bres, localres; // norm of b, curr residual.
  double beta; // beta is the norm of res = b-Ax_0.
  //gmres_struct gmstr; // Passed to solve lls problem, no longer needed: we use gaussian elim.
  int i, j, k, iter;
  inversion_info invif;
  
  // Prepare q, h. These get larger on each iteration.
  q = new double*[max_iter];
  h = new double*[max_iter-1];
  hTh = new double*[max_iter-1];
  //q = (double**)malloc(max_iter*sizeof(double*));
  //h = (double**)malloc((max_iter-1)*sizeof(double*));
  //hTh = (double**)malloc((max_iter-1)*sizeof(double*));
  //gmstr.h = h; // Depreciated, no longer use CG for normal subspace.
  
  // Prepare space for the residual, and a temporary vector.
  res = new double[size];
  tmp = new double[size];
  tmp2 = new double[size];
  y = new double[max_iter];
  bhTy = new double [max_iter];
  //res = (double*)malloc(size*sizeof(double));
  //tmp = (double*)malloc(size*sizeof(double));
  //tmp2 = (double*)malloc(size*sizeof(double));
  //y = (double*)malloc(max_iter*sizeof(double));
  //bhTy = (double*)malloc(max_iter*sizeof(double));
  
  // Compute the residual. Ax' = res gives the solution phi+x'.
  (*matrix_vector)(tmp, phi, extra_info);
  for (i=0;i<size;i++) // res = b - Ax_0
  {
    res[i] = phi0[i]-tmp[i]; 
  }
  
  // Initialize the first q vector. q is res/||res||
  q[0] = new double[size];
  //q[0] = (double*)malloc(size*sizeof(double));
  
  beta = 0.0; bres = 0.0;
  for (i=0;i<size;i++)
  {
    beta = beta + res[i]*res[i];
    bres = bres + phi0[i]*phi0[i];
  }
  beta = sqrt(beta);
  bres = sqrt(bres);
  
  for (i=0;i<size;i++)
  {
    q[0][i] = res[i]/beta;
    //printf("%.8f ", q[0][i]);
  }
  //printf("\n");
  
  
  // Begin Arnoldi iterations. Taken from
  // https://en.wikipedia.org/wiki/Arnoldi_iteration
  for (iter = 1; iter < max_iter; iter++)
  {
    // Allocate a new q, h.
    q[iter] = new double[size];
    //q[iter] = (double*)malloc(size*sizeof(double));
    // Column 'i' of the Hessenberg matrix has
    // 'i+1' non-zero components.
    // To be fancier, we could do malloc((iter+1)...),
    // but then we have to do book-keeping everywhere.
    // That's why I kept having problems.
    h[iter-1] = new double[size];
    hTh[iter-1] = new double[size];
    //h[iter-1] = (double*)malloc(size*sizeof(double));
    //hTh[iter-1] = (double*)malloc(size*sizeof(double));
    for (j=0;j<size;j++)
    {
      q[iter][j] = h[iter-1][j] = hTh[iter-1][j] = 0.0;
    }
    
    // Compute q_{iter} = Aq_{iter-1}
    (*matrix_vector)(q[iter],q[iter-1],extra_info);
    
    // Perform an Arnoldi iteration. In some regards, this is just
    // a Gram-Schmidt process.
    for (j=0;j<iter;j++)
    {
      // Compute q[j] dot q[iter].
      h[iter-1][j] = 0.0;
      for (i=0;i<size;i++)
      {
        h[iter-1][j] = h[iter-1][j] + q[j][i]*q[iter][i]; //q* q
      }
      
      // Subtract off part of q[iter] along q[j].
      for (i=0;i<size;i++)
      {
        q[iter][i] = q[iter][i] - h[iter-1][j]*q[j][i];
      }
    } // Go to next existing q vector.
    
    // Normalize the new q[iter].
    h[iter-1][iter] = 0.0;
    for(i=0;i<size;i++)
    {
      h[iter-1][iter] = h[iter-1][iter] + q[iter][i]*q[iter][i];
    }
    h[iter-1][iter] = sqrt(h[iter-1][iter]);
    
    for(i=0;i<size;i++)
    {
      q[iter][i] = q[iter][i]/h[iter-1][iter];
      //printf("%.8f ", q[iter][i]);
    }
    //printf("\n");
    
    // This is formally the end of the Arnoldi iteration.
    // Right now, H is an (iter+1)x(iter) matrix stored
    // as its transpose.
    // h is iter x iter+1
    // For GMRES, we minimize || H y - beta e_1 ||_2 for y.
    // This becomes the normal equation: H^T H y = beta H^T e_1.
    // We solve this with gaussian elimination, but QR is apparently
    // a better way to do this.
    
    // Compute bhTy = beta H^T e_1, or the first column of H^T.
    //printf("bhTy: ");
    for(i=0;i<iter;i++)
    {
      bhTy[i] = beta*h[i][0];
      //printf("%.8f ", bhTy[i]);
    }
    //printf("\n");
    
    // Need to compute hTh, which is the iter x iter
    // normal matrix.
    // Recall h is iter x (iter+1)
    for (i=0;i<iter;i++)
    {
      for (j=0;j<iter;j++)
      {
        hTh[i][j] = 0.0;
        for (k=0;k<iter+1;k++)
        {
          hTh[i][j] = hTh[i][j] + h[i][k]*h[j][k];
        }
        //printf("%.8f ", hTh[i][j]);
      }
      //printf("\n");
    }
    gaussian_elimination(y, bhTy, hTh, iter);
    
    //printf("Exit gaussian elimination.\n"); fflush(stdout);
    
    // Alternatively, we can solve it with CG.
    // This is bad because we want an "exact" solution, which we
    // can only get in exact arithmetic. In practice, I didn't
    // get any benefit from using CG since the matrix is
    // so dense.
    
    //gmstr.size = iter;
    //minv_vector_cg(y, bhTy, iter, iter, eps*0.001, &gmres_hTh, (void*)(&gmstr));
    
    //double* meh = (double*)malloc(iter*sizeof(double));
    //gmres_hTh(meh, y, (void*)(&gmstr));
    //printf("meh: ");
    //for(i=0;i<iter;i++)
    //{
    //  printf("%.8f ", meh[i]);
    //}
    //printf("\n");
    
    // Use the q's to prolong y to the solution, put in tmp.
    for (i=0;i<size;i++)
    {
      tmp[i] = 0.0;
    }
    
    for (j=0;j<iter;j++)
    {
      for (i=0;i<size;i++)
      {
        tmp[i] = tmp[i]+q[j][i]*y[j];
      }
    }
    
    for (i=0;i<size;i++)
    {
      tmp2[i] = phi[i]+tmp[i];
    }
    
    // Check the residual. In principle, we can do this by looking at
    // the subspace norm, but I'm not sure how to combine that with
    // the fact we have an initial guess. Have to think about that. 
    (*matrix_vector)(res,tmp2,extra_info);
    
    localres = 0.0;
    for (i=0;i<size;i++)
    {
      localres = localres+(phi0[i]-res[i])*(phi0[i]-res[i]);
    }
    localres=sqrt(localres);
    
    if (localres < eps*bres)
    {
      break;
    }
    
    //printf("GMRES: Iter %d Resid %.8e\n", iter, localres/bres);
    
    
  }
  
  for (i=0;i<size;i++)
  {
    phi[i] = tmp2[i];
  }
  
  if(iter == max_iter) {
    //printf("GMRES: Failed to converge iter = %d, rsq = %e\n", iter,localres);
    invif.success = false;
    //return 0;// Failed convergence 
  }
  else
  {
     invif.success = true;
     //printf("GMRES: Converged in %d iterations.\n", iter);
  }
  
  // Check q's. Passes.
  /*
  for (i=0;i<iter;i++)
  {
    for (j=0;j<iter;j++)
    {
      bres = 0.0;
      for (k=0;k<size;k++)
      {
        bres = bres + q[i][k]*q[j][k];
      }
      printf("%.8f ", bres);
    }
    printf("\n");
  }
  */
  
  // Check A Q_n = Q_{n+1} \tilde H_n.
  // This requires a size x (iter-1) space.
  // We can check this column by column (of size 'size').
  // Works!
  /*for (i=0;i<iter-1;i++)
  {
    // Column of A Q_n
    (*matrix_vector)(res, q[i], extra_info);
    
    for (j=0;j<size;j++)
    {
      for (k=0;k<iter;k++)
      {
        res[j] = res[j] - q[k][j]*h[i][k];
      }
    }
    
    bres = 0.0;
    // Find magnitude of res.
    for (j=0;j<size;j++)
    {
      bres = bres + res[j]*res[j];
    }
    printf("Mag column %d: %.8f\n", i, sqrt(bres));
  }*/
    
    
  
  
  // Clean up.
  for (i=0;i<iter-1;i++)
  {
    delete[] h[i];
    delete[] q[i];
    delete[] hTh[i];
    //free(h[i]);
    //free(q[i]);
    //free(hTh[i]);
  }
  delete[] q[iter-1];
  delete[] h;
  delete[] q;
  delete[] res;
  delete[] tmp;
  delete[] tmp2;
  delete[] y;
  delete[] bhTy;
  //free(q[iter-1]);
  //free(h);
  //free(q);
  //free(res);
  //free(tmp);
  //free(tmp2);
  //free(y);
  //free(bhTy);
  
  invif.resSq = localres*localres;
  invif.iter = iter;
  if (invif.success == false) { invif.iter--; } // For loop has an extra incr. at the end.
  invif.name = "GMRES";
  return invif; // Convergence 

} 


// Performs GMRES with restarts when restart_freq is hit.
// This may be sloppy, but it works.
inversion_info minv_vector_gmres_restart(double  *phi, double  *phi0, int size, int max_iter, double res, int restart_freq, void (*matrix_vector)(double*,double*,void*), void* extra_info)
{
  int iter; // counts total number of iterations.
  inversion_info invif;

  iter = 0;  
  do
  {
    invif = minv_vector_gmres_norestart(phi, phi0, size, restart_freq, res, matrix_vector, extra_info);
    iter += invif.iter;
  }
  while (iter < max_iter && invif.success == false && sqrt(invif.resSq) > res);
  
  invif.iter = iter;
  stringstream ss;
  ss << "GMRES(" << restart_freq << ")";
  invif.name = ss.str();
  // invif.resSq is good.
  if (sqrt(invif.resSq) > res)
  {
    invif.success = false;
  }
  else
  {
    invif.success = true;
  }
  
  return invif;
}

inversion_info minv_vector_gmres_norestart(complex<double>  *phi, complex<double>  *phi0, int size, int max_iter, double eps, void (*matrix_vector)(complex<double>*,complex<double>*,void*), void* extra_info)
{
// GMRES solutions to Mphi = b 
//  see https://en.wikipedia.org/wiki/Generalized_minimal_residual_method
// A clearer article: 
// http://www.math.iit.edu/~fass/477577_Chapter_14.pdf
  
  // Solve a counting issue...
  max_iter++;
  
  // Initialize vectors.
  complex<double> **q, **h; // Q matrix, Hessenberg matrix.
  complex<double> **hTh; // stores normal eqn matrix for hTh. 
  complex<double> *res, *tmp, *tmp2;
  complex<double> *y, *bhTy; // Holds the reduced solution. 
  double bres, localres; // norm of b, curr residual.
  double beta; // beta is the norm of res = b-Ax_0.
  //gmres_struct gmstr; // Passed to solve lls problem, no longer needed: we use gaussian elim.
  int i, j, k, iter;
  inversion_info invif;
  
  // Prepare q, h. These get larger on each iteration.
  q = new complex<double>*[max_iter];
  h = new complex<double>*[max_iter-1];
  hTh = new complex<double>*[max_iter-1];
  //q = (double**)malloc(max_iter*sizeof(double*));
  //h = (double**)malloc((max_iter-1)*sizeof(double*));
  //hTh = (double**)malloc((max_iter-1)*sizeof(double*));
  //gmstr.h = h; // Depreciated, no longer use CG for normal subspace.
  
  // Prepare space for the residual, and a temporary vector.
  res = new complex<double>[size];
  tmp = new complex<double>[size];
  tmp2 = new complex<double>[size];
  y = new complex<double>[max_iter];
  bhTy = new complex<double>[max_iter];
  //res = (double*)malloc(size*sizeof(double));
  //tmp = (double*)malloc(size*sizeof(double));
  //tmp2 = (double*)malloc(size*sizeof(double));
  //y = (double*)malloc(max_iter*sizeof(double));
  //bhTy = (double*)malloc(max_iter*sizeof(double));
  
  // Compute the residual. Ax' = res gives the solution phi+x'.
  (*matrix_vector)(tmp, phi, extra_info);
  for (i=0;i<size;i++) // res = b - Ax_0
  {
    res[i] = phi0[i]-tmp[i]; 
  }
  
  // Initialize the first q vector. q is res/||res||
  q[0] = new complex<double>[size];
  //q[0] = (double*)malloc(size*sizeof(double));
  
  beta = sqrt(norm2sq<double>(res, size));
  bres = sqrt(norm2sq<double>(phi0, size));
  
  for (i=0;i<size;i++)
  {
    q[0][i] = res[i]/beta;
    //printf("%.8f ", q[0][i]);
  }
  //printf("\n");
  
  
  // Begin Arnoldi iterations. Taken from
  // https://en.wikipedia.org/wiki/Arnoldi_iteration
  for (iter = 1; iter < max_iter; iter++)
  {
    // Allocate a new q, h.
    q[iter] = new complex<double>[size];
    //q[iter] = (double*)malloc(size*sizeof(double));
    // Column 'i' of the Hessenberg matrix has
    // 'i+1' non-zero components.
    // To be fancier, we could do malloc((iter+1)...),
    // but then we have to do book-keeping everywhere.
    // That's why I kept having problems.
    h[iter-1] = new complex<double>[size];
    hTh[iter-1] = new complex<double>[size];
    //h[iter-1] = (double*)malloc(size*sizeof(double));
    //hTh[iter-1] = (double*)malloc(size*sizeof(double));
    for (j=0;j<size;j++)
    {
      q[iter][j] = h[iter-1][j] = hTh[iter-1][j] = 0.0;
    }
    
    // Compute q_{iter} = Aq_{iter-1}
    (*matrix_vector)(q[iter],q[iter-1],extra_info);
    
    // Perform an Arnoldi iteration. In some regards, this is just
    // a Gram-Schmidt process.
    for (j=0;j<iter;j++)
    {
      // Compute q[j] dot q[iter].
      h[iter-1][j] = dot<double>(q[j],q[iter], size);
      
      // Subtract off part of q[iter] along q[j].
      for (i=0;i<size;i++)
      {
        q[iter][i] = q[iter][i] - h[iter-1][j]*q[j][i];
      }
    } // Go to next existing q vector.
    
    // Normalize the new q[iter].
    h[iter-1][iter] = sqrt(norm2sq<double>(q[iter], size));
    for(i=0;i<size;i++)
    {
      q[iter][i] = q[iter][i]/h[iter-1][iter];
      //printf("%.8f ", q[iter][i]);
    }
    //printf("\n");
    
    // This is formally the end of the Arnoldi iteration.
    // Right now, H is an (iter+1)x(iter) matrix stored
    // as its transpose.
    // h is iter x iter+1
    // For GMRES, we minimize || H y - beta e_1 ||_2 for y.
    // This becomes the normal equation: H^T H y = beta H^T e_1.
    // We solve this with gaussian elimination, but QR is apparently
    // a better way to do this.
    
    // Compute bhTy = beta H^T e_1, or the first column of H^T.
    //printf("bhTy: ");
    for(i=0;i<iter;i++)
    {
      bhTy[i] = beta*h[i][0];
      //printf("%.8f ", bhTy[i]);
    }
    //printf("\n");
    
    // Need to compute hTh, which is the iter x iter
    // normal matrix.
    // Recall h is iter x (iter+1)
    for (i=0;i<iter;i++)
    {
      for (j=0;j<iter;j++)
      {
        hTh[i][j] = dot<double>(h[i],h[j], iter+1);
        //printf("%.8f ", hTh[i][j]);
      }
      //printf("\n");
    }
    gaussian_elimination(y, bhTy, hTh, iter);
    
    //printf("Exit gaussian elimination.\n"); fflush(stdout);
    
    // Alternatively, we can solve it with CG.
    // This is bad because we want an "exact" solution, which we
    // can only get in exact arithmetic. In practice, I didn't
    // get any benefit from using CG since the matrix is
    // so dense.
    
    //gmstr.size = iter;
    //minv_vector_cg(y, bhTy, iter, iter, eps*0.001, &gmres_hTh, (void*)(&gmstr));
    
    //double* meh = (double*)malloc(iter*sizeof(double));
    //gmres_hTh(meh, y, (void*)(&gmstr));
    //printf("meh: ");
    //for(i=0;i<iter;i++)
    //{
    //  printf("%.8f ", meh[i]);
    //}
    //printf("\n");
    
    // Use the q's to prolong y to the solution, put in tmp.
    for (i=0;i<size;i++)
    {
      tmp[i] = 0.0;
    }
    
    for (j=0;j<iter;j++)
    {
      for (i=0;i<size;i++)
      {
        tmp[i] = tmp[i]+q[j][i]*y[j];
      }
    }
    
    for (i=0;i<size;i++)
    {
      tmp2[i] = phi[i]+tmp[i];
    }
    
    // Check the residual. In principle, we can do this by looking at
    // the subspace norm, but I'm not sure how to combine that with
    // the fact we have an initial guess. Have to think about that. 
    (*matrix_vector)(res,tmp2,extra_info);
    
    localres = 0.0;
    for (i=0;i<size;i++)
    {
      localres = localres+pow(abs(phi0[i]-res[i]),2);
    }
    localres=sqrt(localres);
    
    if (localres < eps*bres)
    {
      break;
    }
    
    //printf("GMRES: Iter %d Resid %.8e\n", iter, localres/bres);
    
    
  }
  
  for (i=0;i<size;i++)
  {
    phi[i] = tmp2[i];
  }
  
  if(iter == max_iter) {
    //printf("GMRES: Failed to converge iter = %d, rsq = %e\n", iter,localres);
    invif.success = false;
    //return 0;// Failed convergence 
  }
  else
  {
     invif.success = true;
     //printf("GMRES: Converged in %d iterations.\n", iter);
  }
  
  // Check q's. Passes.
  /*
  complex<double> tmp3;
  for (i=0;i<iter;i++)
  {
    for (j=0;j<iter;j++)
    {
      tmp3 = 0.0;
      for (k=0;k<size;k++)
      {
        tmp3 = tmp3 + conj(q[i][k])*q[j][k];
      }
      cout << tmp3 << " ";
    }
    cout << "\n";
  }*/
  
  
  // Check A Q_n = Q_{n+1} \tilde H_n.
  // This requires a size x (iter-1) space.
  // We can check this column by column (of size 'size').
  // Works!
  /*
  for (i=0;i<iter-1;i++)
  {
    // Column of A Q_n
    (*matrix_vector)(res, q[i], extra_info);
    
    for (j=0;j<size;j++)
    {
      for (k=0;k<iter;k++)
      {
        res[j] = res[j] - q[k][j]*h[i][k];
      }
    }
    
    complex<double> tmp3 = 0.0;
    // Find magnitude of res.
    for (j=0;j<size;j++)
    {
      tmp3 = tmp3 + conj(res[j])*res[j];
    }
    cout << "Mag column " << i << ": " << sqrt(tmp3) << "\n";
  }
  */
    
  
  
  // Clean up.
  for (i=0;i<iter-1;i++)
  {
    delete[] h[i];
    delete[] q[i];
    delete[] hTh[i];
    //free(h[i]);
    //free(q[i]);
    //free(hTh[i]);
  }
  delete[] q[iter-1];
  delete[] h;
  delete[] q;
  delete[] res;
  delete[] tmp;
  delete[] tmp2;
  delete[] y;
  delete[] bhTy;
  //free(q[iter-1]);
  //free(h);
  //free(q);
  //free(res);
  //free(tmp);
  //free(tmp2);
  //free(y);
  //free(bhTy);
  
  invif.resSq = localres*localres;
  invif.iter = iter;
  if (invif.success == false) { invif.iter--; } // For loop has an extra incr. at the end.
  invif.name = "GMRES";
  return invif; // Convergence 

} 


// Performs GMRES with restarts when restart_freq is hit.
// This may be sloppy, but it works.
inversion_info minv_vector_gmres_restart(complex<double>  *phi, complex<double>  *phi0, int size, int max_iter, double res, int restart_freq, void (*matrix_vector)(complex<double>*,complex<double>*,void*), void* extra_info)
{
  int iter; // counts total number of iterations.
  inversion_info invif;

  iter = 0;  
  do
  {
    invif = minv_vector_gmres_norestart(phi, phi0, size, restart_freq, res, matrix_vector, extra_info);
    iter += invif.iter;
  }
  while (iter < max_iter && invif.success == false && sqrt(invif.resSq) > res);
  
  invif.iter = iter;
  stringstream ss;
  ss << "GMRES(" << restart_freq << ")";
  invif.name = ss.str();
  // invif.resSq is good.
  if (sqrt(invif.resSq) > res)
  {
    invif.success = false;
  }
  else
  {
    invif.success = true;
  }
  
  return invif;
}


// This is no longer necessary, we now use
// explicit gaussian elimination instead of CG.
// Solves H^T H phi = phi0. 
// extra_data is really a gmres_struct
// containing h = H^T which is n by (n+1),
// and the size.
/*void gmres_hTh(double* phi, double* phi0, void* extra_data)
{
  gmres_struct gmstr = *(gmres_struct*)(extra_data);
  double* yinter; // Holds the intermediate of H^T H.
  int i,j;
  
  yinter = (double*)malloc((gmstr.size+1)*sizeof(double));
  for (i=0;i<gmstr.size+1;i++)
  {
    yinter[i] = 0.0;
  }
  for (i=0;i<gmstr.size;i++)
  {
    phi[i] = 0.0;
  }
  
  // Do the H phi0 multiply.
  for (i=0;i<gmstr.size+1;i++)
  {
    for (j=0;j<gmstr.size;j++)
    {
      yinter[i] = yinter[i] + gmstr.h[j][i]*phi0[j];
    }
  }
  
  // Do the H^T yinter multiply.
  for (i=0;i<gmstr.size;i++)
  {
    for (j=0;j<gmstr.size+1;j++)
    {
      phi[i] = phi[i] + gmstr.h[i][j]*yinter[j];
    }
  }
  
  free(yinter);
}*/



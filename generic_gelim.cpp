// Fri Apr 29 14:25:39 EDT 2016
// Evan S Weinberg 2016
// A quick implementation of gaussian elimination with row pivoting.

#include <iostream>
#include <cmath>
#include <string>
#include <sstream>
#include <complex>

#include "generic_inverters.h"
#include "generic_vector.h"

using namespace std;

// Perform gaussian elimination on a matrix to solve Ax=b
// If gaussian elimination fails, return 0, else return 1.
int gaussian_elimination(double* x, double* b, double** matrix, int size)
{
  // Initialize variables.
  int i,j,k;
  double** grown_matrix;
  double pivot_space; // For pivot.
  double max_val = 0.0;
  int max_index = -1;
  
  // Declare size of matrix.
  grown_matrix = new double*[size];
  
  for (i=0;i<size;i++)
  {
    grown_matrix[i] = new double[size+1];
  }
  
  
  // Copy things over.
  for (i=0;i<size;i++)
  {
    for (j=0;j<size;j++)
    {
      grown_matrix[i][j] = matrix[i][j];
    }
    grown_matrix[i][size] = b[i];
  }
  
  // Cool! Begin gaussian elimination.
  // Iterate over all rows.
  for (i=0;i<size;i++)
  {
    max_val = 0.0; max_index = -1;
    // First, pivot over self and all lower rows.
    for (j=i;j<size;j++)
    {
      if (abs(grown_matrix[j][i]) > max_val)
      {
        max_index = j;
        max_val = abs(grown_matrix[j][i]);
      }
    }
    
    //printf("Max: %d %.8f\n", max_index, max_val); fflush(stdout);
    if (max_index == -1) // everything's 0!
    {
      return 0; 
    }
    
    // Put maximal row in pivot location.
    // Switch i, max_index.
    if (max_index != i)
    {
      for(j=i;j<size+1;j++)
      {
        pivot_space = grown_matrix[i][j];
        grown_matrix[i][j] = grown_matrix[max_index][j];
        grown_matrix[max_index][j] = pivot_space;
      }
      //printf("Pivoted.\n"); fflush(stdout);
    }
    
    
    // Good! We've safely pivoted. Now, normalize the top row.
    for (j=i+1;j<size+1;j++)
    {
      grown_matrix[i][j] = grown_matrix[i][j]/grown_matrix[i][i];
    }
    grown_matrix[i][i] = 1.0;
    
    // Eliminate the top row from all other rows.
    // This part can get parallelized!
    for (j=0;j<size;j++)
    {
      if (j == i) continue;
      for (k=i+1;k<size+1;k++)
      {
        grown_matrix[j][k] = grown_matrix[j][k] - grown_matrix[j][i]*grown_matrix[i][k];
      }
      grown_matrix[j][i] = 0.0;
    }
    
  } // end loop over rows = i.
  //printf("Exited reduction.\n"); fflush(stdout);
  
  // Copy result in!
  for (i=0;i<size;i++)
  {
    x[i] = grown_matrix[i][size];
  }
  
  //printf("Copied solution.\n"); fflush(stdout);
  
  // Free
  for (i=0;i<size;i++)
  {
    delete[] grown_matrix[i];
  }
  delete[] grown_matrix;
  
  return 1;
}

int gaussian_elimination(complex<double>* x, complex<double>* b, complex<double>** matrix, int size)
{
  // Initialize variables.
  int i,j,k;
  complex<double>** grown_matrix;
  complex<double> pivot_space; // For pivot.
  double max_val = 0.0;
  int max_index = -1;
  
  // Declare size of matrix.
  grown_matrix = new complex<double>*[size];
  for (i=0;i<size;i++)
  {
    grown_matrix[i] = new complex<double>[size+1];
  }
  
  // Copy things over.
  for (i=0;i<size;i++)
  {
    for (j=0;j<size;j++)
    {
      grown_matrix[i][j] = matrix[i][j];
    }
    grown_matrix[i][size] = b[i];
  }
  
  // Cool! Begin gaussian elimination.
  // Iterate over all rows.
  for (i=0;i<size;i++)
  {
    max_val = 0.0; max_index = -1;
    // First, pivot over self and all lower rows.
    for (j=i;j<size;j++)
    {
      if (abs(grown_matrix[j][i]) > max_val)
      {
        max_index = j;
        max_val = abs(grown_matrix[j][i]);
      }
    }
    
    //printf("Max: %d %.8f\n", max_index, max_val); fflush(stdout);
    if (max_index == -1) // everything's 0!
    {
      #ifdef VERBOSE_WARN
      cout << "Singular matrix in Gaussian elimination.\n" << flush;
      #endif
      return 0;
    }
    
    // Put maximal row in pivot location.
    // Switch i, max_index.
    if (max_index != i)
    {
      for(j=i;j<size+1;j++)
      {
        pivot_space = grown_matrix[i][j];
        grown_matrix[i][j] = grown_matrix[max_index][j];
        grown_matrix[max_index][j] = pivot_space;
      }
      //printf("Pivoted.\n"); fflush(stdout);
    }
    
    
    // Good! We've safely pivoted. Now, normalize the top row.
    for (j=i+1;j<size+1;j++)
    {
      grown_matrix[i][j] = grown_matrix[i][j]/grown_matrix[i][i];
    }
    grown_matrix[i][i] = 1.0;
    
    // Eliminate the top row from all other rows.
    // This part can get parallelized!
    for (j=0;j<size;j++)
    {
      if (j == i) continue;
      for (k=i+1;k<size+1;k++)
      {
        grown_matrix[j][k] = grown_matrix[j][k] - grown_matrix[j][i]*grown_matrix[i][k];
      }
      grown_matrix[j][i] = 0.0;
    }
    
  } // end loop over rows = i.
  //printf("Exited reduction.\n"); fflush(stdout);
  
  // Copy result in!
  for (i=0;i<size;i++)
  {
    x[i] = grown_matrix[i][size];
  }
  
  //printf("Copied solution.\n"); fflush(stdout);
  
  // Free
  for (i=0;i<size;i++)
  {
    free(grown_matrix[i]);
  }
  free(grown_matrix);
  
  return 1;
}


// This code tests the gaussian elimination routine.
/*
void test_gaussian()
{
   int m,n;
   int lsize = 20;
   double** test;
   double* b;
   double* x;
   double* b2;
   
   test = (double**)malloc(lsize*sizeof(double*));
   for (m=0;m<lsize;m++)
   {
     test[m] = (double*)malloc(lsize*sizeof(double));
   }
   b = (double*)malloc(lsize*sizeof(double));
   b2 = (double*)malloc(lsize*sizeof(double));
   x = (double*)malloc(lsize*sizeof(double));
   
   for (m=0;m<lsize;m++)
   {
     test[m][m] = 3;
     test[m][(m+1)%lsize] = -1;
     test[m][(m+3)%lsize] = -6;
     test[m][(m-1+lsize)%lsize] = -1;
     test[m][(m-3+lsize)%lsize] = -6;
   }
   
   b[0] = 1;
   gaussian_elimination(x,b,test,lsize);
   
   for (m=0;m<lsize;m++)
   {
     printf("%f ", x[m]);
   }
   printf("\n");
   
   for (m=0;m<lsize;m++)
   {
     b2[m] = 0.0;
     for (n=0;n<lsize;n++)
     {
       b2[m] = b2[m] + test[m][n]*x[n];
     }
   }
   
   for (m=0;m<lsize;m++)
   {
     printf("%f ", b[m]-b2[m]);
   }
   printf("\n");
   

   return;
}
*/



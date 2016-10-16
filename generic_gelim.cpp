// Fri Apr 29 14:25:39 EDT 2016
// Evan S Weinberg 2016
// A quick implementation of gaussian elimination with row pivoting.

#include <iostream>
#include <cmath>
#include <string>
#include <sstream>
#include <complex>

#include "generic_vector.h"

#include "generic_gelim.h"

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

// Perform gaussian elimination on a matrix to solve Ax=b, where there are multiple b.
// If gaussian elimination fails, return 0, else return 1.
int gaussian_elimination_multi_rhs(double** x, double** b, double** matrix, int n_rhs, int size)
{
  // Initialize variables.
  int i,j,k;
  double** grown_matrix;
  double* pivot_space; // For pivot.
  double max_val = 0.0;
  int max_index = -1;
  
  // Declare size of matrix.
  grown_matrix = new double*[size];
  
  for (i=0;i<size;i++)
  {
    grown_matrix[i] = new double[size+n_rhs];
  }
  
  
  // Copy things over.
  for (i=0;i<size;i++)
  {
    for (j=0;j<size;j++)
    {
      grown_matrix[i][j] = matrix[i][j];
    }
    for (j=0;j<n_rhs;j++)
    {
      grown_matrix[i][size+j] = b[j][i];
    }
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
      pivot_space = grown_matrix[i];
      grown_matrix[i] = grown_matrix[max_index];
      grown_matrix[max_index] = pivot_space;
      //printf("Pivoted.\n"); fflush(stdout);
    }
    
    
    // Eliminate the top row from lower rows.
    for (j=i+1;j<size;j++)
    {
      double factor = -grown_matrix[j][i]/grown_matrix[i][i];
      for (k=i;k<size+n_rhs;k++)
      {
        grown_matrix[j][k] = grown_matrix[j][k] + grown_matrix[i][k]*factor;
      }
      //grown_matrix[j][i] = 0.0;
    }
    /*for (k=i+1;k<size+n_rhs;k++)
    {
      grown_matrix[i][k] /= grown_matrix[i][i];
    }
    grown_matrix[i][i] = 1.0;*/
    
    for (j=0;j<size;j++)
    {
      for (k=0;k<size+n_rhs;k++)
      {
        cout << grown_matrix[j][k] << " ";
      }
      cout << "\n";
    }
    cout << "\n";
    
    /*
    // Good! We've safely pivoted. Now, normalize the top row.
    for (j=i+1;j<size+n_rhs;j++)
    {
      grown_matrix[i][j] = grown_matrix[i][j]/grown_matrix[i][i];
    }
    grown_matrix[i][i] = 1.0;
    
    // Eliminate the top row from all other rows.
    // This part can get parallelized!
    for (j=0;j<size;j++)
    {
      if (j == i) continue;
      for (k=i+1;k<size+n_rhs;k++)
      {
        grown_matrix[j][k] = grown_matrix[j][k] - grown_matrix[j][i]*grown_matrix[i][k];
      }
      grown_matrix[j][i] = 0.0;
    }*/
    
  } // end loop over rows = i.
  //printf("Exited reduction.\n"); fflush(stdout);
  
  // Now we back substitute!
  for (i=size-1;i>=0;i--)
  {
    for (j=i-1;j>=0;j--)
    {
      double factor = -grown_matrix[j][i]/grown_matrix[i][i];
      for (k=i;k<size+n_rhs;k++)
      {
        grown_matrix[j][k] = grown_matrix[j][k]+grown_matrix[i][k]*factor;
      }
      //grown_matrix[j][i] = 0.0;
    }
    for (k=i+1;k<size+n_rhs;k++)
    {
      grown_matrix[i][k] /= grown_matrix[i][i];
    }
    grown_matrix[i][i] = 1.0;
    
    for (j=0;j<size;j++)
    {
      for (k=0;k<size+n_rhs;k++)
      {
        cout << grown_matrix[j][k] << " ";
      }
      cout << "\n";
    }
    cout << "\n";
  }
  
  // Copy result in!
  for (k=0;k<n_rhs;k++)
  {
    for (i=0;i<size;i++)
    {
      x[k][i] = grown_matrix[i][size+k];
    }
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

int gaussian_elimination_multi_rhs(complex<double>** x, complex<double>** b, complex<double>** matrix, int n_rhs, int size)
{
  // Initialize variables.
  int i,j,k;
  complex<double>** grown_matrix;
  complex<double>* pivot_space; // For pivot.
  double max_val = 0.0;
  int max_index = -1;
  
  // Declare size of matrix.
  grown_matrix = new complex<double>*[size];
  
  for (i=0;i<size;i++)
  {
    grown_matrix[i] = new complex<double>[size+n_rhs];
  }
  
  
  // Copy things over.
  for (i=0;i<size;i++)
  {
    for (j=0;j<size;j++)
    {
      grown_matrix[i][j] = matrix[i][j];
    }
    for (j=0;j<n_rhs;j++)
    {
      grown_matrix[i][size+j] = b[j][i];
    }
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
      pivot_space = grown_matrix[i];
      grown_matrix[i] = grown_matrix[max_index];
      grown_matrix[max_index] = pivot_space;
      //printf("Pivoted.\n"); fflush(stdout);
    }
    
    
    // Eliminate the top row from lower rows.
    for (j=i+1;j<size;j++)
    {
      complex<double> factor = -grown_matrix[j][i]/grown_matrix[i][i];
      for (k=i;k<size+n_rhs;k++)
      {
        grown_matrix[j][k] = grown_matrix[j][k] + grown_matrix[i][k]*factor;
      }
      //grown_matrix[j][i] = 0.0;
    }
    /*for (k=i+1;k<size+n_rhs;k++)
    {
      grown_matrix[i][k] /= grown_matrix[i][i];
    }
    grown_matrix[i][i] = 1.0;*/
    
    for (j=0;j<size;j++)
    {
      for (k=0;k<size+n_rhs;k++)
      {
        cout << grown_matrix[j][k] << " ";
      }
      cout << "\n";
    }
    cout << "\n";
    
    /*
    // Good! We've safely pivoted. Now, normalize the top row.
    for (j=i+1;j<size+n_rhs;j++)
    {
      grown_matrix[i][j] = grown_matrix[i][j]/grown_matrix[i][i];
    }
    grown_matrix[i][i] = 1.0;
    
    // Eliminate the top row from all other rows.
    // This part can get parallelized!
    for (j=0;j<size;j++)
    {
      if (j == i) continue;
      for (k=i+1;k<size+n_rhs;k++)
      {
        grown_matrix[j][k] = grown_matrix[j][k] - grown_matrix[j][i]*grown_matrix[i][k];
      }
      grown_matrix[j][i] = 0.0;
    }*/
    
  } // end loop over rows = i.
  //printf("Exited reduction.\n"); fflush(stdout);
  
  // Now we back substitute!
  for (i=size-1;i>=0;i--)
  {
    for (j=i-1;j>=0;j--)
    {
      complex<double> factor = -grown_matrix[j][i]/grown_matrix[i][i];
      for (k=i;k<size+n_rhs;k++)
      {
        grown_matrix[j][k] = grown_matrix[j][k]+grown_matrix[i][k]*factor;
      }
      //grown_matrix[j][i] = 0.0;
    }
    for (k=i+1;k<size+n_rhs;k++)
    {
      grown_matrix[i][k] /= grown_matrix[i][i];
    }
    grown_matrix[i][i] = 1.0;
    
    for (j=0;j<size;j++)
    {
      for (k=0;k<size+n_rhs;k++)
      {
        cout << grown_matrix[j][k] << " ";
      }
      cout << "\n";
    }
    cout << "\n";
  }
  
  // Copy result in!
  for (k=0;k<n_rhs;k++)
  {
    for (i=0;i<size;i++)
    {
      x[k][i] = grown_matrix[i][size+k];
    }
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

// Compute a matrix inverse using gaussian elimination (uses multirhs under the hood). minv and matrix can be the same pointer.
int gaussian_elimination_matrix_inverse(double** minv, double** matrix, int size)
{
  int i,j;
  double tmp;
  
  // Allocate space for unit vectors.
  double** identity_transpose = new double*[size]; // technically, 'x' and 'b' are stored as the transpose in the context of a matrix inversion.
  for (i = 0; i < size; i++)
  {
    identity_transpose[i] = new double[size]; zero<double>(identity_transpose[i], size); identity_transpose[i][i] = 1.0;
  }
  
  // minv = (matrix)^(-1) identity_transpose
  int retval = gaussian_elimination_multi_rhs(minv, identity_transpose, matrix, size, size);
  
  // Need to transpose minv.
  for (i = 0; i < size-1; i++)
  {
    for (j = i+1; j < size; j++)
    {
      tmp = minv[i][j];
      minv[i][j] = minv[j][i];
      minv[j][i] = tmp;
    }
  }
  
  for (j = 0; j < size; j++)
    {
      for (int k = 0; k < size; k++)
      {
        cout << minv[j][k] << " ";
      }
      cout << "\n";
    }
    cout << "\n";
  
  for (i = 0; i < size; i++)
  {
    delete[] identity_transpose[i];
  }
  delete[] identity_transpose;
  
  return retval;
}

int gaussian_elimination_matrix_inverse(complex<double>** minv, complex<double>** matrix, int size)
{
  int i,j;
  complex<double> tmp;
  
  // Allocate space for unit vectors.
  complex<double>** identity_transpose = new complex<double>*[size]; // technically, 'x' and 'b' are stored as the transpose in the context of a matrix inversion.
  for (i = 0; i < size; i++)
  {
    identity_transpose[i] = new complex<double>[size]; zero<double>(identity_transpose[i], size); identity_transpose[i][i] = 1.0;
  }
  
  // minv = (matrix)^(-1) identity_transpose
  int retval = gaussian_elimination_multi_rhs(minv, identity_transpose, matrix, size, size);
  
  // Need to transpose minv.
  for (i = 0; i < size-1; i++)
  {
    for (j = i+1; j < size; j++)
    {
      tmp = minv[i][j];
      minv[i][j] = minv[j][i];
      minv[j][i] = tmp;
    }
  }
  
  for (i = 0; i < size; i++)
  {
    delete[] identity_transpose[i];
  }
  delete[] identity_transpose;
  
  return retval;
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



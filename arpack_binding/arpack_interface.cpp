// ESW 2015-09-07 Get eigenvalues+eigenvectors!

#include <stdio.h>
#include <stdlib.h>
#include <complex>

using std::complex;

#include "arpack_extern.h"
#include "arpack_interface.h"

#include "generic_vector.h"

/************************
* Double Real Symmetric *
************************/

arpack_drs_t* arpack_drs_init(int maxn, int maxnev, int maxncv)
{
   // Temporary things to mathc dssimp.f.
   int ldv; 

   // Initialize and allocate structure.
   arpack_drs_t* ap_str;
   ap_str = new arpack_drs_t;
   
   // Set inputs.
   ap_str->maxn = maxn;
   ap_str->maxnev = maxnev;
   ap_str->maxncv = maxncv;
   
   // To look like dssimp.
   ldv = maxn;
   
   // Allocate fixed memory arrays.
   ap_str->v = new double[ldv*maxncv];
   ap_str->d = new double[maxncv];
   ap_str->workl = new double[maxncv*(maxncv+8)];
   ap_str->workd = new double[3*maxn];
   ap_str->resid = new double[maxn];
   ap_str->select = new int[maxn];
   
   // And we're allocated!
   ap_str->is_allocated = 1;
	
	return ap_str; 
}

arpack_solve_t arpack_drs_getev(arpack_drs_t* arpack_str, double* eval, double** evec, int n, int nev, int ncv, int maxitr, char* which, double tol, double sigma, void (*matrix_vector)(double*,double*,void*), void* extra_info)
{
   // To just copy and paste code.
   int maxn, maxnev, maxncv, ldv;
   double *v, *workl, *workd, *d, *resid;
   int *select;
   int iparam[11];
   int ipntr[11];
   
   char bmat;
   char* hwmny; // ESW addition: how many eigenvalues to get.
   int ido, lworkl, info, ierr, ishfts, mode1, i, j;
   double zero;
   arpack_solve_t info_solve; // Gets returned.
   int rvec; // Actually a boolean!

   // Initialize values that aren't passed!
   maxn = arpack_str->maxn;
   maxnev = arpack_str->maxnev;
   maxncv = arpack_str->maxncv;
   ldv = maxn;
   
   v = arpack_str->v;
   workl = arpack_str->workl;
   workd = arpack_str->workd;
   d = arpack_str->d;
   resid = arpack_str->resid;
   select = arpack_str->select;
   
   zero = 0.0;
   
   bmat = 'I'; // This is a standard problem, as opposed
               // to a generalized eigenvalue problem.
               
   hwmny = "All"; // Get all eigenvalues.
   
   // Various error checking.
   if (n > maxn)
   {
      printf(" Error: n > maxn.\n");
      info_solve.is_error = -1;
      return info_solve;
   }
   else if (nev > maxnev)
   {
      printf(" Error: nev > maxnev.\n");
      info_solve.is_error = -1;
      return info_solve;
   }
   else if (ncv > maxncv)
   {
      printf(" Error: ncv > maxncv.\n");
      info_solve.is_error = -1;
      return info_solve;
   }
   
   // Stopping rules + initial conditions before calling DSAUPD
   lworkl = ncv*(ncv+8); // Trusting arpack here.
   ido = 0; // This is the reverse communication parameter 
            // from DSAUPD. Each call changes the value of this
            // parameter, and based on its value something must
            // be done. Has to be set to 0 before first call.
   info = 0; // On first use, specifies starting vector. Setting it
             // to 0 means use a random initial vector for arnoldi
             // iterations. Non-zero: pass in starting vector to
             // the array "resid".
             
   // Specify the algorithm mode. 
   ishfts = 1; // use an exact shift strategy. check DSAUPD
               // documentation for what this means. There are
               // options here to shift the matrix. (PARAM(1))
   mode1 = 1; // Use mode 1 of DSAUPD. Check documentation! (PARAM(7))
   
   // Set up iparam.
   iparam[0] = ishfts;
   iparam[2] = maxitr; // On return, gives actual number of iters.
   iparam[3] = 1; // Alexei's code does this.
   iparam[6] = mode1; 
   
   printf("Got to reverse comm.\n"); fflush(stdout);
   
   // Main reverse communication loop!
   while (ido != 99) // 99 means we're complete!
   {
      // Call dsaupd!
      
      ARPACK(dsaupd)(&ido, &bmat, &n, which, &nev, &tol,
                       resid, &ncv, v, &ldv,
                       (int*)iparam, (int*)ipntr, workd,
                       workl, &lworkl, &info);
      
      // Check for errors!                   
      if (info != 0)
      {
         printf(" Error in ARPACK dsaupd: %d.\n", info);
         info_solve.is_error = -1;
         return info_solve;
      }
      
      // See if we need to iterate more!
      if (ido == -1 || ido == 1)
      {
         // Perform the y = A*x, where 'x' starts at workd[ipntr[0]-1]
         // and y starts at workd[ipntr[1]-1], where the -1 is
         // because C, not FORTRAN.
         
         double* rhs = &(workd[ipntr[0]-1]);
         double* lhs = &(workd[ipntr[1]-1]);
         
         // Call rhs = A*lhs. 
         (*matrix_vector)(lhs, rhs, extra_info);
      }
      
   
   } // Loop back.
   
   printf("Finished iterating.\n"); fflush(stdout);
   
   // We're good! Post procecess to get eigenvalues and, if desired,
   // eigenvectors by rvec = true.
   
   rvec = 1; 
   
   // Important: d is eigenvalues, v is eigenvectors.
   ARPACK(dseupd)(&rvec, hwmny, select, d,
                   v,  &ldv, &sigma, &bmat, &n, which,
                   &nev, &tol, resid, &ncv, v,
                   &ldv, (int*)iparam, (int*)ipntr, workd,
                   workl, &lworkl, &ierr);

   printf("Got eigens.\n"); fflush(stdout);

   // Check for errors!
   if (ierr != 0)
   {
      printf("Error in ARPACK dseupd: %d.\n", info);
      info_solve.is_error = -1;
      return info_solve;
   }
   else
   {
      // No errors!
      info_solve.is_error = 0;
   
      // Number of converged eigenvalues.
      info_solve.nconv = iparam[4];
      
      // Number of update iterations.
      info_solve.niter = iparam[2];
      
      // Number of OP*x.
      info_solve.nops = iparam[8];
      
      /*
      // Get the residual norm:
      // || A*x - lambda*x ||
      // for the nconv accurately computed eigenvals/vecs.
      // Put the residual here:
      double* resid_arr = (d+maxncv);
      
      // Loop over all good eigenvalues.
      for (int i=0;i<nconv;i++)
      {
         // Get pointer to start of i'th eigenvector.
         double *evec = (v+i*ldv);
         
         // If you want to print it out, uncomment this.
         for (int j=0;j<NDIM;j++)
         {
            printf("%d\t%.8e\n", i, evec[j]);
         }
         
         // Get A*x, where x is an eigenvector.
         matrix_vector((double*)ax, evec);
         
         resid_arr[i] = 0;
         // Get the norm squared.
         for (int j=0;j<NDIM;j++)
         {
            resid_arr[i] += (ax[j]-d[i]*evec[j])*(ax[j]-d[i]*evec[j]);
         }
         resid_arr[i] = sqrt(resid_arr[i])/d[i];
      }
      */
      
      // Copy eval, evec into place.
      // This is sloppy...
      for (i=0;i<nev;i++)
      {
         eval[i]=d[i];
		  copy<double>(evec[i], (v+i*ldv), n);
      }
      
      return info_solve;
      
   }
   
   // The code shouldn't get here, so something went wrong.
   info_solve.is_error = -1;
   return info_solve;
}

// And free it up.
void arpack_drs_free(arpack_drs_t** arpack_str)
{
   arpack_drs_t* ap_str = *arpack_str;
   if (ap_str != NULL && ap_str->is_allocated == 1)
   {
      delete[] ap_str->v;
      delete[] ap_str->d;
      delete[] ap_str->workl;
      delete[] ap_str->workd;
      delete[] ap_str->resid;
      delete[] ap_str->select;
      delete ap_str;
      ap_str = NULL;
   }
}

/*************************
* Double Complex General *
*************************/


arpack_dcn_t* arpack_dcn_init(int maxn, int maxnev, int maxncv)
{
   // Temporary things to mathc dssimp.f.
   int ldv; 

   // Initialize and allocate structure.
   arpack_dcn_t* ap_str;
   ap_str = new arpack_dcn_t; 
   
   // Set inputs.
   ap_str->maxn = maxn;
   ap_str->maxnev = maxnev;
   ap_str->maxncv = maxncv;
   
   // To look like dssimp.
   ldv = maxn;
   
   // Allocate fixed memory arrays.
   ap_str->v = new complex<double>[ldv*maxncv];
   ap_str->d = new complex<double>[maxncv];
   ap_str->workl = new complex<double>[(3*maxncv*maxncv+5*maxncv)];
   ap_str->workd = new complex<double>[3*maxn];
   ap_str->workev = new complex<double>[2*maxncv];
   ap_str->resid = new complex<double>[maxn];
   ap_str->rwork = new double[maxncv];
   ap_str->rd = new double[maxncv*3];
   ap_str->select = new int[maxn];
   
   // And we're allocated!
   ap_str->is_allocated = 1;
	
	return ap_str; 
}

arpack_solve_t arpack_dcn_getev(arpack_dcn_t* arpack_str, complex<double>* eval, complex<double>** evec, int n, int nev, int ncv, int maxitr, char* which, double tol, complex<double> sigma, void (*matrix_vector)(complex<double>*,complex<double>*,void*), void* extra_info)
{
   // To just copy and paste code.
   int maxn, maxnev, maxncv, ldv;
   complex<double> *v, *workl, *workd, *workev, *d, *resid;
   double *rwork, *rd;
   int *select;
   int iparam[11];
   int ipntr[14];
   
   char bmat;
   char hwmny; // ESW addition: how many eigenvalues to get.
   int ido, lworkl, info, ierr, ishfts, mode1, i, j;
   double zero;
   arpack_solve_t info_solve; // Gets returned.
   int rvec; // Actually a boolean!

   // Initialize values that aren't passed!
   maxn = arpack_str->maxn;
   maxnev = arpack_str->maxnev;
   maxncv = arpack_str->maxncv;
   ldv = maxn;
   
   v = arpack_str->v;
   workl = arpack_str->workl;
   workd = arpack_str->workd;
   d = arpack_str->d;
   workev = arpack_str->workev;
   rwork = arpack_str->rwork;
   rd = arpack_str->rd;
   resid = arpack_str->resid;
   select = arpack_str->select;
   
   zero = 0.0;
   
   bmat = 'I'; // This is a standard problem, as opposed
               // to a generalized eigenvalue problem.
               
   hwmny = 'A'; // Get all eigenvalues.
   
   // Various error checking.
   if (n > maxn)
   {
      printf(" Error: n > maxn.\n");
      info_solve.is_error = -1;
      return info_solve;
   }
   else if (nev > maxnev)
   {
      printf(" Error: nev > maxnev.\n");
      info_solve.is_error = -1;
      return info_solve;
   }
   else if (ncv > maxncv)
   {
      printf(" Error: ncv > maxncv.\n");
      info_solve.is_error = -1;
      return info_solve;
   }
   
   // Stopping rules + initial conditions before calling DSAUPD
   lworkl = 3*ncv*ncv+5*ncv; // Trusting arpack here.
   ido = 0; // This is the reverse communication parameter 
            // from DSAUPD. Each call changes the value of this
            // parameter, and based on its value something must
            // be done. Has to be set to 0 before first call.
   info = 0; // On first use, specifies starting vector. Setting it
             // to 0 means use a random initial vector for arnoldi
             // iterations. Non-zero: pass in starting vector to
             // the array "resid".
             
   // Specify the algorithm mode. 
   ishfts = 1; // use an exact shift strategy. check DSAUPD
               // documentation for what this means. There are
               // options here to shift the matrix. (PARAM(1))
   mode1 = 1; // Use mode 1 of DSAUPD. Check documentation! (PARAM(7))
   
   // Set up iparam.
   iparam[0] = ishfts;
   iparam[2] = maxitr; // On return, gives actual number of iters.
   iparam[3] = 1; // Alexei's code does this.
   iparam[6] = mode1; 
   
   printf("Got to reverse comm.\n"); fflush(stdout);
   
   // Main reverse communication loop!
   while (ido != 99) // 99 means we're complete!
   {
      // Call dsaupd!
      
      ARPACK(znaupd)(&ido, &bmat, &n, which, &nev, &tol,
                       resid, &ncv, v, &ldv,
                       (int*)iparam, (int*)ipntr, workd,
                       workl, &lworkl, (double*)rwork, &info);
      
      // Check for errors!                   
      if (info != 0)
      {
         printf(" Error in ARPACK znaupd: %d.\n", info);
         info_solve.is_error = -1;
         return info_solve;
      }
      
      // See if we need to iterate more!
      if (ido == -1 || ido == 1)
      {
         // Perform the y = A*x, where 'x' starts at workd[ipntr[0]-1]
         // and y starts at workd[ipntr[1]-1], where the -1 is
         // because C, not FORTRAN.
         
         complex<double>* rhs = &(workd[ipntr[0]-1]);
         complex<double>* lhs = &(workd[ipntr[1]-1]);
         
         // Call rhs = A*lhs. 
         
         (*matrix_vector)(lhs, rhs, extra_info);
         
      }
      
   
   } // Loop back.
   
   printf("Finished iterating.\n"); fflush(stdout);
   
   // We're good! Post procecess to get eigenvalues and, if desired,
   // eigenvectors by rvec = true.
   
   rvec = 1; 
   
   // Important: d is eigenvalues, v is eigenvectors.
   ARPACK(zneupd)(&rvec, &hwmny, select, d,
                   v,  &ldv, &sigma, (complex<double>*)workev,
				   &bmat, &n, which,
                   &nev, &tol, resid, &ncv, v,
                   &ldv, (int*)iparam, (int*)ipntr, workd,
                   workl, &lworkl, rwork, &ierr);

   printf("Got eigens.\n"); fflush(stdout);

   // Check for errors!
   if (ierr != 0)
   {
      printf("Error in ARPACK zneupd: %d.\n", info);
      info_solve.is_error = -1;
      return info_solve;
   }
   else
   {
      // No errors!
      info_solve.is_error = 0;
   
      // Number of converged eigenvalues.
      info_solve.nconv = iparam[4];
      
      // Number of update iterations.
      info_solve.niter = iparam[2];
      
      // Number of OP*x.
      info_solve.nops = iparam[8];
      
      /*
      // Get the residual norm:
      // || A*x - lambda*x ||
      // for the nconv accurately computed eigenvals/vecs.
      // Put the residual here:
      double* resid_arr = (d+maxncv);
      
      // Loop over all good eigenvalues.
      for (int i=0;i<nconv;i++)
      {
         // Get pointer to start of i'th eigenvector.
         complex<double> *evec = (v+i*ldv);
         
         // If you want to print it out, uncomment this.
         for (int j=0;j<NDIM;j++)
         {
            printf("%d\t%.8e\n", i, evec[j]);
         }
         
         // Get A*x, where x is an eigenvector.
         matrix_vector((complex<double>*)ax, evec);
         
         resid_arr[i] = 0;
         complex<double> z,zbar;
         
         // Get the norm squared.
         for (int j=0;j<NDIM;j++)
         {
            z = (ax[j]-d[i]*evec[j]);
            zbar = creal(z)-cimag(z)*I;
            resid_arr[i] += z*zbar;
         }
         z = d[i];
         zbar = creal(z)-cimag(z)*I;
         resid_arr[i] = sqrt(resid_arr[i])/sqrt(z*zbar);
      }
      */
      
      // Copy eval, evec into place.
      // This is sloppy...
      for (i=0;i<nev;i++)
      {
         eval[i]=d[i];
		 copy<double>(evec[i], (v+i*ldv), n);
      }
      
      return info_solve;
      
   }
   
   // The code shouldn't get here, so something went wrong.
   info_solve.is_error = -1;
   return info_solve;
}


	// A shifted versions of the index op for rayleigh quotients.
struct shift_struct
{
    void* extra_data;
    void (*op)(complex<double>* a, complex<double>* b, void* extra_data);
    complex<double> shift;
    int length; 
};

void shift_op(complex<double>* lhs, complex<double>* rhs, void* extra_data)
{
    shift_struct* sop = (shift_struct*)extra_data;
    (*sop->op)(lhs, rhs, sop->extra_data);
    for (int i = 0; i < sop->length; i++)
    {
        lhs[i] -= sop->shift*rhs[i];
    }
}

arpack_solve_t arpack_dcn_getev_sinv(arpack_dcn_t* arpack_str, complex<double>* eval, complex<double>** evec, int n, int nev, int ncv, int maxitr, char* which, double tol, complex<double> sigma, void (*matrix_vector)(complex<double>*,complex<double>*,void*), inversion_info (*minv_vector)(complex<double>*,complex<double>*,int,int,double,matrix_vector_p,void*,inversion_verbose_struct*), int maxitr_cg, double tol_cg ,void* extra_info)
{
   // To just copy and paste code.
   int maxn, maxnev, maxncv, ldv;
   complex<double> *v, *workl, *workd, *workev, *d, *resid;
   double *rwork, *rd;
   int *select;
   int iparam[11];
   int ipntr[14];
   
   char bmat;
   char hwmny; // ESW addition: how many eigenvalues to get.
   int ido, lworkl, info, ierr, ishfts, mode1, i, j;
   double zero;
   arpack_solve_t info_solve; // Gets returned.
   int rvec; // Actually a boolean!
	
   // Create a shift invert structure.
	shift_struct shifter;
	shifter.extra_data = extra_info;
	shifter.op = matrix_vector;
	shifter.shift = sigma;
	shifter.length = n; 

   // Initialize values that aren't passed!
   maxn = arpack_str->maxn;
   maxnev = arpack_str->maxnev;
   maxncv = arpack_str->maxncv;
   ldv = maxn;
   
   v = arpack_str->v;
   workl = arpack_str->workl;
   workd = arpack_str->workd;
   d = arpack_str->d;
   workev = arpack_str->workev;
   rwork = arpack_str->rwork;
   rd = arpack_str->rd;
   resid = arpack_str->resid;
   select = arpack_str->select;
   
   zero = 0.0;
   
   bmat = 'I'; // This is a standard problem, as opposed
               // to a generalized eigenvalue problem.
               
   hwmny = 'A'; // Get all eigenvalues.
   
   // Various error checking.
   if (n > maxn)
   {
      printf(" Error: n > maxn.\n");
      info_solve.is_error = -1;
      return info_solve;
   }
   else if (nev > maxnev)
   {
      printf(" Error: nev > maxnev.\n");
      info_solve.is_error = -1;
      return info_solve;
   }
   else if (ncv > maxncv)
   {
      printf(" Error: ncv > maxncv.\n");
      info_solve.is_error = -1;
      return info_solve;
   }
   
   // Stopping rules + initial conditions before calling DSAUPD
   lworkl = 3*ncv*ncv+5*ncv; // Trusting arpack here.
   ido = 0; // This is the reverse communication parameter 
            // from DSAUPD. Each call changes the value of this
            // parameter, and based on its value something must
            // be done. Has to be set to 0 before first call.
   info = 0; // On first use, specifies starting vector. Setting it
             // to 0 means use a random initial vector for arnoldi
             // iterations. Non-zero: pass in starting vector to
             // the array "resid".
             
   // Specify the algorithm mode. 
   ishfts = 1; // use an exact shift strategy. check DSAUPD
               // documentation for what this means. There are
               // options here to shift the matrix. (PARAM(1))
   mode1 = 3; // Use mode 3 of DSAUPD. Check documentation! (PARAM(7))
   
   // Set up iparam.
   iparam[0] = ishfts;
   iparam[2] = maxitr; // On return, gives actual number of iters.
   iparam[3] = 1; // Alexei's code does this.
   iparam[6] = mode1; 
   
   printf("Got to reverse comm.\n"); fflush(stdout);
   
   // Main reverse communication loop!
   while (ido != 99) // 99 means we're complete!
   {
      // Call dsaupd!
      
      ARPACK(znaupd)(&ido, &bmat, &n, which, &nev, &tol,
                       resid, &ncv, v, &ldv,
                       (int*)iparam, (int*)ipntr, workd,
                       workl, &lworkl, (double*)rwork, &info);
      
      // Check for errors!                   
      if (info != 0)
      {
         printf(" Error in ARPACK znaupd: %d.\n", info);
         info_solve.is_error = -1;
         return info_solve;
      }
      
      // See if we need to iterate more!
      if (ido == -1 || ido == 1)
      {
         // Perform the y = A*x, where 'x' starts at workd[ipntr[0]-1]
         // and y starts at workd[ipntr[1]-1], where the -1 is
         // because C, not FORTRAN.
         
         complex<double>* rhs = &(workd[ipntr[0]-1]);
         complex<double>* lhs = &(workd[ipntr[1]-1]);
         
         // Call rhs = (A-sigma*I)^-1 *lhs. 
         
		 
         (*minv_vector)(lhs, rhs, n, maxitr_cg, tol_cg, shift_op, (void*)&shifter,0);
         //minv_vector_bicgstab_shift(rhs, lhs, n, maxitr_cg, 1e-7, &matrix_vector, sigma, extra_info);
         //(*matrix_vector)(lhs, rhs, extra_info);
         
      }
      
   
   } // Loop back.
   
   printf("Finished iterating.\n"); fflush(stdout);
   
   // We're good! Post procecess to get eigenvalues and, if desired,
   // eigenvectors by rvec = true.
   
   rvec = 1; 
   
   // Important: d is eigenvalues, v is eigenvectors.
   ARPACK(zneupd)(&rvec, &hwmny, select, d,
                   v,  &ldv, &sigma, (complex<double>*)workev,
				   &bmat, &n, which,
                   &nev, &tol, resid, &ncv, v,
                   &ldv, (int*)iparam, (int*)ipntr, workd,
                   workl, &lworkl, rwork, &ierr);

   printf("Got eigens.\n"); fflush(stdout);

   // Check for errors!
   if (ierr != 0)
   {
      printf("Error in ARPACK zneupd: %d.\n", info);
      info_solve.is_error = -1;
      return info_solve;
   }
   else
   {
      // No errors!
      info_solve.is_error = 0;
   
      // Number of converged eigenvalues.
      info_solve.nconv = iparam[4];
      
      // Number of update iterations.
      info_solve.niter = iparam[2];
      
      // Number of OP*x.
      info_solve.nops = iparam[8];
      
      /*
      // Get the residual norm:
      // || A*x - lambda*x ||
      // for the nconv accurately computed eigenvals/vecs.
      // Put the residual here:
      double* resid_arr = (d+maxncv);
      
      // Loop over all good eigenvalues.
      for (int i=0;i<nconv;i++)
      {
         // Get pointer to start of i'th eigenvector.
         complex<double> *evec = (v+i*ldv);
         
         // If you want to print it out, uncomment this.
         for (int j=0;j<NDIM;j++)
         {
            printf("%d\t%.8e\n", i, evec[j]);
         }
         
         // Get A*x, where x is an eigenvector.
         matrix_vector((complex<double>*)ax, evec);
         
         resid_arr[i] = 0;
         complex<double> z,zbar;
         
         // Get the norm squared.
         for (int j=0;j<NDIM;j++)
         {
            z = (ax[j]-d[i]*evec[j]);
            zbar = creal(z)-cimag(z)*I;
            resid_arr[i] += z*zbar;
         }
         z = d[i];
         zbar = creal(z)-cimag(z)*I;
         resid_arr[i] = sqrt(resid_arr[i])/sqrt(z*zbar);
      }
      */
      
      // Copy eval, evec into place.
      // This is sloppy...
      for (i=0;i<nev;i++)
      {
         eval[i]=d[i];
		  copy<double>(evec[i], (v+i*ldv), n);
      }
      
      return info_solve;
      
   }
   
   // The code shouldn't get here, so something went wrong.
   info_solve.is_error = -1;
   return info_solve;
}


// And free it up.
void arpack_dcn_free(arpack_dcn_t** arpack_str)
{
   arpack_dcn_t* ap_str = *arpack_str;
   if (ap_str != NULL && ap_str->is_allocated == 1)
   {
      delete[] ap_str->v;
      delete[] ap_str->d;
      delete[] ap_str->workl;
      delete[] ap_str->workd;
	  delete[] ap_str->workev;
	  delete[] ap_str->rwork;
	  delete[] ap_str->rd;
      delete[] ap_str->resid;
      delete[] ap_str->select;
      delete ap_str;
      ap_str = NULL;
   }
}

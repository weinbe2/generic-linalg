// Based on the quda arpack interface Alexei supplied to me.

#ifndef ARPACK_EXTERN
#define ARPACK_EXTERN

#include <complex>

using std::complex;

extern "C" {

#define ARPACK(s) s ## _

// iterate

// double precision real
extern int ARPACK(dsaupd) (int *ido, char *bmat, int *n, char *which, int *nev, double *tol,
                         double *resid, int *ncv, double *v, int *ldv,
                         int *iparam, int *ipntr, double *workd, double *workl,
                         int *lworkl, int *info);

// single precision complex
extern int ARPACK(cnaupd) (int *ido, char *bmat, int *n, char *which, int *nev, float *tol,
                         complex<float> *resid, int *ncv, complex<float> *v, int *ldv,
                         int *iparam, int *ipntr, complex<float> *workd, complex<float> *workl,
                         int *lworkl, float *rwork, int *info);

// double precision complex
extern int ARPACK(znaupd) (int *ido, char *bmat, int *n, char *which, int *nev, double *tol,
                         complex<double> *resid, int *ncv, complex<double> *v, int *ldv, 
                         int *iparam, int *ipntr, complex<double> *workd, complex<double> *workl, 
                         int *lworkl, double *rwork, int *info);

// get eigenvalues when done

// double precision real
extern int ARPACK(dseupd) (int *comp_evecs, char *howmany, int *select, double *evals, 
			 double *v, int *ldv, double *sigma,
			 char *bmat, int *n, char *which, int *nev, double *tol, double *resid, 
                         int *ncv, double *v1, int *ldv1, int *iparam, int *ipntr, 
                         double *workd, double *workl, int *lworkl, int *info);			

// single precision complex
extern int ARPACK(cneupd) (int *comp_evecs, char *howmany, int *select, complex<float> *evals, 
			 complex<float> *v, int *ldv, complex<float> *sigma, complex<float> *workev, 
			 char *bmat, int *n, char *which, int *nev, float *tol, complex<float> *resid, 
                         int *ncv, complex<float> *v1, int *ldv1, int *iparam, int *ipntr, 
                         complex<float> *workd, complex<float> *workl, int *lworkl, float *rwork, int *info);			

// double precision complex
extern int ARPACK(zneupd) (int *comp_evecs, char *howmany, int *select, complex<double> *evals, 
			 complex<double> *v, int *ldv, complex<double> *sigma, complex<double> *workev, 
			 char *bmat, int *n, char *which, int *nev, double *tol, complex<double> *resid, 
                         int *ncv, complex<double> *v1, int *ldv1, int *iparam, int *ipntr, 
                         complex<double> *workd, complex<double> *workl, int *lworkl, double *rwork, int *info);
	
}

#endif


// Based on the quda arpack interface Alexei supplied to me.

#ifndef ARPACK_EXTERN
#define ARPACK_EXTERN

#define ARPACK(s) s ## _

// iterate

// double precision real
extern int ARPACK(dsaupd) (int *ido, char *bmat, int *n, char *which, int *nev, double *tol,
                         double *resid, int *ncv, double *v, int *ldv,
                         int *iparam, int *ipntr, double *workd, double *workl,
                         int *lworkl, int *info);

// single precision complex
extern int ARPACK(cnaupd) (int *ido, char *bmat, int *n, char *which, int *nev, float *tol,
                         _Complex float *resid, int *ncv, _Complex float *v, int *ldv,
                         int *iparam, int *ipntr, _Complex float *workd, _Complex float *workl,
                         int *lworkl, float *rwork, int *info);

// double precision complex
extern int ARPACK(znaupd) (int *ido, char *bmat, int *n, char *which, int *nev, double *tol,
                         _Complex double *resid, int *ncv, _Complex double *v, int *ldv, 
                         int *iparam, int *ipntr, _Complex double *workd, _Complex double *workl, 
                         int *lworkl, double *rwork, int *info);

// get eigenvalues when done

// double precision real
extern int ARPACK(dseupd) (int *comp_evecs, char *howmany, int *select, double *evals, 
			 double *v, int *ldv, double *sigma,
			 char *bmat, int *n, char *which, int *nev, double *tol, double *resid, 
                         int *ncv, double *v1, int *ldv1, int *iparam, int *ipntr, 
                         double *workd, double *workl, int *lworkl, int *info);			

// single precision complex
extern int ARPACK(cneupd) (int *comp_evecs, char *howmany, int *select, _Complex float *evals, 
			 _Complex float *v, int *ldv, _Complex float *sigma, _Complex float *workev, 
			 char *bmat, int *n, char *which, int *nev, float *tol, _Complex float *resid, 
                         int *ncv, _Complex float *v1, int *ldv1, int *iparam, int *ipntr, 
                         _Complex float *workd, _Complex float *workl, int *lworkl, float *rwork, int *info);			

// double precision complex
extern int ARPACK(zneupd) (int *comp_evecs, char *howmany, int *select, _Complex double *evals, 
			 _Complex double *v, int *ldv, _Complex double *sigma, _Complex double *workev, 
			 char *bmat, int *n, char *which, int *nev, double *tol, _Complex double *resid, 
                         int *ncv, _Complex double *v1, int *ldv1, int *iparam, int *ipntr, 
                         _Complex double *workd, _Complex double *workl, int *lworkl, double *rwork, int *info);

#endif


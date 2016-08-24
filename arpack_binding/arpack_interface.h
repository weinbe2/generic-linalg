// ESW 2015-09-07
// A simple wrapper around arpack functionality.

// struct that gives information about the solve.
typedef struct
{
   int is_error; // is non-zero if there's an error.
   int nconv; // number of converged eigenvalues (IPARAM(5))
   int niter; // number of update iterations (IPARAM(3))
   int nops; // number of OP*x applications (IPARAM(9))
} arpack_solve_t;

/************************
* Double Real Symmetric *
************************/

// double precision real sym struct.
typedef struct
{
   // Scalar types which set memory allocation.
   int ldv; // leading dimension of the eigenvectors in V.
            // This doesn't have to be the dimension of the
            // matrix, but it can be.
            
   int maxn; // Maximum possible dimension of the matrix.
             // This doesn't have to equal N, but you might for
             // alignment.
             
   int maxnev; // Maximum number of eigenvalues you'll want.
   
   int maxncv; // Maximum number of basis vectors within IRAP.
   
   // Vectors.
   double* v; // 2d array of size [ldv][maxncv]. Holds eigenvectors.
   double* d; // 1d array of size [maxncv]. Holds eigenvalues.
   double* workl; // 1d array of size maxncv*(maxncv+8). 
                  // Workspace. Size suggested by dssimp.f.
   double* workd; // 1d array of size 3*maxn.
                  // This space gets used in reverse feedback to
                  // hold rhs = A*lhs.
   double* resid; // 1d array of size maxn.
                  // This holds internal residual vectors.       
                  // This is normally zero at the start, but it
                  // can be set to non-zero to choose an initial
                  // guess vector.
   int* select; // Should be a logical. 
   
   // Are we allocated? 0 no, 1 yes.
   int is_allocated; 
   
} arpack_drs_t;

// Initialize the arpack_drs_t struct. This largely allocates
// memory and sets internal values. 
// Arg 00: int maxn. Maximum possible dimension of the matrix.
// Arg 01: int maxnev. Maximum possible requested eigenvalues.
// Arg 02: int maxncv. Maximum possible internal values.

arpack_drs_t* arpack_drs_init(int maxn, int maxnev, int maxncv); 

// Get the requested number of eigenvalues and eigenvectors.
// Returns information about the eigensolve. 
// Arg 00: arpack_drs_t* arpack_str. [In/Out]. Initialized arpack.
// Arg 01: double* eval. [In/Out] An nev length array of eigenvalues.
// Arg 02: double* evec. [In/Out] An nev*n length array of eigenvectors.
// Arg 03: int n. Dimension of matrix.
// Arg 04: int nev. Number of requested eigenvalues.
// Arg 05: int ncv. Number of internal vectors.
// Arg 06: int maxitr. Maximum number of iterations.
// Arg 07: char* which. What part of spectrum to get.
// Arg 08: double tol. What tolerance. 0 = machine precision.
// Arg 09: double sigma. Set to 0.0 for now.
// Arg 10: void (*matrix_vector)(double*,double*,void*). Callback
//           function. Arguments are rhs, lhs, extra_info.
// Arg 11: void* extra_info passed to functions.
arpack_solve_t arpack_drs_getev(arpack_drs_t* arpack_str, double* eval, double* evec, int n, int nev, int ncv, int maxitr, char* which, double tol, double sigma, void (*matrix_vector)(double*,double*,void*), void* extra_info);

// Clean up the arpack_drs_t structure and free it. 
void arpack_drs_free(arpack_drs_t** arpack_str);

/*************************
* Double Complex General *
*************************/

// double precision complex general struct.
typedef struct
{
   // Scalar types which set memory allocation.
   int ldv; // leading dimension of the eigenvectors in V.
            // This doesn't have to be the dimension of the
            // matrix, but it can be.
            
   int maxn; // Maximum possible dimension of the matrix.
             // This doesn't have to equal N, but you might for
             // alignment.
             
   int maxnev; // Maximum number of eigenvalues you'll want.
   
   int maxncv; // Maximum number of basis vectors within IRAP.
   
   // Vectors.
   _Complex double* v; // 2d array of size [ldv][maxncv]. Holds eigenvectors.
   _Complex double* d; // 1d array of size [maxncv]. Holds eigenvalues.
   _Complex double* workl; // 1d array of size 3*maxncv**2
				  // + 5*maxncv.
                  // Workspace. Size suggested by znsimp.f.
   _Complex double* workd; // 1d array of size 3*maxn.
                  // This space gets used in reverse feedback to
                  // hold rhs = A*lhs.
   _Complex double* resid; // 1d array of size maxn.
                  // This holds internal residual vectors.       
                  // This is normally zero at the start, but it
                  // can be set to non-zero to choose an initial
                  // guess vector.
   _Complex double* workev; // Working space. 1d array of size 2*maxncv.
   double* rwork; // Working space. 1d array of size maxncv.
   double* rd; // Working space. 2d array of size [maxncv][3];
   int* select; // Should be a logical. 
   
   // Are we allocated? 0 no, 1 yes.
   int is_allocated; 
   
} arpack_dcn_t;

// Initialize the arpack_drs_t struct. This largely allocates
// memory and sets internal values. 
// Arg 00: int maxn. Maximum possible dimension of the matrix.
// Arg 01: int maxnev. Maximum possible requested eigenvalues.
// Arg 02: int maxncv. Maximum possible internal values.

arpack_dcn_t* arpack_dcn_init(int maxn, int maxnev, int maxncv); 

// Get the requested number of eigenvalues and eigenvectors.
// Returns information about the eigensolve. 
// Arg 00: arpack_dcn_t* arpack_str. [In/Out]. Initialized arpack.
// Arg 01: _Complex double* eval. [In/Out] An nev length array of eigenvalues.
// Arg 02: _Complex double* evec. [In/Out] An nev*n length array of eigenvectors.
// Arg 03: int n. Dimension of matrix.
// Arg 04: int nev. Number of requested eigenvalues.
// Arg 05: int ncv. Number of internal vectors.
// Arg 06: int maxitr. Maximum number of iterations.
// Arg 07: char* which. What part of spectrum to get.
// Arg 08: double tol. What tolerance. 0 = machine precision.
// Arg 09: _Complex double sigma. Set to 0.0 for now.
// Arg 10: void (*matrix_vector)(double*,double*,void*). Callback
//           function. Arguments are rhs, lhs, extra_info.
// Arg 11: void* extra_info passed to functions.
arpack_solve_t arpack_dcn_getev(arpack_dcn_t* arpack_str, _Complex double* eval, _Complex double* evec, int n, int nev, int ncv, int maxitr, char* which, double tol, _Complex double sigma, void (*matrix_vector)(_Complex double*,_Complex double*,void*), void* extra_info);

// Get the requested number of eigenvalues and eigenvectors.
// Returns information about the eigensolve. 
// Uses shift-invert.
// Arg 00: arpack_dcn_t* arpack_str. [In/Out]. Initialized arpack.
// Arg 01: _Complex double* eval. [In/Out] An nev length array of eigenvalues.
// Arg 02: _Complex double* evec. [In/Out] An nev*n length array of eigenvectors.
// Arg 03: int n. Dimension of matrix.
// Arg 04: int nev. Number of requested eigenvalues.
// Arg 05: int ncv. Number of internal vectors.
// Arg 06: int maxitr. Maximum number of iterations.
// Arg 07: char* which. What part of spectrum to get.
// Arg 08: double ev tol. What tolerance. 0 = machine precision.
// Arg 09: _Complex double sigma. Set to 0.0 for now.
// Arg 10: void (*matrix_vector)(double*,double*,void*). Callback
//           function. Arguments are rhs, lhs, extra_info.
// Arg 11: void (*minv_vector)(double*,double*,int,int,double,matrix_vector_p,double,void*). Callback function for matrix inversion.
// Arg 12: int maxitr_cg. Number of iterations for cg.
// Arg 13: double tol_cg. CG precision.
// Arg 14: void* extra_info passed to functions.
typedef void(*matrix_vector_p)(_Complex double*,_Complex double*,void*); 
arpack_solve_t arpack_dcn_getev_sinv(arpack_dcn_t* arpack_str, _Complex double* eval, _Complex double* evec, int n, int nev, int ncv, int maxitr, char* which, double tol, _Complex double sigma, void (*matrix_vector)(_Complex double*,_Complex double*,void*), double (*minv_vector)(_Complex double*,_Complex double*,int,int,double,matrix_vector_p,_Complex double,void*), int maxitr_cg, double tol_cg ,void* extra_info);

// Clean up the arpack_drs_t structure and free it. 
void arpack_dcn_free(arpack_dcn_t** arpack_str);

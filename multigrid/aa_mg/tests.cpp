// Various test routines. Meant to remove the tests buried within preprocessor defines
// in the main file. 

#ifndef MG_TEST_FILE
#define MG_TEST_FILE

#include <complex>
#include <cstring>
using namespace std;

#include "mg.h"
#include "mg_complex.h"
#include "generic_vector.h"
#include "coarse_stencil.h"
#include "tests.h"
#include "operators.h"

#include "generic_gcr.h"
#include "arpack_interface.h"

// Test eigenvalue overlaps. 
void test_eigenvalue_overlap(mg_operator_struct_complex* mgstruct_ptr, staggered_u1_op* stagif, op_type opt, int set_eigen, int set_cv)
{
    mg_operator_struct_complex mgstruct = *mgstruct_ptr; 
    
    int i, j, k;
    inversion_info invif; 

    complex<double>** evals = new complex<double>*[mgstruct.n_refine+1];
    complex<double>*** evecs = new complex<double>**[mgstruct.n_refine+1];
    int* n_eigen = new int[mgstruct.n_refine];
    int lev = 0; 
    int n_cv = 0; 

    for (lev = 0; lev < mgstruct.n_refine; lev++)
    {
        n_eigen[lev] = 0;
        n_cv = 0;

        if (set_eigen == -1 && set_cv == -1) // generate all eigenvalues, eigenvectors. 
        {
            // Allocate space for all eigenvalues, eigenvectors. 
            n_eigen[lev] = mgstruct.curr_fine_size;
            n_cv = mgstruct.curr_fine_size; 
            cout << "[L" << lev+1 << "_ARPACK]: Number of eigenvalues: " << n_eigen[lev] << " Number of cv: " << n_cv << "\n";

            evals[lev] = new complex<double>[mgstruct.curr_fine_size];
            evecs[lev] = new complex<double>*[mgstruct.curr_fine_size];
            for (i = 0; i < n_eigen[lev]; i++)
            {
                evecs[lev][i] = new complex<double>[mgstruct.curr_fine_size];
            }

            // Get low mag half
            arpack_dcn_t* ar_strc = arpack_dcn_init(mgstruct.curr_fine_size, n_eigen[lev]/2, n_cv); // max eigenvectors, internal vecs
            char eigtype[3]; strcpy(eigtype, "SM"); // Smallest magnitude eigenvalues.
            arpack_solve_t info_solve = arpack_dcn_getev(ar_strc, evals[lev], evecs[lev], mgstruct.curr_fine_size, n_eigen[lev]/2, n_cv, 4000, eigtype, 1e-7, 0.0, fine_square_staggered, (void*)&mgstruct); 
            //arpack_dcn_free(&ar_strc);

            // Print info about the eigensolve.
            cout << "[L" << lev+1 << "_ARPACK]: Number of converged eigenvalues: " << info_solve.nconv << "\n";
            cout << "[L" << lev+1 << "_ARPACK]: Number of iteration steps: " << info_solve.niter << "\n";
            cout << "[L" << lev+1 << "_ARPACK]: Number of matrix multiplies: " << info_solve.nops << "\n";

            // Get high mag half
            //arpack_dcn_t* ar_strc = arpack_dcn_init(mgstruct.curr_fine_size, n_eigen, n_cv); // max eigenvectors, internal vecs
            strcpy(eigtype, "LM"); // Smallest magnitude eigenvalues.
            info_solve = arpack_dcn_getev(ar_strc, evals[lev]+(mgstruct.curr_fine_size/2), evecs[lev]+(mgstruct.curr_fine_size/2), mgstruct.curr_fine_size, n_eigen[lev]/2, n_cv, 4000, eigtype, 1e-7, 0.0, fine_square_staggered, (void*)&mgstruct); 
            arpack_dcn_free(&ar_strc);

            // Print info about the eigensolve.
            cout << "[L" << lev+1 << "_ARPACK]: Number of converged eigenvalues: " << info_solve.nconv << "\n";
            cout << "[L" << lev+1 << "_ARPACK]: Number of iteration steps: " << info_solve.niter << "\n";
            cout << "[L" << lev+1 << "_ARPACK]: Number of matrix multiplies: " << info_solve.nops << "\n";

            // End of arpack bindings!   
        }
        else if (set_eigen != -1 && set_cv == -1) // generate n_eigen eigenvalues, min(mgstruct.curr_fine_size, 2.5 n_eigen) cv.
        {
            n_eigen[lev] = set_eigen;
            n_cv = min(mgstruct.curr_fine_size, 2*n_eigen[lev] + n_eigen[lev]/2);
            cout << "[L" << lev+1 << "_ARPACK]: Number of eigenvalues: " << n_eigen[lev] << " Number of cv: " << n_cv << "\n";

            evals[lev] = new complex<double>[n_eigen[lev]];
            evecs[lev] = new complex<double>*[n_eigen[lev]];
            for (i = 0; i < n_eigen[lev]; i++)
            {
                evecs[lev][i] = new complex<double>[mgstruct.curr_fine_size];
            }

            // Get low mag half
            arpack_dcn_t* ar_strc = arpack_dcn_init(mgstruct.curr_fine_size, n_eigen[lev], n_cv); // max eigenvectors, internal vecs
            char eigtype[3]; strcpy(eigtype, "SM"); // Smallest magnitude eigenvalues.
            arpack_solve_t info_solve = arpack_dcn_getev(ar_strc, evals[lev], evecs[lev], mgstruct.curr_fine_size, n_eigen[lev], n_cv, 4000, eigtype, 1e-7, 0.0, fine_square_staggered, (void*)&mgstruct); 
            arpack_dcn_free(&ar_strc);

            // Print info about the eigensolve.
            cout << "[L" << lev+1 << "_ARPACK]: Number of converged eigenvalues: " << info_solve.nconv << "\n";
            cout << "[L" << lev+1 << "_ARPACK]: Number of iteration steps: " << info_solve.niter << "\n";
            cout << "[L" << lev+1 << "_ARPACK]: Number of matrix multiplies: " << info_solve.nops << "\n";

        }
        else // generate n_eigen eigenvalues, min(mgstruct.curr_fine_size, n_cv) cv.
        {
            n_eigen[lev] = set_eigen;
            n_cv = set_cv;

            cout << "[L" << lev+1 << "_ARPACK]: Number of eigenvalues: " << n_eigen[lev] << " Number of cv: " << n_cv << "\n";

            evals[lev] = new complex<double>[n_eigen[lev]];
            evecs[lev] = new complex<double>*[n_eigen[lev]];
            for (i = 0; i < n_eigen[lev]; i++)
            {
                evecs[lev][i] = new complex<double>[mgstruct.curr_fine_size];
            }

            // Get low mag half
            arpack_dcn_t* ar_strc = arpack_dcn_init(mgstruct.curr_fine_size, n_eigen[lev], n_cv); // max eigenvectors, internal vecs
            char eigtype[3]; strcpy(eigtype, "SM"); // Smallest magnitude eigenvalues.
            arpack_solve_t info_solve = arpack_dcn_getev(ar_strc, evals[lev], evecs[lev], mgstruct.curr_fine_size, n_eigen[lev], n_cv, 4000, eigtype, 1e-7, 0.0, fine_square_staggered, (void*)&mgstruct); 
            arpack_dcn_free(&ar_strc);

            // Print info about the eigensolve.
            cout << "[L" << lev+1 << "_ARPACK]: Number of converged eigenvalues: " << info_solve.nconv << "\n";
            cout << "[L" << lev+1 << "_ARPACK]: Number of iteration steps: " << info_solve.niter << "\n";
            cout << "[L" << lev+1 << "_ARPACK]: Number of matrix multiplies: " << info_solve.nops << "\n";
        }

        // Sort eigenvalues (done differently depending on the operator).
        for (i = 0; i < n_eigen[lev]; i++)
        {
            for (j = 0; j < n_eigen[lev]-1; j++)
            {
                switch (opt)
                {
                    case STAGGERED:
                        if (abs(imag(evals[lev][j])) > abs(imag(evals[lev][j+1])))
                        {
                            complex<double> teval = evals[lev][j]; evals[lev][j] = evals[lev][j+1]; evals[lev][j+1] = teval;
                            complex<double>* tevec = evecs[lev][j]; evecs[lev][j] = evecs[lev][j+1]; evecs[lev][j+1] = tevec;
                        }
                        break;
                    case LAPLACE:
                    case LAPLACE_NC2:
                    case G5_STAGGERED:
                    case STAGGERED_NORMAL:
                    case STAGGERED_INDEX:
                        if (abs(real(evals[lev][j])) > abs(real(evals[lev][j+1])))
                        {
                            complex<double> teval = evals[lev][j]; evals[lev][j] = evals[lev][j+1]; evals[lev][j+1] = teval;
                            complex<double>* tevec = evecs[lev][j]; evecs[lev][j] = evecs[lev][j+1]; evecs[lev][j+1] = tevec;
                        }
                        break; 
                }
            }
        }


        cout << "\n\nAll eigenvalues:\n";
        for (i = 0; i < n_eigen[lev]; i++)
        {
            cout << "[L" << lev+1 << "_FINEVAL]: Mass " << stagif->mass << " Num " << i << " Eval " << evals[lev][i] << "\n";
            normalize<double>(evecs[lev][i], mgstruct.curr_fine_size);
        }

        if (lev < mgstruct.n_refine-1)
        {
            level_down(&mgstruct);
        }
    }

    // Something special for the coarsest level
    lev = mgstruct.n_refine; 
    n_eigen[lev] = 0;
    n_cv = 0;

    if (set_eigen == -1 && set_cv == -1) // generate all eigenvalues, eigenvectors. 
    {
        // Allocate space for all eigenvalues, eigenvectors. 
        n_eigen[lev] = mgstruct.curr_coarse_size;
        n_cv = mgstruct.curr_coarse_size; 
        cout << "[L" << lev+1 << "_ARPACK]: Number of eigenvalues: " << n_eigen[lev] << " Number of cv: " << n_cv << "\n";

        evals[lev] = new complex<double>[mgstruct.curr_coarse_size];
        evecs[lev] = new complex<double>*[mgstruct.curr_coarse_size];
        for (i = 0; i < n_eigen[lev]; i++)
        {
            evecs[lev][i] = new complex<double>[mgstruct.curr_coarse_size];
        }

        // Get low mag half
        arpack_dcn_t* ar_strc = arpack_dcn_init(mgstruct.curr_coarse_size, n_eigen[lev]/2, n_cv); // max eigenvectors, internal vecs
        char eigtype[3]; strcpy(eigtype, "SM"); // Smallest magnitude eigenvalues.
        arpack_solve_t info_solve = arpack_dcn_getev(ar_strc, evals[lev], evecs[lev], mgstruct.curr_coarse_size, n_eigen[lev]/2, n_cv, 4000, eigtype, 1e-7, 0.0, coarse_square_staggered, (void*)&mgstruct); 
        //arpack_dcn_free(&ar_strc);

        // Print info about the eigensolve.
        cout << "[L" << lev+1 << "_ARPACK]: Number of converged eigenvalues: " << info_solve.nconv << "\n";
        cout << "[L" << lev+1 << "_ARPACK]: Number of iteration steps: " << info_solve.niter << "\n";
        cout << "[L" << lev+1 << "_ARPACK]: Number of matrix multiplies: " << info_solve.nops << "\n";

        // Get high mag half
        //arpack_dcn_t* ar_strc = arpack_dcn_init(mgstruct.curr_coarse_size, n_eigen, n_cv); // max eigenvectors, internal vecs
        strcpy(eigtype, "LM"); // Smallest magnitude eigenvalues.
        info_solve = arpack_dcn_getev(ar_strc, evals[lev]+(mgstruct.curr_coarse_size/2), evecs[lev]+(mgstruct.curr_coarse_size/2), mgstruct.curr_coarse_size, n_eigen[lev]/2, n_cv, 4000, eigtype, 1e-7, 0.0, coarse_square_staggered, (void*)&mgstruct); 
        arpack_dcn_free(&ar_strc);

        // Print info about the eigensolve.
        cout << "[L" << lev+1 << "_ARPACK]: Number of converged eigenvalues: " << info_solve.nconv << "\n";
        cout << "[L" << lev+1 << "_ARPACK]: Number of iteration steps: " << info_solve.niter << "\n";
        cout << "[L" << lev+1 << "_ARPACK]: Number of matrix multiplies: " << info_solve.nops << "\n";

        // End of arpack bindings!   
    }
    else if (set_eigen != -1 && set_cv == -1) // generate n_eigen eigenvalues, min(mgstruct.curr_coarse_size, 2.5 n_eigen) cv.
    {
        n_eigen[lev] = set_eigen;
        n_cv = min(mgstruct.curr_coarse_size, 2*n_eigen[lev] + n_eigen[lev]/2);
        cout << "[L" << lev+1 << "_ARPACK]: Number of eigenvalues: " << n_eigen[lev] << " Number of cv: " << n_cv << "\n";

        evals[lev] = new complex<double>[n_eigen[lev]];
        evecs[lev] = new complex<double>*[n_eigen[lev]];
        for (i = 0; i < n_eigen[lev]; i++)
        {
            evecs[lev][i] = new complex<double>[mgstruct.curr_coarse_size];
        }

        // Get low mag half
        arpack_dcn_t* ar_strc = arpack_dcn_init(mgstruct.curr_coarse_size, n_eigen[lev], n_cv); // max eigenvectors, internal vecs
        char eigtype[3]; strcpy(eigtype, "SM"); // Smallest magnitude eigenvalues.
        arpack_solve_t info_solve = arpack_dcn_getev(ar_strc, evals[lev], evecs[lev], mgstruct.curr_coarse_size, n_eigen[lev], n_cv, 4000, eigtype, 1e-7, 0.0, coarse_square_staggered, (void*)&mgstruct); 
        arpack_dcn_free(&ar_strc);

        // Print info about the eigensolve.
        cout << "[L" << lev+1 << "_ARPACK]: Number of converged eigenvalues: " << info_solve.nconv << "\n";
        cout << "[L" << lev+1 << "_ARPACK]: Number of iteration steps: " << info_solve.niter << "\n";
        cout << "[L" << lev+1 << "_ARPACK]: Number of matrix multiplies: " << info_solve.nops << "\n";

    }
    else // generate n_eigen eigenvalues, min(mgstruct.curr_coarse_size, n_cv) cv.
    {
        n_eigen[lev] = set_eigen;
        n_cv = set_cv;

        cout << "[L" << lev+1 << "_ARPACK]: Number of eigenvalues: " << n_eigen[lev] << " Number of cv: " << n_cv << "\n";

        evals[lev] = new complex<double>[n_eigen[lev]];
        evecs[lev] = new complex<double>*[n_eigen[lev]];
        for (i = 0; i < n_eigen[lev]; i++)
        {
            evecs[lev][i] = new complex<double>[mgstruct.curr_coarse_size];
        }

        // Get low mag half
        arpack_dcn_t* ar_strc = arpack_dcn_init(mgstruct.curr_coarse_size, n_eigen[lev], n_cv); // max eigenvectors, internal vecs
        char eigtype[3]; strcpy(eigtype, "SM"); // Smallest magnitude eigenvalues.
        arpack_solve_t info_solve = arpack_dcn_getev(ar_strc, evals[lev], evecs[lev], mgstruct.curr_coarse_size, n_eigen[lev], n_cv, 4000, eigtype, 1e-7, 0.0, coarse_square_staggered, (void*)&mgstruct); 
        arpack_dcn_free(&ar_strc);

        // Print info about the eigensolve.
        cout << "[L" << lev+1 << "_ARPACK]: Number of converged eigenvalues: " << info_solve.nconv << "\n";
        cout << "[L" << lev+1 << "_ARPACK]: Number of iteration steps: " << info_solve.niter << "\n";
        cout << "[L" << lev+1 << "_ARPACK]: Number of matrix multiplies: " << info_solve.nops << "\n";
    }

    // Sort eigenvalues (done differently depending on the operator).
    for (i = 0; i < n_eigen[lev]; i++)
    {
        for (j = 0; j < n_eigen[lev]-1; j++)
        {
            switch (opt)
            {
                case STAGGERED:
                    if (abs(imag(evals[lev][j])) > abs(imag(evals[lev][j+1])))
                    {
                        complex<double> teval = evals[lev][j]; evals[lev][j] = evals[lev][j+1]; evals[lev][j+1] = teval;
                        complex<double>* tevec = evecs[lev][j]; evecs[lev][j] = evecs[lev][j+1]; evecs[lev][j+1] = tevec;
                    }
                    break;
                case LAPLACE:
                case LAPLACE_NC2:
                case G5_STAGGERED:
                case STAGGERED_NORMAL:
                case STAGGERED_INDEX:
                    if (abs(real(evals[lev][j])) > abs(real(evals[lev][j+1])))
                    {
                        complex<double> teval = evals[lev][j]; evals[lev][j] = evals[lev][j+1]; evals[lev][j+1] = teval;
                        complex<double>* tevec = evecs[lev][j]; evecs[lev][j] = evecs[lev][j+1]; evecs[lev][j+1] = tevec;
                    }
                    break; 
            }
        }
    }


    cout << "\n\nAll eigenvalues:\n";
    for (i = 0; i < n_eigen[lev]; i++)
    {
        cout << "[L" << lev+1 << "_FINEVAL]: Mass " << stagif->mass << " Num " << i << " Eval " << evals[lev][i] << "\n";
        normalize<double>(evecs[lev][i], mgstruct.curr_coarse_size);
    }

    // End generating coarse level. 



    for (lev = mgstruct.n_refine-2; lev >= 0; lev--)
    {
        level_up(&mgstruct);
    }

    for (lev = 0; lev < mgstruct.n_refine; lev++)
    {

        complex<double>* evec_Pdag = new complex<double>[mgstruct.curr_coarse_size];
        complex<double>* evec_Pdag2 = new complex<double>[mgstruct.curr_coarse_size];
        complex<double>* evec_PPdag = new complex<double>[mgstruct.curr_fine_size];

        // Test overlap of null vectors with eigenvectors.
        // Formally, this is looking at the magnitude of (1 - P P^\dag) eigenvector.

        for (i = 0; i < n_eigen[lev]; i++)
        {
            // Zero out.
            zero<double>(evec_Pdag, mgstruct.curr_coarse_size);
            zero<double>(evec_PPdag, mgstruct.curr_fine_size);

            // Restrict eigenvector.
            restrict(evec_Pdag, evecs[lev][i], &mgstruct);

            // Prolong.
            prolong(evec_PPdag, evec_Pdag, &mgstruct);

            // Subtract off eigenvector, take norm.
            for (j = 0; j < mgstruct.curr_fine_size; j++)
            {
                evec_PPdag[j] -= evecs[lev][i][j];
            }

            cout << "[L" << lev+1 << "_1mPPDAG]: Num " << i << " Overlap " << sqrt(norm2sq<double>(evec_PPdag, mgstruct.curr_fine_size)) << "\n"; 
        }

        // Test how good of a preconditioner the coarse operator is.
        // Formally, this is looking at the magnitude of (1 - P ( P^\dag A P )^(-1) P^\dag A) eigenvector.
        for (i = 0; i < n_eigen[lev]; i++)
        {
            // Zero out.
            zero<double>(evec_Pdag, mgstruct.curr_coarse_size);
            zero<double>(evec_Pdag2, mgstruct.curr_coarse_size);
            zero<double>(evec_PPdag, mgstruct.curr_fine_size);

            // Apply A.
            fine_square_staggered(evec_PPdag, evecs[lev][i], (void*)&mgstruct);

            // Restrict.
            restrict(evec_Pdag, evec_PPdag, &mgstruct);

            // Try a deflation preconditioned solve...
            for (j = 0; j < n_eigen[lev+1]; j++)
            {
                complex<double> def_dot = dot<double>(evecs[lev+1][j], evec_Pdag, mgstruct.curr_coarse_size);
                for (k = 0; k < mgstruct.curr_coarse_size; k++)
                {
                    evec_Pdag2[k] += 1.0/(evals[lev+1][j])*def_dot*evecs[lev+1][j][k];
                }
            }

            // Invert A_coarse against it.
            invif = minv_vector_gcr_restart(evec_Pdag2, evec_Pdag, mgstruct.curr_coarse_size, 10000, 1e-7, 64, coarse_square_staggered, (void*)&mgstruct);
            //cout << "[L" << lev+1 << "_DEFLATE]: Num " << i << " Iter " << invif.iter << "\n";

            // Prolong.
            zero<double>(evec_PPdag, mgstruct.curr_coarse_size);
            prolong(evec_PPdag, evec_Pdag2, &mgstruct);

            // Subtract off eigenvector, take norm.
            for (j = 0; j < mgstruct.curr_fine_size; j++)
            {
                evec_PPdag[j] -= evecs[lev][i][j];
            }

            cout << "[L" << lev+1 << "_1mP_Ac_PDAG_A]: Num " << i << " Overlap " << sqrt(norm2sq<double>(evec_PPdag, mgstruct.curr_fine_size)) << "\n";
        }


        delete[] evec_Pdag;
        delete[] evec_PPdag;
        delete[] evec_Pdag2;

        if (lev < mgstruct.n_refine-1)
        {
            level_down(&mgstruct);
        }
    }

    for (lev = mgstruct.n_refine-2; lev >= 0; lev--)
    {
        level_up(&mgstruct);
    }

    for (lev = 0; lev <= mgstruct.n_refine; lev++)
    {
        for (i = 0; i < n_eigen[lev]; i++)
        {
            delete[] evecs[lev][i];
        }
        delete[] evals[lev];
        delete[] evecs[lev];
    }

    delete[] evecs;
    delete[] evals;
    delete[] n_eigen;


}

// Test constructing the coarse operator, comparing to using prolong/restrict of fine.
void test_stencil_construct(mg_operator_struct_complex* mgstruct, int level, int stencil_size)
{
    // A few sanity checks...
    if (level > mgstruct->n_refine+1)
    {
        cout << "[TEST_ERROR]: Stencil construct level " << level << " is larger than the number of mg levels.\n" << flush;
        return;
    }
    
    if (stencil_size > 2)
    {
        cout << "[TEST_ERROR]: Stencil distances greater than 2 are not supported.\n" << flush; 
        return;
    }
    
    stencil_2d stenc(mgstruct->latt[level], stencil_size);
    
    // Save mgstruct level state.
    int save_level = mgstruct->curr_level;
    if (save_level > 0)
    {
        for (int i = 0; i < save_level; i++)
        {
            level_up(mgstruct);
        }
    }
    
    // Push down an appropriate number of levels.
    if (level > 0)
    {
        for (int i = 0; i < level-1; i++)
        {
            level_down(mgstruct);
        }
    }
    
    // Carry on. 
    
    
    if (level == 0)
    {
        generate_stencil_2d(&stenc, fine_square_staggered, (void*)mgstruct);
    }
    else
    {
        generate_stencil_2d(&stenc, coarse_square_staggered, (void*)mgstruct);
    }
    
    cout << "[TEST]: Level " << level+1 << " Generated stencil.\n" << flush; 
    
    complex<double>* tmp_rhs = new complex<double>[mgstruct->latt[level]->get_lattice_size()];
    complex<double>* tmp_lhs = new complex<double>[mgstruct->latt[level]->get_lattice_size()];
    
    // Whelp, it's something. Let's test it.
    //zero<double>(tmp_rhs, mgstruct->latt[level]->get_lattice_size());
    //tmp_rhs[0] = 1.0;
    std::mt19937 generator (1339u); // RNG, 1339u is the seed. 
    gaussian<double>(tmp_rhs, mgstruct->latt[level]->get_lattice_size(), generator); 
    apply_stencil_2d(tmp_lhs, tmp_rhs, (void*)&stenc);
    
    complex<double>* tmp_lhs2 = new complex<double>[mgstruct->latt[level]->get_lattice_size()];
    if (level == 0)
    {
        fine_square_staggered(tmp_lhs2, tmp_rhs, (void*)mgstruct);
    }
    else
    {
        coarse_square_staggered(tmp_lhs2, tmp_rhs, (void*)mgstruct);
    }
    
    // Get squared difference.
    cout << "[TEST]: Level " << level+1 << " Squared difference: " << diffnorm2sq<double>(tmp_lhs, tmp_lhs2, mgstruct->latt[level]->get_lattice_size()) << "\n" << flush; 
    
    delete[] tmp_rhs; 
    delete[] tmp_lhs;
    delete[] tmp_lhs2;
    
    // Restore level
    if (level > 0)
    {
        for (int i = 0; i < level-1; i++)
        {
            level_up(mgstruct);
        }
    }
    
    // Restore level.
    if (save_level > 0)
    {
        for (int i = 0; i < save_level; i++)
        {
            level_down(mgstruct);
        }
    }
    
}


// Test constructing the coarse operator, comparing applying the full coarse stencil to the piece-by-piece stencil.
void test_stencil_piece(mg_operator_struct_complex* mgstruct, int level, int stencil_size)
{
    // A few sanity checks...
    if (level > mgstruct->n_refine+1)
    {
        cout << "[TEST_ERROR]: Stencil construct level " << level << " is larger than the number of mg levels.\n" << flush;
        return;
    }
    
    if (stencil_size > 2)
    {
        cout << "[TEST_ERROR]: Stencil distances greater than 2 are not supported.\n" << flush; 
        return;
    }
    
    stencil_2d stenc(mgstruct->latt[level], stencil_size);
    
    // Save mgstruct level state.
    int save_level = mgstruct->curr_level;
    if (save_level > 0)
    {
        for (int i = 0; i < save_level; i++)
        {
            level_up(mgstruct);
        }
    }
    
    // Push down an appropriate number of levels.
    if (level > 0)
    {
        for (int i = 0; i < level-1; i++)
        {
            level_down(mgstruct);
        }
    }
    
    // Carry on. 
    
    
    if (level == 0)
    {
        generate_stencil_2d(&stenc, fine_square_staggered, (void*)mgstruct);
    }
    else
    {
        generate_stencil_2d(&stenc, coarse_square_staggered, (void*)mgstruct);
    }
    
    cout << "[TEST]: Level " << level+1 << " Generated stencil.\n" << flush; 
    
    complex<double>* tmp_rhs = new complex<double>[mgstruct->latt[level]->get_lattice_size()];
    complex<double>* tmp_lhs = new complex<double>[mgstruct->latt[level]->get_lattice_size()];
    complex<double>* tmp_lhs2 = new complex<double>[mgstruct->latt[level]->get_lattice_size()];
    complex<double>* tmp_save = new complex<double>[mgstruct->latt[level]->get_lattice_size()];
    
    // Whelp, it's something. Let's test it. Use a random rhs.
    // This applies the full operator.
    stenc.sdir = DIR_ALL; 
    std::mt19937 generator (1339u); // RNG, 1339u is the seed. 
    gaussian<double>(tmp_rhs, mgstruct->latt[level]->get_lattice_size(), generator); 
    apply_stencil_2d(tmp_lhs, tmp_rhs, (void*)&stenc);
    
    // Good! Now try applying all of the stencil pieces, one at a time.
    zero<double>(tmp_save, mgstruct->latt[level]->get_lattice_size());
    // A little enum abuse... Loop over all pieces.
    for (int i = 1; i <= (int)(stencil_size == 1 ? DIR_YM1 : DIR_XP1YM1); i++)
    {
        stenc.sdir = (stencil_dir)i; // Set a direction.
        apply_stencil_2d(tmp_lhs2, tmp_rhs, (void*)&stenc); // Apply that piece of the stencil.
        for (int j = 0; j < mgstruct->latt[level]->get_lattice_size(); j++)
        {
            tmp_save[j] += tmp_lhs2[j];
        }
    }
    
    // And compare via squared difference!
    cout << "[TEST]: Level " << level+1 << " Squared difference: " << diffnorm2sq<double>(tmp_lhs, tmp_save, mgstruct->latt[level]->get_lattice_size()) << "\n" << flush; 
    
    
    delete[] tmp_rhs; 
    delete[] tmp_lhs;
    delete[] tmp_lhs2;
    delete[] tmp_save; 
    
    
    // Restore level
    if (level > 0)
    {
        for (int i = 0; i < level-1; i++)
        {
            level_up(mgstruct);
        }
    }
    
    // Restore level.
    if (save_level > 0)
    {
        for (int i = 0; i < save_level; i++)
        {
            level_down(mgstruct);
        }
    }
}



#endif // MG_TEST_FILE
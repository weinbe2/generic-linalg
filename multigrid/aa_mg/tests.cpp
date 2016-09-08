// Various test routines. Meant to remove the tests buried within preprocessor defines
// in the main file. 

#ifndef MG_TEST_FILE
#define MG_TEST_FILE

#include <complex>
using namespace std;

#include "mg.h"
#include "mg_complex.h"
#include "generic_vector.h"
#include "coarse_stencil.h"
#include "tests.h"

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
    
    stencil_2d stenc;
    stenc.lat = mgstruct->latt[level];
    
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
        generate_stencil_2d(&stenc, stencil_size, fine_square_staggered, (void*)mgstruct);
    }
    else
    {
        generate_stencil_2d(&stenc, stencil_size, coarse_square_staggered, (void*)mgstruct);
    }
    
    cout << "[TEST]: Level " << level+1 << " Generated stencil.\n" << flush; 
    
    complex<double>* tmp_rhs = new complex<double>[mgstruct->latt[level]->get_lattice_size()];
    complex<double>* tmp_lhs = new complex<double>[mgstruct->latt[level]->get_lattice_size()];
    
    // Whelp, it's something. Let's test it.
    zero<double>(tmp_rhs, mgstruct->latt[level]->get_lattice_size());
    tmp_rhs[0] = 1.0;
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
    
    delete[] stenc.clover;
    delete[] stenc.hopping;
    if (stenc.two_link)
    {
        delete[] stenc.two_link;
    }
    
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
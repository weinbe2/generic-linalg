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
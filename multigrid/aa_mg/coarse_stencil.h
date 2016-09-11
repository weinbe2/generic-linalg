// Structures and routines for the coarse stencil.

#ifndef MG_COARSE_STENCIL
#define MG_COARSE_STENCIL

#include <complex>
#include "lattice.h"

using namespace std;

// Enum for all possible stencil directions.
// Largely used for generating additional refinements.
enum stencil_dir
{
    DIR_ALL = 0, // default, full stencil.
    DIR_0 = 1,   // clover
    DIR_XP1 = 2, // +x
    DIR_YP1 = 3, // +y
    DIR_XM1 = 4, // -x
    DIR_YM1 = 5, // -y
    DIR_XP2 = 6, // +2x
    DIR_XP1YP1 = 7, // +x+y
    DIR_YP2 = 8, // +2y
    DIR_XM1YP1 = 9, // -x+y
    DIR_XM2 = 10, // -2x
    DIR_XM1YM1 = 11, // -x-y
    DIR_YM2 = 12, // -2y
    DIR_XP1YM1 = 13, // +x-y
};

struct stencil_2d
{
    // Associated lattice!
    Lattice* lat; 
    
    // What version of the stencil are we using? Default all.
    stencil_dir sdir; 
    
    // Nc x Nc x X x Y
    complex<double>* clover;
    
    // Nc x Nc x X x Y x {+X,+Y,-X,-Y} 
    complex<double>* hopping;
    
    // is there a two-link term? (for normal eqn)
    bool has_two; 
    
    // Nc x Nc x X x Y x {+2X, +X+Y, +2Y, -X+Y, -2X, -X-Y, -2Y, +X-Y}
    complex<double>* two_link; 
    
    // Have we generated a lattice?
    bool generated;
    
    // Base constructor.
    stencil_2d()
    {
        generated = false; 
        lat = 0;
        has_two = false;
        
        clover = 0;
        hopping = 0;
        two_link = 0; 
        
        sdir = DIR_ALL; 
    }
    
    ~stencil_2d()
    {
        /*if (generated)
        {
            if (clover != 0) { delete[] clover; clover = 0; }
            if (hopping != 0) { delete[] hopping; hopping = 0; }
            if (has_two)
            {
                if (two_link != 0) { delete[] two_link; two_link = 0; }
            }
            generated = false; 
        }*/
    }
    
};

// Applies an (up to two link) stencil as defined by a stencil_2d object.
void apply_stencil_2d(complex<double>* lhs, complex<double>* rhs, void* extra_data);

// Generate a stencil operator given a function and a coarsening distance (max 2...)
void generate_stencil_2d(stencil_2d* stenc, int stencil_distance, void (*matrix_vector)(complex<double>*,complex<double>*,void*), void* extra_data);

#endif // MG_COARSE_STENCIL
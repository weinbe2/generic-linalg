// Structures and routines for the coarse stencil.

#ifndef MG_COARSE_STENCIL
#define MG_COARSE_STENCIL

#include <complex>
#include "lattice.h"

using namespace std;

struct stencil_2d
{
    // Associated lattice!
    Lattice* lat; 
    
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
    }
    
    ~stencil_2d()
    {
        if (generated)
        {
            delete[] clover;
            delete[] hopping;
            if (has_two)
            {
                delete[] two_link;
            }
        }
    }
    
};

// Applies an (up to two link) stencil as defined by a stencil_2d object.
void apply_stencil_2d(complex<double>* lhs, complex<double>* rhs, void* extra_data);

// Generate a stencil operator given a function and a coarsening distance (max 2...)
void generate_stencil_2d(stencil_2d* stenc, int stencil_distance, void (*matrix_vector)(complex<double>*,complex<double>*,void*), void* extra_data);

#endif // MG_COARSE_STENCIL
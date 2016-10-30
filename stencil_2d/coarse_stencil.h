// Structures and routines for the coarse stencil.

#ifndef COARSE_STENCIL
#define COARSE_STENCIL

#include <iostream>
#include <complex>
#include "lattice.h"
#include "generic_vector.h"

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
    
    // What is the stencil size?
    int stencil_size; 
    
    // Nc x Nc x X x Y x {+2X, +X+Y, +2Y, -X+Y, -2X, -X-Y, -2Y, +X-Y}
    complex<double>* two_link; 
    
    // Have we generated a lattice?
    bool generated;
    
    // Identity shift, think local mass term. 
    complex<double> shift; 
    
    // Even/odd shift (minus sign on odd), think mass term in g5_staggered
    complex<double> eo_shift; 
    
    // Top/bottom half shift (minus sign on bottom half of color dof, think mass term on coarse g5_staggered.
    // i.e., top half of dof gets a lhs[i] += dof_shift*rhs[i], bottom half gets a rhs[i] -= dof_shift*rhs[i]
    complex<double> dof_shift;
    
    // Base constructor.
    stencil_2d(Lattice* in_lat, int in_stencil_size, complex<double> in_shift = 0.0, complex<double> in_eo_shift = 0.0, complex<double> in_dof_shift = 0.0)
        : lat(in_lat), stencil_size(in_stencil_size), shift(in_shift), eo_shift(in_eo_shift), dof_shift(in_dof_shift)
    {
        generated = false; 
            
        if (stencil_size == 2)
        {
            has_two = true;
        }
        else
        {
            has_two = false;
        }
        
        clover = new complex<double>[lat->get_volume()*lat->get_nc()*lat->get_nc()];
        zero<double>(clover, lat->get_volume()*lat->get_nc()*lat->get_nc());
        hopping = new complex<double>[lat->get_volume()*lat->get_nc()*lat->get_nc()*4]; // for the 4 directions.
        zero<double>(hopping, lat->get_volume()*lat->get_nc()*lat->get_nc()*4);
        if (has_two)
        {
            two_link = new complex<double>[lat->get_volume()*lat->get_nc()*lat->get_nc()*8]; // for 4 two-link + 4 corner.
            zero<double>(two_link, lat->get_volume()*lat->get_nc()*lat->get_nc()*8);
        }
        else
        {
            two_link = 0;
        }
        
        sdir = DIR_ALL; 
    }
    
    // Copy constructor
    stencil_2d(const stencil_2d &obj)
    {
        lat = obj.lat;
        sdir = obj.sdir;
        has_two = obj.has_two;
        stencil_size = obj.stencil_size;
        shift = obj.shift; 
        eo_shift = obj.eo_shift;
        dof_shift = obj.dof_shift; 
        
        clover = new complex<double>[lat->get_volume()*lat->get_nc()*lat->get_nc()];
        copy<double>(clover, obj.clover, lat->get_volume()*lat->get_nc()*lat->get_nc());
        hopping = new complex<double>[lat->get_volume()*lat->get_nc()*lat->get_nc()*4]; // for the 4 directions.
        copy<double>(hopping, obj.hopping, lat->get_volume()*lat->get_nc()*lat->get_nc()*4);
        if (has_two)
        {
            two_link = new complex<double>[lat->get_volume()*lat->get_nc()*lat->get_nc()*8]; // for 4 two-link + 4 corner.
            copy<double>(two_link, obj.two_link, lat->get_volume()*lat->get_nc()*lat->get_nc()*8);
        }
        else
        {
            two_link = 0;
        }
        
        generated = obj.generated;
    }
    
    ~stencil_2d()
    {
        if (clover != 0) { delete[] clover; clover = 0; }
        if (hopping != 0) { delete[] hopping; hopping = 0; }
        if (has_two)
        {
            if (two_link != 0) { delete[] two_link; two_link = 0; }
        }
        generated = false; 

    }
    
    // Clear out stencils!
    void clear_stencils()
    {
        zero<double>(clover, lat->get_volume()*lat->get_nc()*lat->get_nc());
        zero<double>(hopping, lat->get_volume()*lat->get_nc()*lat->get_nc()*4);
        if (has_two)
        {
            zero<double>(two_link, lat->get_volume()*lat->get_nc()*lat->get_nc()*8);
        }
        generated = false; 
    }
    
    
};

// Applies an (up to two link) stencil as defined by a stencil_2d object.
void apply_stencil_2d(complex<double>* lhs, complex<double>* rhs, void* extra_data);

// Applies the eo piece of an (up to 2 link) stencil as defined by a stencil_2d object.
void apply_stencil_2d_eo(complex<double>* lhs, complex<double>* rhs, void* extra_data);
    
// Applies the oe piece of an (up to 2 link) stencil as defined by a stencil_2d object.
void apply_stencil_2d_oe(complex<double>* lhs, complex<double>* rhs, void* extra_data);

// Applies the tb (takes bottom half of dof to top half of dof) of an (up to 2 link) stencil as defined by a stencil_2d object.
void apply_stencil_2d_tb(complex<double>* lhs, complex<double>* rhs, void* extra_data);

// Applies the bt (takes top half of dof to bottom half of dof) of an (up to 2 link) stencil as defined by a stencil_2d object.
void apply_stencil_2d_bt(complex<double>* lhs, complex<double>* rhs, void* extra_data);

// Generate a stencil operator given a function and a coarsening distance (max 2...)
void generate_stencil_2d(stencil_2d* stenc, void (*matrix_vector)(complex<double>*,complex<double>*,void*), void* extra_data);


#endif // COARSE_STENCIL
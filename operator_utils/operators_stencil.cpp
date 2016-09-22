// Build and return the stencil (look in the stencil directory) for a U(1) operator.
// Depends on the lattice class and the staggered_u1_op struct.

#include <complex>
using std::complex;

#include "operators.h"
#include "lattice.h"
#include "coarse_stencil.h"
#include "operators_stencil.h"


// Given an allocated but otherwise empty stencil, appropriately fill it with the 2d staggered links.
void get_square_staggered_u1_stencil(stencil_2d* stenc, staggered_u1_op* stagif)
{
    if (stenc->generated) // If it's already generated, leave!
    {
        return;
    }
    
    if (stenc->lat->get_nc() != 1) // If it's not an nc = 1 lattice, leave!
    {
        return;
    }
    
    // Loop over all sites.
    int i;
    int lattice_size = stenc->lat->get_lattice_size();
    int coords[2], color; // color is always 0.
    int coords_tmp[2];
    int eta1; 
    
    for (i = 0; i < lattice_size; i++)
    {
        stenc->lat->index_to_coord(i, coords, color);
        eta1 = 1 - 2*(coords[0]%2);
        
        // Clover is the mass.
        stenc->clover[i] = stagif->mass;
        
        // +x
        stenc->hopping[i] = -0.5*stagif->lattice[2*i];
        
        // +y
        stenc->hopping[i+lattice_size] = -0.5*eta1*stagif->lattice[2*i+1];
        
        // -x
        coords_tmp[0] = (coords[0]-1+stenc->lat->get_lattice_dimension(0))%stenc->lat->get_lattice_dimension(0);
        coords_tmp[1] = coords[1];
        stenc->hopping[i+2*lattice_size] = 0.5*conj(stagif->lattice[2*stenc->lat->coord_to_index(coords_tmp, 0)]);
        
        // -y
        coords_tmp[0] = coords[0];
        coords_tmp[1] = (coords[1]-1+stenc->lat->get_lattice_dimension(1))%stenc->lat->get_lattice_dimension(1);
        stenc->hopping[i+3*lattice_size] = 0.5*eta1*conj(stagif->lattice[2*stenc->lat->coord_to_index(coords_tmp, 0)+1]);
    }
    
    stenc->generated = true;
}

void get_square_staggered_gamma5_u1_stencil(stencil_2d* stenc, staggered_u1_op* stagif)
{
    if (stenc->generated) // If it's already generated, leave!
    {
        return;
    }
    
    if (stenc->lat->get_nc() != 1) // If it's not an nc = 1 lattice, leave!
    {
        return;
    }
    
    // Loop over all sites.
    int i;
    int eo_sign;
    int lattice_size = stenc->lat->get_lattice_size();
    int coords[2], color; // color is always 0.
    int coords_tmp[2];
    int eta1; 
    
    for (i = 0; i < lattice_size; i++)
    {
        stenc->lat->index_to_coord(i, coords, color);
        eta1 = 1 - 2*(coords[0]%2);
        eo_sign = stenc->lat->index_is_even(i) ? -1 : 1;
        
        // Clover is the mass.
        stenc->clover[i] = eo_sign*stagif->mass;
        
        // +x
        stenc->hopping[i] = -0.5*eo_sign*stagif->lattice[2*i];
        
        // +y
        stenc->hopping[i+lattice_size] = -0.5*eo_sign*eta1*stagif->lattice[2*i+1];
        
        // -x
        coords_tmp[0] = (coords[0]-1+stenc->lat->get_lattice_dimension(0))%stenc->lat->get_lattice_dimension(0);
        coords_tmp[1] = coords[1];
        stenc->hopping[i+2*lattice_size] = 0.5*eo_sign*conj(stagif->lattice[2*stenc->lat->coord_to_index(coords_tmp, 0)]);
        
        // -y
        coords_tmp[0] = coords[0];
        coords_tmp[1] = (coords[1]-1+stenc->lat->get_lattice_dimension(1))%stenc->lat->get_lattice_dimension(1);
        stenc->hopping[i+3*lattice_size] = 0.5*eo_sign*eta1*conj(stagif->lattice[2*stenc->lat->coord_to_index(coords_tmp, 0)+1]);
    }
    
    stenc->generated = true;
}

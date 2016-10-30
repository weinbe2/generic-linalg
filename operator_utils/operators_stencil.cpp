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
        //stenc->clover[i] = stagif->mass;
        
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
    
    stenc->shift = stagif->mass;
    stenc->eo_shift = 0.0;
    stenc->dof_shift = 0.0;
    
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
        eo_sign = stenc->lat->index_is_even(i) ? 1.0 : -1.0;
        
        // Clover is the mass.
        //stenc->clover[i] = eo_sign*stagif->mass;
        
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
    
    stenc->shift = 0.0;
    stenc->eo_shift = stagif->mass;
    stenc->dof_shift = 0.0;
    
    stenc->generated = true;
}

// Given an allocated but otherwise empty stencil, appropriately fill it with staggered dagger links.
// Same as staggered, but with the negative of the hopping term.
void get_square_staggered_dagger_u1_stencil(stencil_2d* stenc, staggered_u1_op* stagif)
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
        //stenc->clover[i] = stagif->mass;
        
        // +x
        stenc->hopping[i] = 0.5*stagif->lattice[2*i];
        
        // +y
        stenc->hopping[i+lattice_size] = 0.5*eta1*stagif->lattice[2*i+1];
        
        // -x
        coords_tmp[0] = (coords[0]-1+stenc->lat->get_lattice_dimension(0))%stenc->lat->get_lattice_dimension(0);
        coords_tmp[1] = coords[1];
        stenc->hopping[i+2*lattice_size] = -0.5*conj(stagif->lattice[2*stenc->lat->coord_to_index(coords_tmp, 0)]);
        
        // -y
        coords_tmp[0] = coords[0];
        coords_tmp[1] = (coords[1]-1+stenc->lat->get_lattice_dimension(1))%stenc->lat->get_lattice_dimension(1);
        stenc->hopping[i+3*lattice_size] = -0.5*eta1*conj(stagif->lattice[2*stenc->lat->coord_to_index(coords_tmp, 0)+1]);
    }
    
    stenc->shift = stagif->mass;
    stenc->eo_shift = 0.0;
    stenc->dof_shift = 0.0;
    
    stenc->generated = true;
}


/////////////
// Functions for applying specialized operator stencils.
/////////////


// Prepare an even rhs for an even/odd preconditioned solve.
// Takes in rhs_orig, returns rhs_e. 
void apply_square_staggered_eoprec_prepare_stencil(complex<double>* rhs_e, complex<double>* rhs_orig, stencil_2d* stenc)
{
  int lattice_size = stenc->lat->get_lattice_size();
  
  // Prepare rhs: m rhs_e - D_{eo} rhs_o
  apply_stencil_2d_eo(rhs_e, rhs_orig, stenc); // zeroes odd sites in rhs_e.
  for (int i = 0; i < lattice_size; i++)
  {
    if (stenc->lat->index_is_even(i))
    {
      rhs_e[i] = stenc->shift*rhs_orig[i] - rhs_e[i];
    }
  }
}

// Square staggered 2d operator w/ u1 function, m^2 - D_{eo} D_{oe} [zeroes odd explicitly]
void apply_square_staggered_m2mdeodoe_stencil(complex<double>* lhs, complex<double>* rhs, void* extra_data)
{
  stencil_2d* stenc = (stencil_2d*)extra_data;
  int lattice_size = stenc->lat->get_lattice_size();
  
  complex<double>* tmp = new complex<double>[lattice_size];
  
  apply_stencil_2d_oe(tmp, rhs, extra_data);
  apply_stencil_2d_eo(lhs, tmp, extra_data);

  for (int i = 0; i < lattice_size; i++)
  {
    if (stenc->lat->index_is_even(i))
    {
      lhs[i] = stenc->shift*stenc->shift*rhs[i] - lhs[i];
    }
  }

  delete[] tmp;
}

// Reconstruct the full solution for an even/odd preconditioned solve.
// Takes in lhs_e, rhs_o, returns lhs_full (copying over the even part from lhs_e)
void apply_square_staggered_eoprec_reconstruct_stencil(complex<double>* lhs_full, complex<double>* lhs_e, complex<double>* rhs_o, stencil_2d* stenc)
{
  int lattice_size = stenc->lat->get_lattice_size();
  double inv_mass = 1.0/real(stenc->shift);
  
  // Reconstruct odd: m^{-1}*(rhs_o - D_{oe} lhs2_e)
  apply_stencil_2d_oe(lhs_full, lhs_e, stenc);
  for (int i = 0; i < lattice_size; i++)
  {
    if (!stenc->lat->index_is_even(i)) // odd
    {
      lhs_full[i] = inv_mass*(rhs_o[i] - lhs_full[i]);
    }
    else
    {
      lhs_full[i] = lhs_e[i];
    }
  }
}




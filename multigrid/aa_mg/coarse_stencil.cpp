// Routines related to applying an (up to 2 link) coarse stencil

#include <iostream>
#include <complex>
    
#include "generic_vector.h"
#include "coarse_stencil.h"
    
using namespace std;

// Applies an (up to two link) stencil as defined by a stencil_2d object.
void apply_stencil_2d(complex<double>* lhs, complex<double>* rhs, void* extra_data)
{
    int i;
    int c; 
    int lattice_size; 
    int nc;
    int coords[2];
    int coords_tmp[2];
    int color;
    stencil_2d links = *(stencil_2d*)extra_data;
    
    lattice_size = links.lat->get_lattice_size();
    nc = links.lat->get_nc();
    
    // Well, there's no time like today.
    for (i = 0; i < lattice_size; i++)
    {
        lhs[i] = 0.0;
        // Apply the stencil!
        
        // Get the coordinate and color index.
        links.lat->index_to_coord(i, (int*)coords, color);
        
        // Step 1: Clover term!
        for (c = 0; c < nc; c++)
        {
            lhs[i] += links.clover[c + nc*i]*rhs[links.lat->coord_to_index((int*)coords, c)];
        }
        
        // Step 2: Hopping term. 
        
        // +xhat.
        coords_tmp[0] = (coords[0] + 1)%links.lat->get_lattice_dimension(0);
        coords_tmp[1] = coords[1];
        for (c = 0; c < nc; c++)
        {
            lhs[i] += links.hopping[c + nc*i + 0*nc*lattice_size]*rhs[links.lat->coord_to_index((int*)coords_tmp, c)];
        }
        
        // +yhat
        coords_tmp[0] = coords[0];
        coords_tmp[1] = (coords[1] + 1)%links.lat->get_lattice_dimension(1);
        for (c = 0; c < nc; c++)
        {
            lhs[i] += links.hopping[c + nc*i + 1*nc*lattice_size]*rhs[links.lat->coord_to_index((int*)coords_tmp, c)];
        }
        
        // -xhat
        coords_tmp[0] = (coords[0] - 1 + links.lat->get_lattice_dimension(0))%links.lat->get_lattice_dimension(0);
        coords_tmp[1] = coords[1];
        for (c = 0; c < nc; c++)
        {
            lhs[i] += links.hopping[c + nc*i + 2*nc*lattice_size]*rhs[links.lat->coord_to_index((int*)coords_tmp, c)];
        }
        
        // -yhat
        coords_tmp[0] = coords[0];
        coords_tmp[1] = (coords[1] - 1 + links.lat->get_lattice_dimension(1))%links.lat->get_lattice_dimension(1);
        for (c = 0; c < nc; c++)
        {
            lhs[i] += links.hopping[c + nc*i + 3*nc*lattice_size]*rhs[links.lat->coord_to_index((int*)coords_tmp, c)];
        }
        
        // Step 3: Two-link term, if need be.
        
        if (links.has_two)
        {
        
            // +2 xhat.
            coords_tmp[0] = (coords[0] + 2)%links.lat->get_lattice_dimension(0);
            coords_tmp[1] = coords[1];
            for (c = 0; c < nc; c++)
            {
                lhs[i] += links.two_link[c + nc*i+0*nc*links.lat->coord_to_index((int*)coords,0)]*rhs[links.lat->coord_to_index((int*)coords_tmp, c)];
            }

            // +xhat+yhat
            coords_tmp[0] = (coords[0] + 1)%links.lat->get_lattice_dimension(0);
            coords_tmp[1] = (coords[1] + 1)%links.lat->get_lattice_dimension(1);
            for (c = 0; c < nc; c++)
            {
                lhs[i] += links.two_link[c + nc*i+1*nc*links.lat->coord_to_index((int*)coords,0)]*rhs[links.lat->coord_to_index((int*)coords_tmp, c)];
            }

            // +2 yhat
            coords_tmp[0] = coords[0];
            coords_tmp[1] = (coords[1] + 2)%links.lat->get_lattice_dimension(1);
            for (c = 0; c < nc; c++)
            {
                lhs[i] += links.two_link[c + nc*i+2*nc*links.lat->coord_to_index((int*)coords,0)]*rhs[links.lat->coord_to_index((int*)coords_tmp, c)];
            }

            // -xhat + yhat
            coords_tmp[0] = (coords[0] - 1 + 2*links.lat->get_lattice_dimension(0))%links.lat->get_lattice_dimension(0);
            coords_tmp[1] = (coords[1] + 1)%links.lat->get_lattice_dimension(1);
            for (c = 0; c < nc; c++)
            {
                lhs[i] += links.two_link[c + nc*i+3*nc*links.lat->coord_to_index((int*)coords,0)]*rhs[links.lat->coord_to_index((int*)coords_tmp, c)];
            }

            // -2 xhat.
            coords_tmp[0] = (coords[0] - 2 + 2*links.lat->get_lattice_dimension(0))%links.lat->get_lattice_dimension(0);
            coords_tmp[1] = coords[1];
            for (c = 0; c < nc; c++)
            {
                lhs[i] += links.two_link[c + nc*i+4*nc*links.lat->coord_to_index((int*)coords,0)]*rhs[links.lat->coord_to_index((int*)coords_tmp, c)];
            }

            // -xhat-yhat
            coords_tmp[0] = (coords[0] - 1 + 2*links.lat->get_lattice_dimension(0))%links.lat->get_lattice_dimension(0);
            coords_tmp[1] = (coords[1] - 1 + 2*links.lat->get_lattice_dimension(1))%links.lat->get_lattice_dimension(1);
            for (c = 0; c < nc; c++)
            {
                lhs[i] += links.two_link[c + nc*i+5*nc*links.lat->coord_to_index((int*)coords,0)]*rhs[links.lat->coord_to_index((int*)coords_tmp, c)];
            }

            // -2 yhat
            coords_tmp[0] = coords[0];
            coords_tmp[1] = (coords[1] - 2 + 2*links.lat->get_lattice_dimension(1))%links.lat->get_lattice_dimension(1);
            for (c = 0; c < nc; c++)
            {
                lhs[i] += links.two_link[c + nc*i+6*nc*links.lat->coord_to_index((int*)coords,0)]*rhs[links.lat->coord_to_index((int*)coords_tmp, c)];
            }

            // +xhat - yhat
            coords_tmp[0] = (coords[0] + 1)%links.lat->get_lattice_dimension(0);
            coords_tmp[1] = (coords[1] - 1 + 2*links.lat->get_lattice_dimension(1))%links.lat->get_lattice_dimension(1);
            for (c = 0; c < nc; c++)
            {
                lhs[i] += links.two_link[c + nc*i+7*nc*links.lat->coord_to_index((int*)coords,0)]*rhs[links.lat->coord_to_index((int*)coords_tmp, c)];
            }

        }
        
    }
}


// Generate a stencil operator given a function and a coarsening distance (max 2...)
// Currently not working...
void generate_stencil_2d(stencil_2d* stenc, int stencil_distance, void (*matrix_vector)(complex<double>*,complex<double>*,void*), void* extra_data)
{
    if (stenc->generated || stencil_distance > 2)
    {
        return;
    }
    
    Lattice* latt = stenc->lat;
    
    int i;
    int coord[2];
    int coord_tmp[2];
    int color;
    
    int latt_size = latt->get_lattice_size();
    int nc = latt->get_nc();
    
    stenc->clover = new complex<double>(latt_size*nc);
    stenc->hopping = new complex<double>(latt_size*nc*4); // 4 for the 4 directions
    if (stencil_distance == 2)
    {
        stenc->has_two = true;
        stenc->two_link = new complex<double>(latt_size*nc*8);
    }
    
    complex<double>* tmp_rhs = new complex<double>[latt_size];
    complex<double>* tmp_lhs = new complex<double>[latt_size];
    
    // Stencil. Currently inefficient---we apply a Dslash for every site on the lattice. Bad!
    for (i = 0; i < latt_size; i++)
    {
        latt->index_to_coord(i, (int*)coord, color);
        
        // Place a site!
        zero<double>(tmp_rhs, latt_size);
        tmp_rhs[i] = 1.0;
        
        zero<double>(tmp_lhs, latt_size);
        (*matrix_vector)(tmp_lhs, tmp_rhs, extra_data);
        
        for (int c = 0; c < nc; c++)
        {
            // clover.
            stenc->clover[color+nc*latt->coord_to_index((int*)coord, c)] = tmp_lhs[latt->coord_to_index((int*)coord, c)];
            
            
            // hopping.
            
            // +x 
            coord_tmp[0] = (coord[0]-1+latt->get_lattice_dimension(0))%latt->get_lattice_dimension(0);
            coord_tmp[1] = coord[1];
            
            stenc->hopping[color+nc*latt->coord_to_index((int*)coord_tmp,c)+0*nc*latt_size] = tmp_lhs[latt->coord_to_index((int*)coord_tmp, c)];
            
            cout << "Built hop +x " << i << "\n" << flush;
            
            // +y
            coord_tmp[0] = coord[0];
            coord_tmp[1] = (coord[1]-1+latt->get_lattice_dimension(1))%latt->get_lattice_dimension(1);
            
            stenc->hopping[color+nc*latt->coord_to_index((int*)coord_tmp,c)+1*nc*latt_size] = tmp_lhs[latt->coord_to_index((int*)coord_tmp, c)];
            
            cout << "Built hop +y " << i << "\n" << flush;
            
            // -x 
            coord_tmp[0] = (coord[0]+1)%latt->get_lattice_dimension(0);
            coord_tmp[1] = coord[1];
            stenc->hopping[color+nc*latt->coord_to_index((int*)coord_tmp,c)+2*nc*latt_size] = tmp_lhs[latt->coord_to_index((int*)coord_tmp, c)];
            
            cout << "Built hop -x " << i << "\n" << flush;
            
            // -y
            coord_tmp[0] = coord[0];
            coord_tmp[1] = (coord[1]+1)%latt->get_lattice_dimension(1);
            stenc->hopping[color+nc*latt->coord_to_index((int*)coord_tmp,c)+3*nc*latt_size] = tmp_lhs[latt->coord_to_index((int*)coord_tmp, c)];
            
            cout << "Built hop -y " << i << "\n" << flush;
            
            // two link.
            if (stenc->has_two)
            {
                // +2x 
                coord_tmp[0] = (coord[0]-2+2*latt->get_lattice_dimension(0))%latt->get_lattice_dimension(0);
                coord_tmp[1] = coord[1];
                stenc->two_link[color+nc*latt->coord_to_index((int*)coord_tmp,c)+0*nc*latt_size] = tmp_lhs[latt->coord_to_index((int*)coord_tmp, c)];

                // +x+y
                coord_tmp[0] = (coord[0]-1+2*latt->get_lattice_dimension(0))%latt->get_lattice_dimension(0);
                coord_tmp[1] = (coord[1]-1+2*latt->get_lattice_dimension(1))%latt->get_lattice_dimension(1);
                stenc->two_link[color+nc*latt->coord_to_index((int*)coord_tmp,c)+1*nc*latt_size] = tmp_lhs[latt->coord_to_index((int*)coord_tmp, c)];

                // +2y
                coord_tmp[0] = coord[0];
                coord_tmp[1] = (coord[1]-2+2*latt->get_lattice_dimension(1))%latt->get_lattice_dimension(1);
                stenc->two_link[color+nc*latt->coord_to_index((int*)coord_tmp,c)+2*nc*latt_size] = tmp_lhs[latt->coord_to_index((int*)coord_tmp, c)];

                // -x+y
                coord_tmp[0] = (coord[0]+1)%latt->get_lattice_dimension(0);
                coord_tmp[1] = (coord[1]-1+2*latt->get_lattice_dimension(1))%latt->get_lattice_dimension(1);
                stenc->two_link[color+nc*latt->coord_to_index((int*)coord_tmp,c)+3*nc*latt_size] = tmp_lhs[latt->coord_to_index((int*)coord_tmp, c)];
                
                // -2x 
                coord_tmp[0] = (coord[0]+2)%latt->get_lattice_dimension(0);
                coord_tmp[1] = coord[1];
                stenc->two_link[color+nc*latt->coord_to_index((int*)coord_tmp,c)+4*nc*latt_size] = tmp_lhs[latt->coord_to_index((int*)coord_tmp, c)];

                // -x-y
                coord_tmp[0] = (coord[0]+1)%latt->get_lattice_dimension(0);
                coord_tmp[1] = (coord[1]+1)%latt->get_lattice_dimension(1);
                stenc->two_link[color+nc*latt->coord_to_index((int*)coord_tmp,c)+5*nc*latt_size] = tmp_lhs[latt->coord_to_index((int*)coord_tmp, c)];

                // -2y
                coord_tmp[0] = coord[0];
                coord_tmp[1] = (coord[1]+2)%latt->get_lattice_dimension(1);
                stenc->two_link[color+nc*latt->coord_to_index((int*)coord_tmp,c)+6*nc*latt_size] = tmp_lhs[latt->coord_to_index((int*)coord_tmp, c)];

                // +x-y
                coord_tmp[0] = (coord[0]-1+2*latt->get_lattice_dimension(0))%latt->get_lattice_dimension(0);
                coord_tmp[1] = (coord[1]+1)%latt->get_lattice_dimension(1);
                stenc->two_link[color+nc*latt->coord_to_index((int*)coord_tmp,c)+7*nc*latt_size] = tmp_lhs[latt->coord_to_index((int*)coord_tmp, c)];
            }
        }
    }
    
    delete[] tmp_lhs;
    delete[] tmp_rhs; 
    
    stenc->generated = true; 
    
}

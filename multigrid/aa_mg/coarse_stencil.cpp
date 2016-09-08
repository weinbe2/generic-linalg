// Routines related to applying an (up to 2 link) coarse stencil

#include <complex>
    
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
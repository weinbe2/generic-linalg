// Routines related to applying an (up to 2 link) coarse stencil

#include <iostream>
#include <complex>
    
#include "generic_vector.h"
#include "coarse_stencil.h"
#include "mg_complex.h"
    
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
    
    // Are we applying the full stencil, or just one direction?
    // While I don't like essentially copy+pasting code, 
    // I'd rather just do one "if" statement.
    if (links.sdir == DIR_ALL)
    {
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
                    lhs[i] += links.two_link[c + nc*i+0*nc*lattice_size]*rhs[links.lat->coord_to_index((int*)coords_tmp, c)];
                }

                // +xhat+yhat
                coords_tmp[0] = (coords[0] + 1)%links.lat->get_lattice_dimension(0);
                coords_tmp[1] = (coords[1] + 1)%links.lat->get_lattice_dimension(1);
                for (c = 0; c < nc; c++)
                {
                    lhs[i] += links.two_link[c + nc*i+1*nc*lattice_size]*rhs[links.lat->coord_to_index((int*)coords_tmp, c)];
                }

                // +2 yhat
                coords_tmp[0] = coords[0];
                coords_tmp[1] = (coords[1] + 2)%links.lat->get_lattice_dimension(1);
                for (c = 0; c < nc; c++)
                {
                    lhs[i] += links.two_link[c + nc*i+2*nc*lattice_size]*rhs[links.lat->coord_to_index((int*)coords_tmp, c)];
                }

                // -xhat + yhat
                coords_tmp[0] = (coords[0] - 1 + 2*links.lat->get_lattice_dimension(0))%links.lat->get_lattice_dimension(0);
                coords_tmp[1] = (coords[1] + 1)%links.lat->get_lattice_dimension(1);
                for (c = 0; c < nc; c++)
                {
                    lhs[i] += links.two_link[c + nc*i+3*nc*lattice_size]*rhs[links.lat->coord_to_index((int*)coords_tmp, c)];
                }

                // -2 xhat.
                coords_tmp[0] = (coords[0] - 2 + 2*links.lat->get_lattice_dimension(0))%links.lat->get_lattice_dimension(0);
                coords_tmp[1] = coords[1];
                for (c = 0; c < nc; c++)
                {
                    lhs[i] += links.two_link[c + nc*i+4*nc*lattice_size]*rhs[links.lat->coord_to_index((int*)coords_tmp, c)];
                }

                // -xhat-yhat
                coords_tmp[0] = (coords[0] - 1 + 2*links.lat->get_lattice_dimension(0))%links.lat->get_lattice_dimension(0);
                coords_tmp[1] = (coords[1] - 1 + 2*links.lat->get_lattice_dimension(1))%links.lat->get_lattice_dimension(1);
                for (c = 0; c < nc; c++)
                {
                    lhs[i] += links.two_link[c + nc*i+5*nc*lattice_size]*rhs[links.lat->coord_to_index((int*)coords_tmp, c)];
                }

                // -2 yhat
                coords_tmp[0] = coords[0];
                coords_tmp[1] = (coords[1] - 2 + 2*links.lat->get_lattice_dimension(1))%links.lat->get_lattice_dimension(1);
                for (c = 0; c < nc; c++)
                {
                    lhs[i] += links.two_link[c + nc*i+6*nc*lattice_size]*rhs[links.lat->coord_to_index((int*)coords_tmp, c)];
                }

                // +xhat - yhat
                coords_tmp[0] = (coords[0] + 1)%links.lat->get_lattice_dimension(0);
                coords_tmp[1] = (coords[1] - 1 + 2*links.lat->get_lattice_dimension(1))%links.lat->get_lattice_dimension(1);
                for (c = 0; c < nc; c++)
                {
                    lhs[i] += links.two_link[c + nc*i+7*nc*lattice_size]*rhs[links.lat->coord_to_index((int*)coords_tmp, c)];
                }

            }

        }
    }
    else // just one piece of the stencil. Allllll the copy and pasting. 
    {
        switch (links.sdir)
        {
            case DIR_ALL: // Well, this can't happen.
                break;
            case DIR_0: // clover term
                for (i = 0; i < lattice_size; i++)
                {
                    lhs[i] = 0.0; // zero out.
                    links.lat->index_to_coord(i, (int*)coords, color); // get the coordinate and color index.
                    for (c = 0; c < nc; c++)
                    {
                        lhs[i] += links.clover[c + nc*i]*rhs[links.lat->coord_to_index((int*)coords, c)];
                    }
                }
                break;
            case DIR_XP1: // +x
                for (i = 0; i < lattice_size; i++)
                {
                    lhs[i] = 0.0; // zero out.
                    links.lat->index_to_coord(i, (int*)coords, color); // get the coordinate and color index.
                    coords_tmp[0] = (coords[0] + 1)%links.lat->get_lattice_dimension(0);
                    coords_tmp[1] = coords[1];
                    for (c = 0; c < nc; c++)
                    {
                        lhs[i] += links.hopping[c + nc*i + 0*nc*lattice_size]*rhs[links.lat->coord_to_index((int*)coords_tmp, c)];
                    }
                }
                break;
            case DIR_YP1: // +y
                for (i = 0; i < lattice_size; i++)
                {
                    lhs[i] = 0.0; // zero out.
                    links.lat->index_to_coord(i, (int*)coords, color); // get the coordinate and color index.
                    coords_tmp[0] = coords[0];
                    coords_tmp[1] = (coords[1] + 1)%links.lat->get_lattice_dimension(1);
                    for (c = 0; c < nc; c++)
                    {
                        lhs[i] += links.hopping[c + nc*i + 1*nc*lattice_size]*rhs[links.lat->coord_to_index((int*)coords_tmp, c)];
                    }
                }
                break;
            case DIR_XM1: // -x
                for (i = 0; i < lattice_size; i++)
                {
                    lhs[i] = 0.0; // zero out.
                    links.lat->index_to_coord(i, (int*)coords, color); // get the coordinate and color index.
                    coords_tmp[0] = (coords[0] - 1 + links.lat->get_lattice_dimension(0))%links.lat->get_lattice_dimension(0);
                    coords_tmp[1] = coords[1];
                    for (c = 0; c < nc; c++)
                    {
                        lhs[i] += links.hopping[c + nc*i + 2*nc*lattice_size]*rhs[links.lat->coord_to_index((int*)coords_tmp, c)];
                    }
                }
                break;
            case DIR_YM1: //-y
                for (i = 0; i < lattice_size; i++)
                {
                    lhs[i] = 0.0; // zero out.
                    links.lat->index_to_coord(i, (int*)coords, color); // get the coordinate and color index.
                    coords_tmp[0] = coords[0];
                    coords_tmp[1] = (coords[1] - 1 + links.lat->get_lattice_dimension(1))%links.lat->get_lattice_dimension(1);
                    for (c = 0; c < nc; c++)
                    {
                        lhs[i] += links.hopping[c + nc*i + 3*nc*lattice_size]*rhs[links.lat->coord_to_index((int*)coords_tmp, c)];
                    }
                }
                break;
            case DIR_XP2: // +2x
                if (links.has_two)
                {
                    for (i = 0; i < lattice_size; i++)
                    {
                        lhs[i] = 0.0; // zero out.
                        links.lat->index_to_coord(i, (int*)coords, color); // get the coordinate and color index.
                        coords_tmp[0] = (coords[0] + 2)%links.lat->get_lattice_dimension(0);
                        coords_tmp[1] = coords[1];
                        for (c = 0; c < nc; c++)
                        {
                            lhs[i] += links.two_link[c + nc*i+0*nc*lattice_size]*rhs[links.lat->coord_to_index((int*)coords_tmp, c)];
                        }
                    }
                }
                break;
            case DIR_XP1YP1: // +x+y
                if (links.has_two)
                {
                    for (i = 0; i < lattice_size; i++)
                    {
                        lhs[i] = 0.0; // zero out.
                        links.lat->index_to_coord(i, (int*)coords, color); // get the coordinate and color index.
                        coords_tmp[0] = (coords[0] + 1)%links.lat->get_lattice_dimension(0);
                        coords_tmp[1] = (coords[1] + 1)%links.lat->get_lattice_dimension(1);
                        for (c = 0; c < nc; c++)
                        {
                            lhs[i] += links.two_link[c + nc*i+1*nc*lattice_size]*rhs[links.lat->coord_to_index((int*)coords_tmp, c)];
                        }
                    }
                }
                break;
            case DIR_YP2: // +2y
                if (links.has_two)
                {
                    for (i = 0; i < lattice_size; i++)
                    {
                        lhs[i] = 0.0; // zero out.
                        links.lat->index_to_coord(i, (int*)coords, color); // get the coordinate and color index.
                        coords_tmp[0] = coords[0];
                        coords_tmp[1] = (coords[1] + 2)%links.lat->get_lattice_dimension(1);
                        for (c = 0; c < nc; c++)
                        {
                            lhs[i] += links.two_link[c + nc*i+2*nc*lattice_size]*rhs[links.lat->coord_to_index((int*)coords_tmp, c)];
                        }
                    }
                }
                break;
            case DIR_XM1YP1: // -x+y
                if (links.has_two)
                {
                    for (i = 0; i < lattice_size; i++)
                    {
                        lhs[i] = 0.0; // zero out.
                        links.lat->index_to_coord(i, (int*)coords, color); // get the coordinate and color index.
                        coords_tmp[0] = (coords[0] - 1 + 2*links.lat->get_lattice_dimension(0))%links.lat->get_lattice_dimension(0);
                        coords_tmp[1] = (coords[1] + 1)%links.lat->get_lattice_dimension(1);
                        for (c = 0; c < nc; c++)
                        {
                            lhs[i] += links.two_link[c + nc*i+3*nc*lattice_size]*rhs[links.lat->coord_to_index((int*)coords_tmp, c)];
                        }

                    }
                }
                break;
            case DIR_XM2: // -2x
                if (links.has_two)
                {
                    for (i = 0; i < lattice_size; i++)
                    {
                        lhs[i] = 0.0; // zero out.
                        links.lat->index_to_coord(i, (int*)coords, color); // get the coordinate and color index.
                        coords_tmp[0] = (coords[0] - 2 + 2*links.lat->get_lattice_dimension(0))%links.lat->get_lattice_dimension(0);
                        coords_tmp[1] = coords[1];
                        for (c = 0; c < nc; c++)
                        {
                            lhs[i] += links.two_link[c + nc*i+4*nc*lattice_size]*rhs[links.lat->coord_to_index((int*)coords_tmp, c)];
                        }
                    }
                }
                break;
            case DIR_XM1YM1: // -x-y
                if (links.has_two)
                {
                    for (i = 0; i < lattice_size; i++)
                    {
                        lhs[i] = 0.0; // zero out.
                        links.lat->index_to_coord(i, (int*)coords, color); // get the coordinate and color index.
                        coords_tmp[0] = (coords[0] - 1 + 2*links.lat->get_lattice_dimension(0))%links.lat->get_lattice_dimension(0);
                        coords_tmp[1] = (coords[1] - 1 + 2*links.lat->get_lattice_dimension(1))%links.lat->get_lattice_dimension(1);
                        for (c = 0; c < nc; c++)
                        {
                            lhs[i] += links.two_link[c + nc*i+5*nc*lattice_size]*rhs[links.lat->coord_to_index((int*)coords_tmp, c)];
                        }
                    }
                }
                break;
            case DIR_YM2: // -2y
                if (links.has_two)
                {
                    for (i = 0; i < lattice_size; i++)
                    {
                        lhs[i] = 0.0; // zero out.
                        links.lat->index_to_coord(i, (int*)coords, color); // get the coordinate and color index.
                        coords_tmp[0] = coords[0];
                        coords_tmp[1] = (coords[1] - 2 + 2*links.lat->get_lattice_dimension(1))%links.lat->get_lattice_dimension(1);
                        for (c = 0; c < nc; c++)
                        {
                            lhs[i] += links.two_link[c + nc*i+6*nc*lattice_size]*rhs[links.lat->coord_to_index((int*)coords_tmp, c)];
                        }
                    }
                }
                break;
            case DIR_XP1YM1: // +x-y
                if (links.has_two)
                {
                    for (i = 0; i < lattice_size; i++)
                    {
                        lhs[i] = 0.0; // zero out.
                        links.lat->index_to_coord(i, (int*)coords, color); // get the coordinate and color index.
                        coords_tmp[0] = (coords[0] + 1)%links.lat->get_lattice_dimension(0);
                        coords_tmp[1] = (coords[1] - 1 + 2*links.lat->get_lattice_dimension(1))%links.lat->get_lattice_dimension(1);
                        for (c = 0; c < nc; c++)
                        {
                            lhs[i] += links.two_link[c + nc*i+7*nc*lattice_size]*rhs[links.lat->coord_to_index((int*)coords_tmp, c)];
                        }
                    }
                }
                break;
        } // end switch
    } // end only one link for the stencil. 
}


// Generate a stencil operator given a function and a coarsening distance (max 2...). Assumes that the physical lengths are longer than coarsening distance*2. 
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
    
    stenc->clover = new complex<double>[latt_size*nc];
    zero<double>(stenc->clover, latt_size*nc);
    
    stenc->hopping = new complex<double>[latt_size*nc*4]; // 4 for the 4 directions
    zero<double>(stenc->hopping, latt_size*nc*4);
    if (stencil_distance == 2)
    {
        stenc->has_two = true;
        stenc->two_link = new complex<double>[latt_size*nc*8];
        zero<double>(stenc->two_link, latt_size*nc*8);
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
            stenc->clover[color+nc*latt->coord_to_index((int*)coord, c)] += tmp_lhs[latt->coord_to_index((int*)coord, c)];
            
            
            // hopping.
            
            // +x 
            coord_tmp[0] = (coord[0]-1+latt->get_lattice_dimension(0))%latt->get_lattice_dimension(0);
            coord_tmp[1] = coord[1];
            stenc->hopping[color+nc*latt->coord_to_index((int*)coord_tmp,c)+0*nc*latt_size] += tmp_lhs[latt->coord_to_index((int*)coord_tmp, c)];
            
            // +y
            coord_tmp[0] = coord[0];
            coord_tmp[1] = (coord[1]-1+latt->get_lattice_dimension(1))%latt->get_lattice_dimension(1);
            
            stenc->hopping[color+nc*latt->coord_to_index((int*)coord_tmp,c)+1*nc*latt_size] += tmp_lhs[latt->coord_to_index((int*)coord_tmp, c)];
            
            // -x 
            coord_tmp[0] = (coord[0]+1)%latt->get_lattice_dimension(0);
            coord_tmp[1] = coord[1];
            stenc->hopping[color+nc*latt->coord_to_index((int*)coord_tmp,c)+2*nc*latt_size] += tmp_lhs[latt->coord_to_index((int*)coord_tmp, c)];

            // -y
            coord_tmp[0] = coord[0];
            coord_tmp[1] = (coord[1]+1)%latt->get_lattice_dimension(1);
            stenc->hopping[color+nc*latt->coord_to_index((int*)coord_tmp,c)+3*nc*latt_size] += tmp_lhs[latt->coord_to_index((int*)coord_tmp, c)];
            
            // two link.
            if (stenc->has_two)
            {
                // +2x 
                coord_tmp[0] = (coord[0]-2+2*latt->get_lattice_dimension(0))%latt->get_lattice_dimension(0);
                coord_tmp[1] = coord[1];
                stenc->two_link[color+nc*latt->coord_to_index((int*)coord_tmp,c)+0*nc*latt_size] += tmp_lhs[latt->coord_to_index((int*)coord_tmp, c)];

                // +x+y
                coord_tmp[0] = (coord[0]-1+2*latt->get_lattice_dimension(0))%latt->get_lattice_dimension(0);
                coord_tmp[1] = (coord[1]-1+2*latt->get_lattice_dimension(1))%latt->get_lattice_dimension(1);
                stenc->two_link[color+nc*latt->coord_to_index((int*)coord_tmp,c)+1*nc*latt_size] += tmp_lhs[latt->coord_to_index((int*)coord_tmp, c)];

                // +2y
                coord_tmp[0] = coord[0];
                coord_tmp[1] = (coord[1]-2+2*latt->get_lattice_dimension(1))%latt->get_lattice_dimension(1);
                stenc->two_link[color+nc*latt->coord_to_index((int*)coord_tmp,c)+2*nc*latt_size] += tmp_lhs[latt->coord_to_index((int*)coord_tmp, c)];

                // -x+y
                coord_tmp[0] = (coord[0]+1)%latt->get_lattice_dimension(0);
                coord_tmp[1] = (coord[1]-1+2*latt->get_lattice_dimension(1))%latt->get_lattice_dimension(1);
                stenc->two_link[color+nc*latt->coord_to_index((int*)coord_tmp,c)+3*nc*latt_size] += tmp_lhs[latt->coord_to_index((int*)coord_tmp, c)];
                
                // -2x 
                coord_tmp[0] = (coord[0]+2)%latt->get_lattice_dimension(0);
                coord_tmp[1] = coord[1];
                stenc->two_link[color+nc*latt->coord_to_index((int*)coord_tmp,c)+4*nc*latt_size] += tmp_lhs[latt->coord_to_index((int*)coord_tmp, c)];

                // -x-y
                coord_tmp[0] = (coord[0]+1)%latt->get_lattice_dimension(0);
                coord_tmp[1] = (coord[1]+1)%latt->get_lattice_dimension(1);
                stenc->two_link[color+nc*latt->coord_to_index((int*)coord_tmp,c)+5*nc*latt_size] += tmp_lhs[latt->coord_to_index((int*)coord_tmp, c)];

                // -2y
                coord_tmp[0] = coord[0];
                coord_tmp[1] = (coord[1]+2)%latt->get_lattice_dimension(1);
                stenc->two_link[color+nc*latt->coord_to_index((int*)coord_tmp,c)+6*nc*latt_size] += tmp_lhs[latt->coord_to_index((int*)coord_tmp, c)];

                // +x-y
                coord_tmp[0] = (coord[0]-1+2*latt->get_lattice_dimension(0))%latt->get_lattice_dimension(0);
                coord_tmp[1] = (coord[1]+1)%latt->get_lattice_dimension(1);
                stenc->two_link[color+nc*latt->coord_to_index((int*)coord_tmp,c)+7*nc*latt_size] += tmp_lhs[latt->coord_to_index((int*)coord_tmp, c)];
            }
        }
    }
    
    delete[] tmp_lhs;
    delete[] tmp_rhs; 
    
    stenc->generated = true; 
    
}


// Generate a coarse stencil from a fine stencil. This takes advantage of the prolong and restrict functions
// explicitly, and also depends on the "sdir" variable the stencil object includes.
void generate_coarse_from_fine_stencil(stencil_2d* stenc_coarse, stencil_2d* stenc_fine, mg_operator_struct_complex* mgstruct, int stencil_distance)
{
    if (stenc_coarse->generated || stencil_distance > 2)
    {
        return;
    }
    
    Lattice* latt_fine = stenc_fine->lat;
    Lattice* latt_coarse = stenc_coarse->lat; 
    
    int i,x,y;
    int coord[2];
    int color, c;
    
    int fine_size = latt_fine->get_lattice_size(); // Only required to allocate temporary vectors. 
    int coarse_size = latt_coarse->get_lattice_size();
    int coarse_nc = latt_coarse->get_nc();
    
    
    // Allocate the coarse stencil. 
    stenc_coarse->clover = new complex<double>[coarse_size*coarse_nc];
    zero<double>(stenc_coarse->clover, coarse_size*coarse_nc);
    
    stenc_coarse->hopping = new complex<double>[coarse_size*coarse_nc*4]; // 4 for the 4 directions
    zero<double>(stenc_coarse->hopping, coarse_size*coarse_nc*4);
    if (stencil_distance == 2)
    {
        stenc_coarse->has_two = true;
        stenc_coarse->two_link = new complex<double>[coarse_size*coarse_nc*8];
        zero<double>(stenc_coarse->two_link, coarse_size*coarse_nc*8);
    }
    
    // Coarse vectors.
    complex<double>* tmp_rhs = new complex<double>[coarse_size];
    complex<double>* tmp_lhs = new complex<double>[coarse_size];
    
    // Fine vectors.
    complex<double>* tmp_Prhs = new complex<double>[fine_size];
    complex<double>* tmp_APrhs = new complex<double>[fine_size];
    
    
    // Save the state of the fine stencil.
    stencil_dir saved_dir = stenc_fine->sdir; 
    
    // We build the coarse stencil one piece at a time: the fine clover, the fine hopping terms, the fine two-link terms.
    
    // Piece 1: Fine Clover.
    // For the fine clover, we learn about some of the coarse clover. This takes one application per coarse color.
    stenc_fine->sdir = DIR_0;
    for (color = 0; color < coarse_nc; color++)
    {
        zero<double>(tmp_rhs, coarse_size);
        
        // Set a 1 at each coarse site.
        for (i = 0; i < latt_coarse->get_volume(); i++)
        {
            tmp_rhs[i*coarse_nc+color] = 1.0;
        }
        
        // Apply P^\dag A_0 P
        prolong(tmp_Prhs, tmp_rhs, mgstruct);
        apply_stencil_2d(tmp_APrhs, tmp_Prhs, (void*)stenc_fine);
        restrict(tmp_lhs, tmp_APrhs, mgstruct);
        
        // Loop over each coarse site, update the clover.
        for (i = 0; i < latt_coarse->get_lattice_size(); i++)
        {
            latt_coarse->index_to_coord(i, (int*)coord, c);
            
            stenc_coarse->clover[color+coarse_nc*latt_coarse->coord_to_index((int*)coord, c)] += tmp_lhs[latt_coarse->coord_to_index(coord, c)];
        }
    }
    
    // Piece 2: The Hopping Term.
    // For the fine hopping term, we learn about more of the coarse clover, and some of the coarse hopping term.
    // In the case of a distance = 1 fine stencil, we learn about the entire coarse hopping term,
    // and also finish learning about the coarse clover.
    // We learn about the hopping term in two hits per direction: one for even sites, one for odd sites.
    for (color = 0; color < coarse_nc; color++)
    {   
        // First step: Learn about the +x part of the stencil.
        // One for the even x sites, one for the odd x sites.
        {
            stenc_fine->sdir = DIR_XP1;


            // Even sites first.
            zero<double>(tmp_rhs, coarse_size);
            for (i = 0; i < latt_coarse->get_volume(); i+=2)
            {
                tmp_rhs[i*coarse_nc+color] = 1.0;
            }

            // Apply P^\dag A_+x P
            prolong(tmp_Prhs, tmp_rhs, mgstruct);
            apply_stencil_2d(tmp_APrhs, tmp_Prhs, (void*)stenc_fine);
            restrict(tmp_lhs, tmp_APrhs, mgstruct);

            // Loop over each coarse site, update the clover if it's an even site, hopping if it's an odd site.
            for (i = 0; i < latt_coarse->get_lattice_size(); i++)
            {
                latt_coarse->index_to_coord(i, (int*)coord, c);

                if (coord[0] % 2 == 0) // even, clover
                {
                    stenc_coarse->clover[color+coarse_nc*latt_coarse->coord_to_index(coord, c)] += tmp_lhs[latt_coarse->coord_to_index((int*)coord, c)];
                }
                else // odd, hopping
                {
                    stenc_coarse->hopping[color+coarse_nc*latt_coarse->coord_to_index(coord,c)+0*coarse_nc*coarse_size] += tmp_lhs[latt_coarse->coord_to_index(coord, c)];
                }
            }

            // Odd sites second.
            zero<double>(tmp_rhs, coarse_size);
            for (i = 1; i < latt_coarse->get_volume(); i+=2)
            {
                tmp_rhs[i*coarse_nc+color] = 1.0;
            }

            // Apply P^\dag A_+x P
            prolong(tmp_Prhs, tmp_rhs, mgstruct);
            apply_stencil_2d(tmp_APrhs, tmp_Prhs, (void*)stenc_fine);
            restrict(tmp_lhs, tmp_APrhs, mgstruct);

            // Loop over each coarse site, update the clover if it's an even site, hopping if it's an odd site.
            for (i = 0; i < latt_coarse->get_lattice_size(); i++)
            {
                latt_coarse->index_to_coord(i, (int*)coord, c);

                if (coord[0] % 2 == 1) // odd, clover
                {
                    stenc_coarse->clover[color+coarse_nc*latt_coarse->coord_to_index((int*)coord, c)] += tmp_lhs[latt_coarse->coord_to_index((int*)coord, c)];
                }
                else // even, hopping
                {
                    stenc_coarse->hopping[color+coarse_nc*latt_coarse->coord_to_index(coord,c)+0*coarse_nc*coarse_size] += tmp_lhs[latt_coarse->coord_to_index(coord, c)];
                }
            }
        }
        // Done with +x.
        cout << "Generated +x stencil.\n" << flush; 
        
        // Second step: Learn about the -x part of the stencil.
        // One for the even x sites, one for the odd x sites.
        {
            stenc_fine->sdir = DIR_XM1;


            // Even sites first.
            zero<double>(tmp_rhs, coarse_size);
            for (i = 0; i < latt_coarse->get_volume(); i+=2)
            {
                tmp_rhs[i*coarse_nc+color] = 1.0;
            }

            // Apply P^\dag A_-x P
            prolong(tmp_Prhs, tmp_rhs, mgstruct);
            apply_stencil_2d(tmp_APrhs, tmp_Prhs, (void*)stenc_fine);
            restrict(tmp_lhs, tmp_APrhs, mgstruct);

            // Loop over each coarse site, update the clover if it's an even site, hopping if it's an odd site.
            for (i = 0; i < latt_coarse->get_lattice_size(); i++)
            {
                latt_coarse->index_to_coord(i, (int*)coord, c);

                if (coord[0] % 2 == 0) // even, clover
                {
                    stenc_coarse->clover[color+coarse_nc*latt_coarse->coord_to_index((int*)coord, c)] += tmp_lhs[latt_coarse->coord_to_index((int*)coord, c)];
                }
                else // odd, hopping
                {
                    stenc_coarse->hopping[color+coarse_nc*latt_coarse->coord_to_index((int*)coord,c)+2*coarse_nc*coarse_size] += tmp_lhs[latt_coarse->coord_to_index((int*)coord, c)];
                }
            }

            // Odd sites second.
            zero<double>(tmp_rhs, coarse_size);
            for (i = 1; i < latt_coarse->get_volume(); i+=2)
            {
                tmp_rhs[i*coarse_nc+color] = 1.0;
            }

            // Apply P^\dag A_-x P
            prolong(tmp_Prhs, tmp_rhs, mgstruct);
            apply_stencil_2d(tmp_APrhs, tmp_Prhs, (void*)stenc_fine);
            restrict(tmp_lhs, tmp_APrhs, mgstruct);

            // Loop over each coarse site, update the clover if it's an even site, hopping if it's an odd site.
            for (i = 0; i < latt_coarse->get_lattice_size(); i++)
            {
                latt_coarse->index_to_coord(i, (int*)coord, c);

                if (coord[0] % 2 == 1) // odd, clover
                {
                    stenc_coarse->clover[color+coarse_nc*latt_coarse->coord_to_index((int*)coord, c)] += tmp_lhs[latt_coarse->coord_to_index((int*)coord, c)];
                }
                else // even, hopping
                {
                    stenc_coarse->hopping[color+coarse_nc*latt_coarse->coord_to_index((int*)coord,c)+2*coarse_nc*coarse_size] += tmp_lhs[latt_coarse->coord_to_index((int*)coord, c)];
                }
            }
        }
        // Done with -x.
        cout << "Generated -x stencil.\n" << flush; 
        
        // Third step: Learn about the +y part of the stencil.
        // One for the even y sites, one for the odd y sites.
        {
            stenc_fine->sdir = DIR_YP1;


            // Even sites first.
            zero<double>(tmp_rhs, coarse_size);
            for (i = 0; i < latt_coarse->get_volume()-latt_coarse->get_lattice_dimension(0); i++)
            {
                if ((i/latt_coarse->get_lattice_dimension(0)) % 2 == 1) // if it has an odd y...
                {
                    i += latt_coarse->get_lattice_dimension(0);
                }
                tmp_rhs[i*coarse_nc+color] = 1.0;
            }

            // Apply P^\dag A_+y P
            prolong(tmp_Prhs, tmp_rhs, mgstruct);
            apply_stencil_2d(tmp_APrhs, tmp_Prhs, (void*)stenc_fine);
            restrict(tmp_lhs, tmp_APrhs, mgstruct);

            // Loop over each coarse site, update the clover if it's an even site, hopping if it's an odd site.
            for (i = 0; i < latt_coarse->get_lattice_size(); i++)
            {
                latt_coarse->index_to_coord(i, (int*)coord, c);

                if (coord[1] % 2 == 0) // even, clover
                {
                    stenc_coarse->clover[color+coarse_nc*latt_coarse->coord_to_index((int*)coord, c)] += tmp_lhs[latt_coarse->coord_to_index((int*)coord, c)];
                }
                else // odd, hopping
                {
                    stenc_coarse->hopping[color+coarse_nc*latt_coarse->coord_to_index((int*)coord,c)+1*coarse_nc*coarse_size] += tmp_lhs[latt_coarse->coord_to_index((int*)coord, c)];
                }
            }

            // Odd sites second.
            zero<double>(tmp_rhs, coarse_size);
            for (i = latt_coarse->get_lattice_dimension(0); i < latt_coarse->get_volume(); i++)
            {
                if ((i/latt_coarse->get_lattice_dimension(0)) % 2 == 0) // if it has an even y...
                {
                    i += latt_coarse->get_lattice_dimension(0);
                }
                tmp_rhs[i*coarse_nc+color] = 1.0;
            }

            // Apply P^\dag A_+y P
            prolong(tmp_Prhs, tmp_rhs, mgstruct);
            apply_stencil_2d(tmp_APrhs, tmp_Prhs, (void*)stenc_fine);
            restrict(tmp_lhs, tmp_APrhs, mgstruct);

            // Loop over each coarse site, update the clover if it's an even site, hopping if it's an odd site.
            for (i = 0; i < latt_coarse->get_lattice_size(); i++)
            {
                latt_coarse->index_to_coord(i, (int*)coord, c);

                if (coord[1] % 2 == 1) // odd, clover
                {
                    stenc_coarse->clover[color+coarse_nc*latt_coarse->coord_to_index((int*)coord, c)] += tmp_lhs[latt_coarse->coord_to_index((int*)coord, c)];
                }
                else // even, hopping
                {
                    stenc_coarse->hopping[color+coarse_nc*latt_coarse->coord_to_index((int*)coord,c)+1*coarse_nc*coarse_size] += tmp_lhs[latt_coarse->coord_to_index((int*)coord, c)];
                }
            }
        }
        // Done with +y.
        cout << "Generated +y stencil.\n" << flush; 
        
        
        // Second step: Learn about the -y part of the stencil.
        // One for the even y sites, one for the odd y sites.
        {
            stenc_fine->sdir = DIR_YM1;


            // Even sites first.
            zero<double>(tmp_rhs, coarse_size);
            for (i = 0; i < latt_coarse->get_volume()-latt_coarse->get_lattice_dimension(0); i++)
            {
                if ((i/latt_coarse->get_lattice_dimension(0)) % 2 == 1) // if it has an odd y...
                {
                    i += latt_coarse->get_lattice_dimension(0);
                }
                tmp_rhs[i*coarse_nc+color] = 1.0;
            }

            // Apply P^\dag A_-y P
            prolong(tmp_Prhs, tmp_rhs, mgstruct);
            apply_stencil_2d(tmp_APrhs, tmp_Prhs, (void*)stenc_fine);
            restrict(tmp_lhs, tmp_APrhs, mgstruct);

            // Loop over each coarse site, update the clover if it's an even site, hopping if it's an odd site.
            for (i = 0; i < latt_coarse->get_lattice_size(); i++)
            {
                latt_coarse->index_to_coord(i, (int*)coord, c);

                if (coord[1] % 2 == 0) // even, clover
                {
                    stenc_coarse->clover[color+coarse_nc*latt_coarse->coord_to_index((int*)coord, c)] += tmp_lhs[latt_coarse->coord_to_index((int*)coord, c)];
                }
                else // odd, hopping
                {
                    stenc_coarse->hopping[color+coarse_nc*latt_coarse->coord_to_index((int*)coord,c)+3*coarse_nc*coarse_size] += tmp_lhs[latt_coarse->coord_to_index((int*)coord, c)];
                }
            }

            // Odd sites second.
            zero<double>(tmp_rhs, coarse_size);
            for (i = latt_coarse->get_lattice_dimension(0); i < latt_coarse->get_volume(); i++)
            {
                if ((i/latt_coarse->get_lattice_dimension(0)) % 2 == 0) // if it has an even y...
                {
                    i += latt_coarse->get_lattice_dimension(0);
                }
                tmp_rhs[i*coarse_nc+color] = 1.0;
            }

            // Apply P^\dag A_-y P
            prolong(tmp_Prhs, tmp_rhs, mgstruct);
            apply_stencil_2d(tmp_APrhs, tmp_Prhs, (void*)stenc_fine);
            restrict(tmp_lhs, tmp_APrhs, mgstruct);

            // Loop over each coarse site, update the clover if it's an even site, hopping if it's an odd site.
            for (i = 0; i < latt_coarse->get_lattice_size(); i++)
            {
                latt_coarse->index_to_coord(i, (int*)coord, c);

                if (coord[1] % 2 == 1) // odd, clover
                {
                    stenc_coarse->clover[color+coarse_nc*latt_coarse->coord_to_index((int*)coord, c)] += tmp_lhs[latt_coarse->coord_to_index((int*)coord, c)];
                }
                else // even, hopping
                {
                    stenc_coarse->hopping[color+coarse_nc*latt_coarse->coord_to_index((int*)coord,c)+3*coarse_nc*coarse_size] += tmp_lhs[latt_coarse->coord_to_index((int*)coord, c)];
                }
            }
        }
        // Done with -y.
        cout << "Generated -y stencil.\n" << flush;
        
        cout << "Finished color " << color << "\n" << flush; 
    }
         
    delete[] tmp_lhs;
    delete[] tmp_rhs; 
    
    delete[] tmp_Prhs;
    delete[] tmp_APrhs; 
    
    stenc_fine->generated = true; 
    
    // Restore the state of the fine stencil.
    stenc_fine->sdir = saved_dir; 
    
}




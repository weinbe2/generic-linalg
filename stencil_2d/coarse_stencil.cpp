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
    stencil_2d* links = (stencil_2d*)extra_data;
    
    lattice_size = links->lat->get_lattice_size();
    nc = links->lat->get_nc();
    
    // Are we applying the full stencil, or just one direction?
    // While I don't like essentially copy+pasting code, 
    // I'd rather just do one "if" statement.
    if (links->sdir == DIR_ALL)
    {
        // Well, there's no time like today.
        for (i = 0; i < lattice_size; i++)
        {
            lhs[i] = 0.0;
            // Apply the stencil!

            // Get the coordinate and color index.
            links->lat->index_to_coord(i, (int*)coords, color);


            // Step 1: Clover term!
            for (c = 0; c < nc; c++)
            {
                lhs[i] += links->clover[c + nc*i]*rhs[links->lat->coord_to_index((int*)coords, c)];
            }

            // Step 2: Hopping term. 

            // +xhat.
            coords_tmp[0] = (coords[0] + 1)%links->lat->get_lattice_dimension(0);
            coords_tmp[1] = coords[1];
            for (c = 0; c < nc; c++)
            {
                lhs[i] += links->hopping[c + nc*i + 0*nc*lattice_size]*rhs[links->lat->coord_to_index((int*)coords_tmp, c)];
            }

            // +yhat
            coords_tmp[0] = coords[0];
            coords_tmp[1] = (coords[1] + 1)%links->lat->get_lattice_dimension(1);
            for (c = 0; c < nc; c++)
            {
                lhs[i] += links->hopping[c + nc*i + 1*nc*lattice_size]*rhs[links->lat->coord_to_index((int*)coords_tmp, c)];
            }

            // -xhat
            coords_tmp[0] = (coords[0] - 1 + links->lat->get_lattice_dimension(0))%links->lat->get_lattice_dimension(0);
            coords_tmp[1] = coords[1];
            for (c = 0; c < nc; c++)
            {
                lhs[i] += links->hopping[c + nc*i + 2*nc*lattice_size]*rhs[links->lat->coord_to_index((int*)coords_tmp, c)];
            }

            // -yhat
            coords_tmp[0] = coords[0];
            coords_tmp[1] = (coords[1] - 1 + links->lat->get_lattice_dimension(1))%links->lat->get_lattice_dimension(1);
            for (c = 0; c < nc; c++)
            {
                lhs[i] += links->hopping[c + nc*i + 3*nc*lattice_size]*rhs[links->lat->coord_to_index((int*)coords_tmp, c)];
            }

            // Step 3: Two-link term, if need be.

            if (links->has_two)
            {

                // +2 xhat.
                coords_tmp[0] = (coords[0] + 2)%links->lat->get_lattice_dimension(0);
                coords_tmp[1] = coords[1];
                for (c = 0; c < nc; c++)
                {
                    lhs[i] += links->two_link[c + nc*i+0*nc*lattice_size]*rhs[links->lat->coord_to_index((int*)coords_tmp, c)];
                }

                // +xhat+yhat
                coords_tmp[0] = (coords[0] + 1)%links->lat->get_lattice_dimension(0);
                coords_tmp[1] = (coords[1] + 1)%links->lat->get_lattice_dimension(1);
                for (c = 0; c < nc; c++)
                {
                    lhs[i] += links->two_link[c + nc*i+1*nc*lattice_size]*rhs[links->lat->coord_to_index((int*)coords_tmp, c)];
                }

                // +2 yhat
                coords_tmp[0] = coords[0];
                coords_tmp[1] = (coords[1] + 2)%links->lat->get_lattice_dimension(1);
                for (c = 0; c < nc; c++)
                {
                    lhs[i] += links->two_link[c + nc*i+2*nc*lattice_size]*rhs[links->lat->coord_to_index((int*)coords_tmp, c)];
                }

                // -xhat + yhat
                coords_tmp[0] = (coords[0] - 1 + 2*links->lat->get_lattice_dimension(0))%links->lat->get_lattice_dimension(0);
                coords_tmp[1] = (coords[1] + 1)%links->lat->get_lattice_dimension(1);
                for (c = 0; c < nc; c++)
                {
                    lhs[i] += links->two_link[c + nc*i+3*nc*lattice_size]*rhs[links->lat->coord_to_index((int*)coords_tmp, c)];
                }

                // -2 xhat.
                coords_tmp[0] = (coords[0] - 2 + 2*links->lat->get_lattice_dimension(0))%links->lat->get_lattice_dimension(0);
                coords_tmp[1] = coords[1];
                for (c = 0; c < nc; c++)
                {
                    lhs[i] += links->two_link[c + nc*i+4*nc*lattice_size]*rhs[links->lat->coord_to_index((int*)coords_tmp, c)];
                }

                // -xhat-yhat
                coords_tmp[0] = (coords[0] - 1 + 2*links->lat->get_lattice_dimension(0))%links->lat->get_lattice_dimension(0);
                coords_tmp[1] = (coords[1] - 1 + 2*links->lat->get_lattice_dimension(1))%links->lat->get_lattice_dimension(1);
                for (c = 0; c < nc; c++)
                {
                    lhs[i] += links->two_link[c + nc*i+5*nc*lattice_size]*rhs[links->lat->coord_to_index((int*)coords_tmp, c)];
                }

                // -2 yhat
                coords_tmp[0] = coords[0];
                coords_tmp[1] = (coords[1] - 2 + 2*links->lat->get_lattice_dimension(1))%links->lat->get_lattice_dimension(1);
                for (c = 0; c < nc; c++)
                {
                    lhs[i] += links->two_link[c + nc*i+6*nc*lattice_size]*rhs[links->lat->coord_to_index((int*)coords_tmp, c)];
                }

                // +xhat - yhat
                coords_tmp[0] = (coords[0] + 1)%links->lat->get_lattice_dimension(0);
                coords_tmp[1] = (coords[1] - 1 + 2*links->lat->get_lattice_dimension(1))%links->lat->get_lattice_dimension(1);
                for (c = 0; c < nc; c++)
                {
                    lhs[i] += links->two_link[c + nc*i+7*nc*lattice_size]*rhs[links->lat->coord_to_index((int*)coords_tmp, c)];
                }

            }

        }
    }
    else // just one piece of the stencil. Allllll the copy and pasting. 
    {
        switch (links->sdir)
        {
            case DIR_ALL: // Well, this can't happen.
                break;
            case DIR_0: // clover term
                for (i = 0; i < lattice_size; i++)
                {
                    lhs[i] = 0.0; // zero out.
                    links->lat->index_to_coord(i, (int*)coords, color); // get the coordinate and color index.
                    for (c = 0; c < nc; c++)
                    {
                        lhs[i] += links->clover[c + nc*i]*rhs[links->lat->coord_to_index((int*)coords, c)];
                    }
                }
                break;
            case DIR_XP1: // +x
                for (i = 0; i < lattice_size; i++)
                {
                    lhs[i] = 0.0; // zero out.
                    links->lat->index_to_coord(i, (int*)coords, color); // get the coordinate and color index.
                    coords_tmp[0] = (coords[0] + 1)%links->lat->get_lattice_dimension(0);
                    coords_tmp[1] = coords[1];
                    for (c = 0; c < nc; c++)
                    {
                        lhs[i] += links->hopping[c + nc*i + 0*nc*lattice_size]*rhs[links->lat->coord_to_index((int*)coords_tmp, c)];
                    }
                }
                break;
            case DIR_YP1: // +y
                for (i = 0; i < lattice_size; i++)
                {
                    lhs[i] = 0.0; // zero out.
                    links->lat->index_to_coord(i, (int*)coords, color); // get the coordinate and color index.
                    coords_tmp[0] = coords[0];
                    coords_tmp[1] = (coords[1] + 1)%links->lat->get_lattice_dimension(1);
                    for (c = 0; c < nc; c++)
                    {
                        lhs[i] += links->hopping[c + nc*i + 1*nc*lattice_size]*rhs[links->lat->coord_to_index((int*)coords_tmp, c)];
                    }
                }
                break;
            case DIR_XM1: // -x
                for (i = 0; i < lattice_size; i++)
                {
                    lhs[i] = 0.0; // zero out.
                    links->lat->index_to_coord(i, (int*)coords, color); // get the coordinate and color index.
                    coords_tmp[0] = (coords[0] - 1 + links->lat->get_lattice_dimension(0))%links->lat->get_lattice_dimension(0);
                    coords_tmp[1] = coords[1];
                    for (c = 0; c < nc; c++)
                    {
                        lhs[i] += links->hopping[c + nc*i + 2*nc*lattice_size]*rhs[links->lat->coord_to_index((int*)coords_tmp, c)];
                    }
                }
                break;
            case DIR_YM1: //-y
                for (i = 0; i < lattice_size; i++)
                {
                    lhs[i] = 0.0; // zero out.
                    links->lat->index_to_coord(i, (int*)coords, color); // get the coordinate and color index.
                    coords_tmp[0] = coords[0];
                    coords_tmp[1] = (coords[1] - 1 + links->lat->get_lattice_dimension(1))%links->lat->get_lattice_dimension(1);
                    for (c = 0; c < nc; c++)
                    {
                        lhs[i] += links->hopping[c + nc*i + 3*nc*lattice_size]*rhs[links->lat->coord_to_index((int*)coords_tmp, c)];
                    }
                }
                break;
            case DIR_XP2: // +2x
                if (links->has_two)
                {
                    for (i = 0; i < lattice_size; i++)
                    {
                        lhs[i] = 0.0; // zero out.
                        links->lat->index_to_coord(i, (int*)coords, color); // get the coordinate and color index.
                        coords_tmp[0] = (coords[0] + 2)%links->lat->get_lattice_dimension(0);
                        coords_tmp[1] = coords[1];
                        for (c = 0; c < nc; c++)
                        {
                            lhs[i] += links->two_link[c + nc*i+0*nc*lattice_size]*rhs[links->lat->coord_to_index((int*)coords_tmp, c)];
                        }
                    }
                }
                break;
            case DIR_XP1YP1: // +x+y
                if (links->has_two)
                {
                    for (i = 0; i < lattice_size; i++)
                    {
                        lhs[i] = 0.0; // zero out.
                        links->lat->index_to_coord(i, (int*)coords, color); // get the coordinate and color index.
                        coords_tmp[0] = (coords[0] + 1)%links->lat->get_lattice_dimension(0);
                        coords_tmp[1] = (coords[1] + 1)%links->lat->get_lattice_dimension(1);
                        for (c = 0; c < nc; c++)
                        {
                            lhs[i] += links->two_link[c + nc*i+1*nc*lattice_size]*rhs[links->lat->coord_to_index((int*)coords_tmp, c)];
                        }
                    }
                }
                break;
            case DIR_YP2: // +2y
                if (links->has_two)
                {
                    for (i = 0; i < lattice_size; i++)
                    {
                        lhs[i] = 0.0; // zero out.
                        links->lat->index_to_coord(i, (int*)coords, color); // get the coordinate and color index.
                        coords_tmp[0] = coords[0];
                        coords_tmp[1] = (coords[1] + 2)%links->lat->get_lattice_dimension(1);
                        for (c = 0; c < nc; c++)
                        {
                            lhs[i] += links->two_link[c + nc*i+2*nc*lattice_size]*rhs[links->lat->coord_to_index((int*)coords_tmp, c)];
                        }
                    }
                }
                break;
            case DIR_XM1YP1: // -x+y
                if (links->has_two)
                {
                    for (i = 0; i < lattice_size; i++)
                    {
                        lhs[i] = 0.0; // zero out.
                        links->lat->index_to_coord(i, (int*)coords, color); // get the coordinate and color index.
                        coords_tmp[0] = (coords[0] - 1 + 2*links->lat->get_lattice_dimension(0))%links->lat->get_lattice_dimension(0);
                        coords_tmp[1] = (coords[1] + 1)%links->lat->get_lattice_dimension(1);
                        for (c = 0; c < nc; c++)
                        {
                            lhs[i] += links->two_link[c + nc*i+3*nc*lattice_size]*rhs[links->lat->coord_to_index((int*)coords_tmp, c)];
                        }

                    }
                }
                break;
            case DIR_XM2: // -2x
                if (links->has_two)
                {
                    for (i = 0; i < lattice_size; i++)
                    {
                        lhs[i] = 0.0; // zero out.
                        links->lat->index_to_coord(i, (int*)coords, color); // get the coordinate and color index.
                        coords_tmp[0] = (coords[0] - 2 + 2*links->lat->get_lattice_dimension(0))%links->lat->get_lattice_dimension(0);
                        coords_tmp[1] = coords[1];
                        for (c = 0; c < nc; c++)
                        {
                            lhs[i] += links->two_link[c + nc*i+4*nc*lattice_size]*rhs[links->lat->coord_to_index((int*)coords_tmp, c)];
                        }
                    }
                }
                break;
            case DIR_XM1YM1: // -x-y
                if (links->has_two)
                {
                    for (i = 0; i < lattice_size; i++)
                    {
                        lhs[i] = 0.0; // zero out.
                        links->lat->index_to_coord(i, (int*)coords, color); // get the coordinate and color index.
                        coords_tmp[0] = (coords[0] - 1 + 2*links->lat->get_lattice_dimension(0))%links->lat->get_lattice_dimension(0);
                        coords_tmp[1] = (coords[1] - 1 + 2*links->lat->get_lattice_dimension(1))%links->lat->get_lattice_dimension(1);
                        for (c = 0; c < nc; c++)
                        {
                            lhs[i] += links->two_link[c + nc*i+5*nc*lattice_size]*rhs[links->lat->coord_to_index((int*)coords_tmp, c)];
                        }
                    }
                }
                break;
            case DIR_YM2: // -2y
                if (links->has_two)
                {
                    for (i = 0; i < lattice_size; i++)
                    {
                        lhs[i] = 0.0; // zero out.
                        links->lat->index_to_coord(i, (int*)coords, color); // get the coordinate and color index.
                        coords_tmp[0] = coords[0];
                        coords_tmp[1] = (coords[1] - 2 + 2*links->lat->get_lattice_dimension(1))%links->lat->get_lattice_dimension(1);
                        for (c = 0; c < nc; c++)
                        {
                            lhs[i] += links->two_link[c + nc*i+6*nc*lattice_size]*rhs[links->lat->coord_to_index((int*)coords_tmp, c)];
                        }
                    }
                }
                break;
            case DIR_XP1YM1: // +x-y
                if (links->has_two)
                {
                    for (i = 0; i < lattice_size; i++)
                    {
                        lhs[i] = 0.0; // zero out.
                        links->lat->index_to_coord(i, (int*)coords, color); // get the coordinate and color index.
                        coords_tmp[0] = (coords[0] + 1)%links->lat->get_lattice_dimension(0);
                        coords_tmp[1] = (coords[1] - 1 + 2*links->lat->get_lattice_dimension(1))%links->lat->get_lattice_dimension(1);
                        for (c = 0; c < nc; c++)
                        {
                            lhs[i] += links->two_link[c + nc*i+7*nc*lattice_size]*rhs[links->lat->coord_to_index((int*)coords_tmp, c)];
                        }
                    }
                }
                break;
        } // end switch
    } // end only one link for the stencil. 
}


// Generate a stencil operator given a function and a coarsening distance (max 2...). Assumes that the physical lengths are longer than coarsening distance*2. 
void generate_stencil_2d(stencil_2d* stenc, void (*matrix_vector)(complex<double>*,complex<double>*,void*), void* extra_data)
{
    
    if (stenc->generated || stenc->stencil_size > 2)
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



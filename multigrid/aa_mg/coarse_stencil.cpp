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
        
        // And let's stop there for now. 
    }
}
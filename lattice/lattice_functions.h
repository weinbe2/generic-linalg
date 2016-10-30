// ESW Sun Oct 30 12:54:24 EDT 2016
// A few functions useful for LQCD calculations.

#ifndef LATTICE_FCN_HEADER
#define LATTICE_FCN_HEADER

#include "lattice.h"

// Apply (-1)^[sum of coordinates], \epsilon in the context of staggered fermions.
// Safe for out and in to be the same vector.
void lattice_epsilon(complex<double>* out, complex<double>* in, Lattice* latt)
{
  int len = latt->get_lattice_size();
  for (int i = 0; i < len; i++)
  {
    if (!latt->index_is_even(i))
    {
      out[i] = -in[i];
    }
    else
    {
      out[i] = in[i];
    }
  }
}

// Apply \sigma_3 if the number of colors is even. This multiplies the bottom half of the dof by -1.
// Safe for out and in to be the same vector.
void lattice_sigma3(complex<double>* out, complex<double>* in, Lattice* latt)
{
  int nc = latt->get_nc();
  int vol = latt->get_volume(); // coordinate volume
  if (nc % 2 != 0)
  {
    for (int i = 0; i < vol*nc; i++)
    {
      out[i] = in[i];
    }
    return; // probably shouldn't silent error...
  }
  for (int i = 0; i < vol; i++)
  {
    for (int c = 0; c < nc/2; c++)
    {
      out[i*nc+c] = in[i*nc+c];
    }
    
    for (int c = nc/2; c < nc; c++)
    {
      out[i*nc+c] = -in[i*nc+c];
    }
  }
}

#endif
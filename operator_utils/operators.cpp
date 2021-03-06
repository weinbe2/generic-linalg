
#include <iostream>
#include <complex>
#include "operators.h"

using namespace std;

int get_stencil_size(op_type opt)
{
  switch (opt)
  {
    case STAGGERED:
    case LAPLACE:
    case LAPLACE_NC2:
    case G5_STAGGERED:
      return 1;
      break;
    case STAGGERED_INDEX: // Has around the corner
    case STAGGERED_NORMAL: // Has around the corner, two link.
      return 2;
  }
  return 0;
}

// Square lattice.
// Kinetic term for a 2D laplace w/ period bc. Applies lhs = A*rhs.
// The unit vectors are e_1 = xhat, e_2 = yhat.
void square_laplace(complex<double>* lhs, complex<double>* rhs, void* extra_data)
{
    // Declare variables.
    int i;
    int x,y,c;
    int tmp; 
    staggered_u1_op* stagif = (staggered_u1_op*)extra_data;
    double mass = stagif->mass; 
    int x_fine = stagif->x_fine;
    int y_fine = stagif->y_fine;
    int Nc = stagif->Nc; // 1 is the trivial laplace. 
    
    // Apply the stencil.
    for (i = 0; i < x_fine*y_fine*Nc; i++)
    {
        lhs[i] = 0.0;
        c = i%Nc;      // get color.
        tmp = (i-c)/Nc;
        x = tmp%x_fine; // integer mod.
        y = tmp/x_fine; // integer divide.

        // + e1.
        lhs[i] = lhs[i]-rhs[y*x_fine*Nc+((x+1)%x_fine)*Nc+c];
        
        // - e1.
        lhs[i] = lhs[i]- rhs[y*x_fine*Nc+((x+x_fine-1)%x_fine)*Nc+c]; // The extra +N is because of the % sign convention.

        // + e2.
        lhs[i] = lhs[i]- rhs[((y+1)%y_fine)*x_fine*Nc+x*Nc+c];

        // - e2.
        lhs[i] = lhs[i]- rhs[((y+y_fine-1)%y_fine)*x_fine*Nc+x*Nc+c];

        // 0
        // Added mass term here.
        lhs[i] = lhs[i]+ (4+mass)*rhs[i];
    }

}

// Square lattice.
// Kinetic term for a 2D laplacian w/ period bc. Applies lhs = A*rhs.
// The unit vectors are e_1 = xhat, e_2 = yhat.
// The "extra_data" is a cast to a complex gauge_field[N*N*2], 
//    loaded by the function read_lattice_u1. 
void square_laplace_u1(complex<double>* lhs, complex<double>* rhs, void* extra_data)
{
    // Declare variables.
    int i;
    int x,y;

    staggered_u1_op* stagif = (staggered_u1_op*)extra_data;
    complex<double>* lattice = stagif->lattice;
    double mass = stagif->mass; 
    int x_fine = stagif->x_fine;
    int y_fine = stagif->y_fine;

    // For a 2D square lattice, the stencil is:
    //     |  0 -1  0 |
    //     | -1 +4 -1 |
    //     |  0 -1  0 |
    //
    // e2 = yhat
    // ^
    // | 
    // |-> e1 = xhat
    
    // Apply the stencil.
    for (i = 0; i < x_fine*y_fine; i++)
    {
        lhs[i] = 0.0;
        x = i%x_fine; // integer mod.
        y = i/x_fine; // integer divide.

        // + e1.
        lhs[i] = lhs[i]-lattice[y*x_fine*2+x*2]*rhs[y*x_fine+((x+1)%x_fine)];

        // - e1.
        lhs[i] = lhs[i]-conj(lattice[y*x_fine*2+((x+x_fine-1)%x_fine)*2])*rhs[y*x_fine+((x+x_fine-1)%x_fine)]; // The extra +N is because of the % sign convention.

        // + e2.
        lhs[i] = lhs[i]-lattice[y*x_fine*2+x*2+1]*rhs[((y+1)%y_fine)*x_fine+x];

        // - e2.
        lhs[i] = lhs[i]-conj(lattice[((y+y_fine-1)%y_fine)*x_fine*2+x*2+1])*rhs[((y+y_fine-1)%y_fine)*x_fine+x];

        // 0
        // Added mass term here.
        lhs[i] = lhs[i]+ (4+mass)*rhs[i];

    }

}

// Square lattice.
// Kinetic term for a 2D staggered w/ period bc. Applies lhs = A*rhs.
// The unit vectors are e_1 = xhat, e_2 = yhat.
// The "extra_data" doesn't include anything.
// Apsi[x][y] = m psi[x,y] - U[y][x,x+1] 
void square_staggered(complex<double>* lhs, complex<double>* rhs, void* extra_data)
{
    // Declare variables.
    int i;
    int x,y;
    double eta1;
    staggered_u1_op* stagif = (staggered_u1_op*)extra_data;
    double mass = stagif->mass; 
    int x_fine = stagif->x_fine;
    int y_fine = stagif->y_fine;

    // For a 2D square lattice, the stencil is:
    //   1 |  0 -eta1  0 |
    //   - | +1    0  -1 |  , where eta1 = (-1)^x
    //   2 |  0 +eta1  0 |
    //
    // e2 = yhat
    // ^
    // | 
    // |-> e1 = xhat

    // Apply the stencil.
    for (i = 0; i < x_fine*y_fine; i++)
    {
        lhs[i] = 0.0;
        x = i%x_fine; // integer mod.
        y = i/x_fine; // integer divide.
        eta1 = 1 - 2*(x%2);

        // + e1.
        lhs[i] = lhs[i]-rhs[y*x_fine+((x+1)%x_fine)];

        // - e1.
        lhs[i] = lhs[i]+ rhs[y*x_fine+((x+x_fine-1)%x_fine)]; // The extra +N is because of the % sign convention.

        // + e2.
        lhs[i] = lhs[i]- eta1*rhs[((y+1)%y_fine)*x_fine+x];

        // - e2.
        lhs[i] = lhs[i]+ eta1*rhs[((y+y_fine-1)%y_fine)*x_fine+x];

        // Normalization.
        lhs[i] = 0.5*lhs[i];

        // 0
        // Added mass term here.
        lhs[i] = lhs[i]+ mass*rhs[i];
    }

}

// Square lattice.
// Kinetic term for a 2D staggered w/ period bc. Applies lhs = A*rhs.
// The unit vectors are e_1 = xhat, e_2 = yhat.
// The "extra_data" is a cast to a complex gauge_field[N*N*2], 
//    loaded by the function read_lattice_u1. 
// Apsi[x][y] = m psi[x,y] - U[y][x,x+1] 
void square_staggered_u1(complex<double>* lhs, complex<double>* rhs, void* extra_data)
{
    // Declare variables.
    int i;
    int x,y;
    double eta1; 
    staggered_u1_op* stagif = (staggered_u1_op*)extra_data;
    complex<double>* lattice = stagif->lattice;
    double mass = stagif->mass; 
    int x_fine = stagif->x_fine;
    int y_fine = stagif->y_fine;

    // For a 2D square lattice, the stencil is:
    //   1 |  0 -eta1  0 |
    //   - | +1    0  -1 |  , where eta1 = (-1)^x
    //   2 |  0 +eta1  0 |
    //
    // e2 = yhat
    // ^
    // | 
    // |-> e1 = xhat

    // Apply the stencil.
    for (i = 0; i < x_fine*y_fine; i++)
    {
        lhs[i] = 0.0;
        x = i%x_fine; // integer mod.
        y = i/x_fine; // integer divide.
        eta1 = 1 - 2*(x%2);

        // + e1.
        lhs[i] = lhs[i]-lattice[y*x_fine*2+x*2]*rhs[y*x_fine+((x+1)%x_fine)];

        // - e1.
        lhs[i] = lhs[i]+ conj(lattice[y*x_fine*2+((x+x_fine-1)%x_fine)*2])*rhs[y*x_fine+((x+x_fine-1)%x_fine)]; // The extra +N is because of the % sign convention.

        // + e2.
        lhs[i] = lhs[i]- eta1*lattice[y*x_fine*2+x*2+1]*rhs[((y+1)%y_fine)*x_fine+x];

        // - e2.
        lhs[i] = lhs[i]+ eta1*conj(lattice[((y+y_fine-1)%y_fine)*x_fine*2+x*2+1])*rhs[((y+y_fine-1)%y_fine)*x_fine+x];

        // Normalization.
        lhs[i] = 0.5*lhs[i];

        // 0
        // Added mass term here.
        lhs[i] = lhs[i]+ mass*rhs[i];

        // Apply a gamma5.
        /*if ((x+y)%2 == 1)
        {
            lhs[i] = -lhs[i];
        }*/
    }
}

// \gamma_5
void gamma_5(complex<double>* lhs, complex<double>* rhs, void* extra_data)
{
    // Declare variables.
    int i;
    int x,y;
    staggered_u1_op* stagif = (staggered_u1_op*)extra_data;
    int x_fine = stagif->x_fine;
    int y_fine = stagif->y_fine;

    // Apply the stencil.
    for (i = 0; i < x_fine*y_fine; i++)
    {
        x = i%x_fine; // integer mod.
        y = i/x_fine; // integer divide.
        
        lhs[i] = ((double)(1 - 2*((x+y)%2)))*rhs[i];
    }
}

// Square \gamma_5 staggered 2d operator w/out u1 function.
void square_staggered_gamma5(complex<double>* lhs, complex<double>* rhs, void* extra_data)
{
    // Declare variables.
    int i;
    int x,y;
    double eta1; 
    staggered_u1_op* stagif = (staggered_u1_op*)extra_data;
    double mass = stagif->mass; 
    int x_fine = stagif->x_fine;
    int y_fine = stagif->y_fine;
    double eo_sign; 

    // Apply the stencil.
    for (i = 0; i < x_fine*y_fine; i++)
    {
        lhs[i] = 0.0;
        x = i%x_fine; // integer mod.
        y = i/x_fine; // integer divide.
        eta1 = 1 - 2*(x%2);
        eo_sign = ((x+y)%2 == 0) ? 1.0 : -1.0;

        // + e1.
        lhs[i] = lhs[i]-eo_sign*rhs[y*x_fine+((x+1)%x_fine)];

        // - e1.
        lhs[i] = lhs[i]+ eo_sign*rhs[y*x_fine+((x+x_fine-1)%x_fine)]; // The extra +N is because of the % sign convention.

        // + e2.
        lhs[i] = lhs[i]- eo_sign*eta1*rhs[((y+1)%y_fine)*x_fine+x];

        // - e2.
        lhs[i] = lhs[i]+ eo_sign*eta1*rhs[((y+y_fine-1)%y_fine)*x_fine+x];

        // Normalization.
        lhs[i] = 0.5*lhs[i];

        // 0
        // Added mass term here.
        lhs[i] = lhs[i]+ eo_sign*mass*rhs[i];

    }
    
    /*    
    staggered_u1_op* stagif = (staggered_u1_op*)extra_data;
    complex<double>* tmp = new complex<double>[stagif->x_fine*stagif->y_fine];
    
    square_staggered(tmp, rhs, extra_data);
    gamma_5(lhs, tmp, extra_data);
    
    delete[] tmp;
    */
}

// Square \gamma_5 staggered 2d operator w/ u1 function.
void square_staggered_gamma5_u1(complex<double>* lhs, complex<double>* rhs, void* extra_data)
{
    // Declare variables.
    int i;
    int x,y;
    double eta1; 
    staggered_u1_op* stagif = (staggered_u1_op*)extra_data;
    complex<double>* lattice = stagif->lattice;
    double mass = stagif->mass; 
    int x_fine = stagif->x_fine;
    int y_fine = stagif->y_fine;
    double eo_sign; 

    // Apply the stencil.
    for (i = 0; i < x_fine*y_fine; i++)
    {
        lhs[i] = 0.0;
        x = i%x_fine; // integer mod.
        y = i/x_fine; // integer divide.
        eta1 = 1 - 2*(x%2);
        eo_sign = ((x+y)%2 == 0) ? 1.0 : -1.0;

        // + e1.
        lhs[i] = lhs[i]-eo_sign*lattice[y*x_fine*2+x*2]*rhs[y*x_fine+((x+1)%x_fine)];

        // - e1.
        lhs[i] = lhs[i]+ eo_sign*conj(lattice[y*x_fine*2+((x+x_fine-1)%x_fine)*2])*rhs[y*x_fine+((x+x_fine-1)%x_fine)]; // The extra +N is because of the % sign convention.

        // + e2.
        lhs[i] = lhs[i]- eo_sign*eta1*lattice[y*x_fine*2+x*2+1]*rhs[((y+1)%y_fine)*x_fine+x];

        // - e2.
        lhs[i] = lhs[i]+ eo_sign*eta1*conj(lattice[((y+y_fine-1)%y_fine)*x_fine*2+x*2+1])*rhs[((y+y_fine-1)%y_fine)*x_fine+x];

        // Normalization.
        lhs[i] = 0.5*lhs[i];

        // 0
        // Added mass term here.
        lhs[i] = lhs[i]+ eo_sign*mass*rhs[i];

    }
    
    /*
    staggered_u1_op* stagif = (staggered_u1_op*)extra_data;
    complex<double>* tmp = new complex<double>[stagif->x_fine*stagif->y_fine];
    
    square_staggered_u1(tmp, rhs, extra_data);
    gamma_5(lhs, tmp, extra_data);
    
    delete[] tmp;*/
    
}

// staggered dagger operator.
// Mass is the same, hopping is flipped.
void square_staggered_dagger_u1(complex<double>* lhs, complex<double>* rhs, void* extra_data)
{
    // Declare variables.
    int i;
    int x,y;
    double eta1; 
    staggered_u1_op* stagif = (staggered_u1_op*)extra_data;
    complex<double>* lattice = stagif->lattice;
    double mass = stagif->mass; 
    int x_fine = stagif->x_fine;
    int y_fine = stagif->y_fine;

    // For a 2D square lattice, the stencil is:
    //   1 |  0 +eta1  0 |
    //   - | -1    0  +1 |  , where eta1 = (-1)^x
    //   2 |  0 -eta1  0 |
    //
    // e2 = yhat
    // ^
    // | 
    // |-> e1 = xhat

    // Apply the stencil.
    for (i = 0; i < x_fine*y_fine; i++)
    {
        lhs[i] = 0.0;
        x = i%x_fine; // integer mod.
        y = i/x_fine; // integer divide.
        eta1 = 1 - 2*(x%2);

        // + e1.
        lhs[i] = lhs[i]+lattice[y*x_fine*2+x*2]*rhs[y*x_fine+((x+1)%x_fine)];

        // - e1.
        lhs[i] = lhs[i]- conj(lattice[y*x_fine*2+((x+x_fine-1)%x_fine)*2])*rhs[y*x_fine+((x+x_fine-1)%x_fine)]; // The extra +N is because of the % sign convention.

        // + e2.
        lhs[i] = lhs[i]+ eta1*lattice[y*x_fine*2+x*2+1]*rhs[((y+1)%y_fine)*x_fine+x];

        // - e2.
        lhs[i] = lhs[i]- eta1*conj(lattice[((y+y_fine-1)%y_fine)*x_fine*2+x*2+1])*rhs[((y+y_fine-1)%y_fine)*x_fine+x];

        // Normalization.
        lhs[i] = 0.5*lhs[i];

        // 0
        // Added mass term here.
        lhs[i] = lhs[i]+ mass*rhs[i];

        // Apply a gamma5.
        /*if ((x+y)%2 == 1)
        {
            lhs[i] = -lhs[i];
        }*/
    }
    
    /*
    staggered_u1_op* stagif = (staggered_u1_op*)extra_data;
    complex<double>* tmp = new complex<double>[stagif->x_fine*stagif->y_fine];
    
    gamma_5(tmp, rhs, extra_data);
    square_staggered_u1(lhs, tmp, extra_data);
    gamma_5(tmp, lhs, extra_data);
    for (int i = 0; i < stagif->x_fine*stagif->y_fine; i++)
    {
        lhs[i] = tmp[i];
    }
    
    delete[] tmp;*/
}

// Staggered normal equations.
void square_staggered_normal_u1(complex<double>* lhs, complex<double>* rhs, void* extra_data)
{
    staggered_u1_op* stagif = (staggered_u1_op*)extra_data;
    complex<double>* tmp = new complex<double>[stagif->x_fine*stagif->y_fine];
    
    square_staggered_u1(tmp, rhs, extra_data);
    square_staggered_dagger_u1(lhs, tmp, extra_data);
    
    delete[] tmp;
}

// Square staggered 2d operator w/ u1 function, D_{eo} only.
void square_staggered_deo_u1(complex<double>* lhs, complex<double>* rhs, void* extra_data)
{
    // Declare variables.
    int i;
    int x,y;
    double eta1; 
    staggered_u1_op* stagif = (staggered_u1_op*)extra_data;
    complex<double>* lattice = stagif->lattice;
    int x_fine = stagif->x_fine;
    int y_fine = stagif->y_fine;

    // Apply the stencil.
    for (i = 0; i < x_fine*y_fine; i++)
    {
        lhs[i] = 0.0;
        x = i%x_fine; // integer mod.
        y = i/x_fine; // integer divide.
        eta1 = 1 - 2*(x%2);
        
        if ((x+y)%2 == 1) // update even only.
        {
            continue; 
        }

        // + e1.
        lhs[i] = lhs[i]-lattice[y*x_fine*2+x*2]*rhs[y*x_fine+((x+1)%x_fine)];

        // - e1.
        lhs[i] = lhs[i]+ conj(lattice[y*x_fine*2+((x+x_fine-1)%x_fine)*2])*rhs[y*x_fine+((x+x_fine-1)%x_fine)]; // The extra +N is because of the % sign convention.
      
        // + e2.
        lhs[i] = lhs[i]- eta1*lattice[y*x_fine*2+x*2+1]*rhs[((y+1)%y_fine)*x_fine+x];

        // - e2.
        lhs[i] = lhs[i]+ eta1*conj(lattice[((y+y_fine-1)%y_fine)*x_fine*2+x*2+1])*rhs[((y+y_fine-1)%y_fine)*x_fine+x];
      
        // Normalization.
        lhs[i] = 0.5*lhs[i];
    }
    
}

// Square staggered 2d operator w/ u1 function, D_{oe} only.
void square_staggered_doe_u1(complex<double>* lhs, complex<double>* rhs, void* extra_data)
{
    // Declare variables.
    int i;
    int x,y;
    double eta1; 
    staggered_u1_op* stagif = (staggered_u1_op*)extra_data;
    complex<double>* lattice = stagif->lattice;
    int x_fine = stagif->x_fine;
    int y_fine = stagif->y_fine;

    // Apply the stencil.
    for (i = 0; i < x_fine*y_fine; i++)
    {
        lhs[i] = 0.0;
        x = i%x_fine; // integer mod.
        y = i/x_fine; // integer divide.
        eta1 = 1 - 2*(x%2);
        
        if ((x+y)%2 == 0) // update odd only.
        {
            continue; 
        }

        // + e1.
        lhs[i] = lhs[i]-lattice[y*x_fine*2+x*2]*rhs[y*x_fine+((x+1)%x_fine)];

        // - e1.
        lhs[i] = lhs[i]+ conj(lattice[y*x_fine*2+((x+x_fine-1)%x_fine)*2])*rhs[y*x_fine+((x+x_fine-1)%x_fine)]; // The extra +N is because of the % sign convention.

        // + e2.
        lhs[i] = lhs[i]- eta1*lattice[y*x_fine*2+x*2+1]*rhs[((y+1)%y_fine)*x_fine+x];

        // - e2.
        lhs[i] = lhs[i]+ eta1*conj(lattice[((y+y_fine-1)%y_fine)*x_fine*2+x*2+1])*rhs[((y+y_fine-1)%y_fine)*x_fine+x];

        // Normalization.
        lhs[i] = 0.5*lhs[i];

    }
}

// Prepare an even rhs for an even/odd preconditioned solve.
// Takes in rhs_orig, returns rhs_e. 
void square_staggered_eoprec_prepare(complex<double>* rhs_e, complex<double>* rhs_orig, void* extra_data)
{
  int i,x,y;
  staggered_u1_op* stagif = (staggered_u1_op*)extra_data;
  int x_fine = stagif->x_fine;
  int y_fine = stagif->y_fine;
  int fine_size = x_fine*y_fine; // *Nc
  
  // Prepare rhs: m rhs_e - D_{eo} rhs_o
  square_staggered_deo_u1(rhs_e, rhs_orig, extra_data); // zeroes odd sites in rhs_e.
  for (i = 0; i < fine_size; i++)
  {
    x = i%x_fine;
    y = i/x_fine;
    
    if ((x+y)%2 == 0) // even
    {
      rhs_e[i] = stagif->mass*rhs_orig[i] - rhs_e[i];
    }
  }
}

// Square staggered 2d operator w/ u1 function, m^2 - D_{eo} D_{oe} [zeroes odd explicitly]
void square_staggered_m2mdeodoe_u1(complex<double>* lhs, complex<double>* rhs, void* extra_data)
{
  int i,x,y;
  staggered_u1_op* stagif = (staggered_u1_op*)extra_data;
  complex<double>* tmp = new complex<double>[stagif->x_fine*stagif->y_fine];
  int x_fine = stagif->x_fine;
  int y_fine = stagif->y_fine;

  square_staggered_doe_u1(tmp, rhs, extra_data);
  square_staggered_deo_u1(lhs, tmp, extra_data);

  for (i = 0; i < x_fine*y_fine; i++)
  {
    x = i%x_fine; // integer mod.
    y = i/x_fine; // integer divide.
    if ((x+y)%2 == 0)
    {
      lhs[i] = stagif->mass*stagif->mass*rhs[i] - lhs[i];
    }
  }

  delete[] tmp;
}

// Reconstruct the full solution for an even/odd preconditioned solve.
// Takes in lhs_e, rhs_o, returns lhs_full (copying over the even part from lhs_e)
void square_staggered_eoprec_reconstruct(complex<double>* lhs_full, complex<double>* lhs_e, complex<double>* rhs_o, void* extra_data)
{
  int i,x,y;
  staggered_u1_op* stagif = (staggered_u1_op*)extra_data;
  int x_fine = stagif->x_fine;
  int y_fine = stagif->y_fine;
  int fine_size = x_fine*y_fine; // *Nc
  double inv_mass = 1.0/stagif->mass;
  
  // Reconstruct odd: m^{-1}*(rhs_o - D_{oe} lhs2_e)
  square_staggered_doe_u1(lhs_full, lhs_e, extra_data);
  for (i = 0; i < fine_size; i++)
  {
    x = i%x_fine;
    y = i/x_fine;
    if ((x+y)%2 == 1) // odd
    {
      lhs_full[i] = inv_mass*(rhs_o[i] - lhs_full[i]);
    }
    else
    {
      lhs_full[i] = lhs_e[i];
    }
  }
}


// Square lattice.
// Kinetic term for a 2D staggered w/ period bc, plus 2-link laplace.
// Applies lhs = A*rhs.
// The unit vectors are e_1 = xhat, e_2 = yhat.
// The "extra_data" is a cast to a complex gauge_field[N*N*2], 
//    loaded by the function read_lattice_u1. 
// Apsi[x][y] = m psi[x,y] - U[y][x,x+1] 
void square_staggered_2linklaplace_u1(complex<double>* lhs, complex<double>* rhs, void* extra_data)
{
  // Declare variables.
  int i;
  int x,y;
  double eta1; 
  staggered_u1_op* stagif = (staggered_u1_op*)extra_data;
  complex<double>* lattice = stagif->lattice;
  double mass = stagif->mass; 
  int x_fine = stagif->x_fine;
  int y_fine = stagif->y_fine;

  // For a 2D square lattice, the stencil is:
  //   1 |  0 -eta1  0 |
  //   - | +1    0  -1 |  , where eta1 = (-1)^x
  //   2 |  0 +eta1  0 |
  //
  // e2 = yhat
  // ^
  // | 
  // |-> e1 = xhat
  // 
  // Plus a two link laplace stencil term. 

  // Apply the stencil.
  for (i = 0; i < x_fine*y_fine; i++)
  {
    lhs[i] = 0.0;
    x = i%x_fine; // integer mod.
    y = i/x_fine; // integer divide.
    eta1 = 1 - 2*(x%2);

    // + e1.
    lhs[i] = lhs[i]-lattice[y*x_fine*2+x*2]*rhs[y*x_fine+((x+1)%x_fine)];

    // - e1.
    lhs[i] = lhs[i]+ conj(lattice[y*x_fine*2+((x+x_fine-1)%x_fine)*2])*rhs[y*x_fine+((x+x_fine-1)%x_fine)]; // The extra +N is because of the % sign convention.

    // + e2.
    lhs[i] = lhs[i]- eta1*lattice[y*x_fine*2+x*2+1]*rhs[((y+1)%y_fine)*x_fine+x];

    // - e2.
    lhs[i] = lhs[i]+ eta1*conj(lattice[((y+y_fine-1)%y_fine)*x_fine*2+x*2+1])*rhs[((y+y_fine-1)%y_fine)*x_fine+x];

    // Normalization.
    lhs[i] = 0.5*lhs[i];

    // 0
    // Added mass term here.
    lhs[i] = lhs[i]+ mass*rhs[i];
    
    // Add the two link laplace term, rescaled by a wilson coeff. 
    
    // 0
    lhs[i] = lhs[i] + stagif->wilson_coeff * 4.0*rhs[i];
    
    // + 2*e1.
    lhs[i] = lhs[i]-stagif->wilson_coeff * lattice[y*x_fine*2+x*2] * lattice[y*x_fine*2+((x+1)%x_fine)*2] * rhs[y*x_fine+((x+2)%x_fine)];

    // - 2*e1.
    lhs[i] = lhs[i]-stagif->wilson_coeff * conj(lattice[y*x_fine*2+((x+x_fine-1)%x_fine)*2]) * conj(lattice[y*x_fine*2+((x+x_fine-2)%x_fine)*2]) * rhs[y*x_fine+((x+x_fine-2)%x_fine)]; 

    // + 2*e2.
    lhs[i] = lhs[i]-stagif->wilson_coeff * lattice[y*x_fine*2+x*2+1] * lattice[((y+1)%y_fine)*x_fine*2+x*2+1] * rhs[((y+2)%y_fine)*x_fine+x];

    // - 2*e2.
    lhs[i] = lhs[i]-stagif->wilson_coeff * conj(lattice[((y+y_fine-1)%y_fine)*x_fine*2+x*2+1]) * conj(lattice[((y+y_fine-2)%y_fine)*x_fine*2+x*2+1]) * rhs[((y+y_fine-2)%y_fine)*x_fine+x];
  }
}

// Operators for symmetric shifts.
void staggered_symmshift_x(complex<double>* lhs, complex<double>* rhs, void* extra_data)
{
    // Declare variables.
    int i;
    int x,y;
    staggered_u1_op* stagif = (staggered_u1_op*)extra_data;
    complex<double>* lattice = stagif->lattice;
    int x_fine = stagif->x_fine;
    int y_fine = stagif->y_fine;

    // For a 2D square lattice, the stencil is:
    //   1 |  0    0   0 |
    //   - | +1    0  +1 |  , where eta1 = (-1)^x
    //   2 |  0    0   0 |
    //
    // e2 = yhat
    // ^
    // | 
    // |-> e1 = xhat

    // Apply the stencil.
    for (i = 0; i < x_fine*y_fine; i++)
    {
        lhs[i] = 0.0;
        x = i%x_fine; // integer mod.
        y = i/x_fine; // integer divide.

        // + e1.
        lhs[i] = lhs[i]+lattice[y*x_fine*2+x*2]*rhs[y*x_fine+((x+1)%x_fine)];

        // - e1.
        lhs[i] = lhs[i]+ conj(lattice[y*x_fine*2+((x+x_fine-1)%x_fine)*2])*rhs[y*x_fine+((x+x_fine-1)%x_fine)]; // The extra +N is because of the % sign convention.


        // Normalization.
        lhs[i] = 0.5*lhs[i];


    }
}


void staggered_symmshift_y(complex<double>* lhs, complex<double>* rhs, void* extra_data)
{
    // Declare variables.
    int i;
    int x,y;
    double eta1; 
    staggered_u1_op* stagif = (staggered_u1_op*)extra_data;
    complex<double>* lattice = stagif->lattice;
    int x_fine = stagif->x_fine;
    int y_fine = stagif->y_fine;

    // For a 2D square lattice, the stencil is:
    //   1 |  0 +eta1  0 |
    //   - |  0    0   0 |  , where eta1 = (-1)^x
    //   2 |  0 +eta1  0 |
    //
    // e2 = yhat
    // ^
    // | 
    // |-> e1 = xhat

    // Apply the stencil.
    for (i = 0; i < x_fine*y_fine; i++)
    {
        lhs[i] = 0.0;
        x = i%x_fine; // integer mod.
        y = i/x_fine; // integer divide.
        eta1 = 1 - 2*(x%2);

        // + e2.
        lhs[i] = lhs[i]+ eta1*lattice[y*x_fine*2+x*2+1]*rhs[((y+1)%y_fine)*x_fine+x];

        // - e2.
        lhs[i] = lhs[i]+ eta1*conj(lattice[((y+y_fine-1)%y_fine)*x_fine*2+x*2+1])*rhs[((y+y_fine-1)%y_fine)*x_fine+x];

        // Normalization.
        lhs[i] = 0.5*lhs[i];

    }
}

// Staggered index operator. (See arXiv 1410.5733, 1203.2560)
// The index operator is i D_st - m \Gamma_5, where \Gamma_5 = i/2 (\Gamma_1 \Gamma_2 - \Gamma_2 \Gamma_1), where those are the symmshifts. 
void staggered_index_operator(complex<double>* lhs, complex<double>* rhs, void* extra_data)
{
    int i;
    complex<double> cplxI = complex<double>(0.0, 1.0);
    complex<double> cplxId2 = complex<double>(0.0, 0.5);
    
    staggered_u1_op* stagif = (staggered_u1_op*)extra_data;
    double mass = stagif->mass; 
    int x_fine = stagif->x_fine;
    int y_fine = stagif->y_fine;
    
    complex<double>* tmp = new complex<double>[stagif->x_fine*stagif->y_fine];
    complex<double>* tmp2 = new complex<double>[stagif->x_fine*stagif->y_fine];
    
    // Zero out lhs.
    for (i = 0; i < x_fine*y_fine; i++)
    {
        lhs[i] = 0.0;
    }
    
    // First, apply i D_st. We need a massless kernel for this.
    staggered_u1_op stag_nomass;
    stag_nomass.x_fine = x_fine;
    stag_nomass.y_fine = y_fine;
    stag_nomass.mass = 0.0;
    stag_nomass.lattice = stagif->lattice;
    stag_nomass.Nc = stagif->Nc;
    
    square_staggered_u1(lhs, rhs, (void*)&stag_nomass);
    for (i = 0; i < x_fine*y_fine; i++)
    {
        lhs[i] *= cplxI; 
    }
    
    
    // Next, form -m i/2 \Gamma_1 \Gamma_2. 
    
    staggered_symmshift_y(tmp, rhs, extra_data);
    staggered_symmshift_x(tmp2, tmp, extra_data);
    for (i = 0; i < x_fine*y_fine; i++)
    {
        lhs[i] -= mass*cplxId2*tmp2[i];
    }
    
    // Last, form + mi/2 \Gamma_2 \Gamma_1.
    
    staggered_symmshift_x(tmp, rhs, extra_data);
    staggered_symmshift_y(tmp2, tmp, extra_data);
    for (i = 0; i < x_fine*y_fine; i++)
    {
        lhs[i] += mass*cplxId2*tmp2[i];
    }
    
    delete[] tmp;
    delete[] tmp2; 
}


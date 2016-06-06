
#include <fstream>
#include <cmath>
#include <complex>
#include <random>

using namespace std; 

#include "u1_utils.h"

#define PI 3.14159265358979323846

// Load complex gauge field from file. 
// Based on Rich Brower's u1 gauge routines. 
// Reads in a U1 phase lattice from file, returns complex fields. 
// Rich's code has 'y' as the fast direction. Need to transpose!
void read_gauge_u1(complex<double>* gauge_field, int x_len, int y_len, string input_file)
{
   double phase_tmp;   
   fstream in_file;
   
   in_file.open(input_file,ios::in); 
   for(int x =0;x< x_len;x++)
   {
      for(int y =0;y< y_len;y++)
      {
         for(int mu=0; mu<2; mu++)
         {
            in_file >> phase_tmp;
            gauge_field[y*2*x_len+x*2+mu] = polar(1.0,phase_tmp);
            //cout << polar(1.0, phase_tmp) << "\n";
         }
      }
   }
   in_file.close(); 
   
   return;
}

// Write complex gauge field to file.
// Based on Rich Brower's u1 gauge routines.
// Writes a U1 phase lattice from file from complex fields.
// Rich's code has 'y' as the fast direction. Need to transpose!
void write_gauge_u1(complex<double>* gauge_field, int x_len, int y_len, string output_file)
{
    double phase_tmp;
    fstream out_file;
    
    out_file.open(output_file, ios::in|ios::out|ios::trunc);
    out_file.setf(ios_base::fixed, ios_base::floatfield);
    out_file.precision(20);
    
    for(int x =0;x< x_len;x++)
    {
        for(int y =0;y< y_len;y++)
        {
            for(int mu=0; mu<2; mu++)
            {
                phase_tmp = arg(gauge_field[y*2*x_len+x*2+mu]);
                out_file << phase_tmp << "\n";
            }
        }
    }
    out_file.close(); 
}

// Create a unit gauge field.
// Just set everything to 1!
void unit_gauge_u1(complex<double>* gauge_field, int x_len, int y_len)
{
    for (int i = 0; i < 2*y_len*x_len; i++)
    {
        gauge_field[i] = 1.0;
    }
}

// Create a hot gauge field, uniformly distributed in -Pi -> Pi.
// mt19937 can be created+seeded as: std::mt19937 generator (seed1);
void rand_gauge_u1(complex<double>* gauge_field, int x_len, int y_len, std::mt19937 &generator)
{
    // Generate a uniform distribution.
    std::uniform_real_distribution<> dist(-PI, PI);
    for (int i = 0; i < 2*y_len*x_len; i++)
    {
        gauge_field[i] = polar(1.0, dist(generator));
    }
}

// Create a gaussian gauge field with variance = 1/beta
// beta -> 0 is a hot start, beta -> inf is a cold start. 
// Based on code by Rich Brower, re-written for C++11.
void gauss_gauge_u1(complex<double>* gauge_field, int x_len, int y_len, std::mt19937 &generator, double beta)
{
    // Take abs value of beta.
    if (beta < 0) { beta = -beta; }
    
    // If beta is 0, just call a hot start.
    if (beta == 0)
    {
        rand_gauge_u1(gauge_field, x_len, y_len, generator);
    }
    
    // Generate a normal distribution with mean 0, stddev 1/sqrt(beta)
    std::normal_distribution<> dist(0.0, 1.0/sqrt(beta));
    
    // Generate fields!
    for (int i = 0; i < 2*y_len*x_len; i++)
    {
        gauge_field[i] = polar(1.0, dist(generator));
    }
    
}

// Generate a random gauge transform.
// mt19937 can be created+seeded as: std::mt19937 generator (seed1);
void rand_trans_u1(complex<double>* gauge_trans, int x_len, int y_len, std::mt19937 &generator)
{
    // Generate a uniform distribution.
    std::uniform_real_distribution<> dist(-PI, PI);
    for (int i = 0; i < y_len*x_len; i++)
    {
        gauge_trans[i] = polar(1.0, dist(generator));
    }
}

// Apply a gauge transform:
// u_i(x) = g(x) u_i(x) g^\dag(x+\hat{i})
void apply_gauge_trans_u1(complex<double>* gauge_field, complex<double>* gauge_trans, int x_len, int y_len)
{
    for (int y = 0; y < y_len; y++)
    {
        for (int x = 0; x < x_len; x++)
        {
            // Update x direction.
            gauge_field[2*x_len*y + 2*x] = gauge_trans[x_len*y + x]*gauge_field[2*x_len*y + 2*x]*conj(gauge_trans[x_len*y + (x+1)%x_len]);
            
            // Update y direction.
            gauge_field[2*x_len*y + 2*x + 1] = gauge_trans[x_len*y + x]*gauge_field[2*x_len*y + 2*x + 1]*conj(gauge_trans[x_len*((y+1)%y_len) + x]);
        }
    }
}

// Get average plaquette. 
complex<double> get_plaquette_u1(complex<double>* gauge_field, int x_len, int y_len)
{
    complex<double> plaq = 0.0;
    complex<double> tmp_plaq = 0.0;
    for (int y = 0; y < y_len; y++)
    {
        for (int x = 0; x < x_len; x++)
        {
            tmp_plaq = gauge_field[2*x_len*y + 2*x]*
                           gauge_field[2*x_len*y + 2*((x+1)%x_len) + 1]*
                      conj(gauge_field[2*x_len*((y+1)%y_len) + 2*x])*
                      conj(gauge_field[2*x_len*y + 2*x + 1]);
            plaq += tmp_plaq;
        }
    }
    return plaq / ((double)(x_len*y_len));
    
}



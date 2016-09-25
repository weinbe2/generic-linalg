// Test file for multi-rhs gaussian elimination and blockCG.

#include <iostream>
#include <iomanip> // to set output precision.
#include <cmath>
#include <string>
#include <sstream>
#include <complex>
#include <random>
#include <cstring> // should be replaced by using sstream

// Things for timing. 
#include <stdlib.h>
#include <unistd.h>
#include <time.h>


#include "generic_cg.h"
#include "generic_block_cg_mrhs.h"
#include "generic_gelim.h"
#include "generic_vector.h"
#include "verbosity.h"
#include "u1_utils.h"
#include "operators.h"

using namespace std; 

// Define pi.
#define PI 3.141592653589793

// Custom routine to load gauge field.
void internal_load_gauge_u1(complex<double>* lattice, int x_fine, int y_fine, double BETA);

// real square laplace
void square_laplace(double* lhs, double* rhs, void* extra_data);

// Timing routine.
timespec diff(timespec start, timespec end)
{
    timespec tmp;
    if ((end.tv_nsec-start.tv_nsec) < 0)
    {
        tmp.tv_sec = end.tv_sec-start.tv_sec-1;
        tmp.tv_nsec = 1000000000+end.tv_nsec-start.tv_nsec;
    }
    else
    {
        tmp.tv_sec = end.tv_sec-start.tv_sec;
        tmp.tv_nsec = end.tv_nsec-start.tv_nsec;
    }
    return tmp;
}

int main(int argc, char** argv)
{  
    // Declare some variables.
    cout << setiosflags(ios::scientific) << setprecision(6);
    int i,j,k;
    double **lhs, **rhs, *check; // For some Kinetic terms.
    double tmp; 
    double* tmp2; 
    std::mt19937 generator (1337u); // RNG, 1337u is the seed. 
    inversion_info invif;
    staggered_u1_op stagif;
    double bnorm; 
    stringstream sstream; 
    
    
    // First test: real multiRHS gaussian elimination.
    // We'll essentially test this by trying a matrix inversion and checking the result.
    // We'll also tag along another rhs. 
    // The matrix is {{1,2,3},{4,5,6},{7,8,10}}
    
    cout << "Real multiRHS gaussian elimination test.\n\n";
    
    {
        double** x = new double*[4];
        double** b = new double*[4];
        double** mat = new double*[3];
        for (i = 0; i < 4; i++)
        {
            x[i] = new double[3];
            b[i] = new double[3];
            if (i < 3)
            {
                mat[i] = new double[3];
            }
        }
        mat[0][0] = 1; mat[0][1] = 2; mat[0][2] = 3;
        mat[1][0] = 4; mat[1][1] = 5; mat[1][2] = 6;
        mat[2][0] = 7; mat[2][1] = 8; mat[2][2] = 10;

        // First rhs is e1, second rhs is e2, third rhs is e3.
        b[0][0] = 1; b[0][1] = 0; b[0][2] = 0; // e1.
        b[1][0] = 0; b[1][1] = 1; b[1][2] = 0; // e2.
        b[2][0] = 0; b[2][1] = 0; b[2][2] = 1; // e3.

        // Extra rhs for lulz.
        b[3][0] = 2; b[3][1] = 4; b[3][2] = 5; 

        // Do the multi-rhs solve.
        gaussian_elimination_multi_rhs(x, b, mat, 4, 3);

        // Check some things! First, the real inverse.
        cout << "Real inverse check.\n"; 
        for (i = 0; i < 3; i++)
        {
            for (j = 0; j < 3; j++)
            {
                double accum = 0;
                for (k = 0; k < 3; k++)
                {
                    accum += mat[i][k]*x[j][k]; 
                }
                cout << accum << " ";
            }
            cout << "\n";
        }

        // Next, the single rhs check.
        cout << "Real single rhs check.\n";
        double tot_resid = 0.0;
        for (i = 0; i < 3; i++)
        {
            double elem = 0.0;
            for (j = 0; j < 3; j++)
            {
                elem += mat[i][j]*x[3][j];
            }
            tot_resid += (elem - b[3][i])*(elem - b[3][i]);
        }
        cout << "The residual of Ax = {{2},{4},{5}} solve is " << sqrt(tot_resid) << ".\n";

        // Matrix inverse test.
        gaussian_elimination_matrix_inverse(mat, mat, 3);
        gaussian_elimination_matrix_inverse(x, mat, 3);
        cout << "Real explicit inverse check.\n";
        for (i = 0; i < 3; i++)
        {
            for (j = 0; j < 3; j++)
            {
                double accum = 0;
                for (k = 0; k < 3; k++)
                {
                    accum += mat[i][k]*x[k][j]; 
                }
                cout << accum << " ";
            }
            cout << "\n";
        }
        
        for (i = 0; i < 4; i++)
        {
            delete[] x[i];
            delete[] b[i];
            if (i < 3)
            {
                delete[] mat[i];
            }
        }
        delete[] x;
        delete[] b;
        delete[] mat; 
    }
    
    // Next test: complex multiRHS gaussian elimination.
    // We'll essentially test this by trying a matrix inversion and checking the result.
    // We'll also tag along another rhs. 
    // The matrix is {{1,2,3},{4,5,6},{7,8,10+i}}
    
    cout << "Real multiRHS gaussian elimination test.\n\n";
    
    {
        complex<double>** x = new complex<double>*[4];
        complex<double>** b = new complex<double>*[4];
        complex<double>** mat = new complex<double>*[3];
        for (i = 0; i < 4; i++)
        {
            x[i] = new complex<double>[3];
            b[i] = new complex<double>[3];
            if (i < 3)
            {
                mat[i] = new complex<double>[3];
            }
        }
        mat[0][0] = 1; mat[0][1] = 2; mat[0][2] = 3;
        mat[1][0] = 4; mat[1][1] = 5; mat[1][2] = 6;
        mat[2][0] = 7; mat[2][1] = 8; mat[2][2] = complex<double>(10,1);

        // First rhs is e1, second rhs is e2, third rhs is e3.
        b[0][0] = 1; b[0][1] = 0; b[0][2] = 0; // e1.
        b[1][0] = 0; b[1][1] = 1; b[1][2] = 0; // e2.
        b[2][0] = 0; b[2][1] = 0; b[2][2] = 1; // e3.

        // Extra rhs for lulz.
        b[3][0] = 2; b[3][1] = complex<double>(4.0,-3.0); b[3][2] = 5; 

        // Do the multi-rhs solve.
        gaussian_elimination_multi_rhs(x, b, mat, 4, 3);

        // Check some things! First, the real inverse.
        cout << "Complex inverse check.\n"; 
        for (i = 0; i < 3; i++)
        {
            for (j = 0; j < 3; j++)
            {
                complex<double> accum = 0;
                for (k = 0; k < 3; k++)
                {
                    accum += mat[i][k]*x[j][k]; 
                }
                cout << accum << " ";
            }
            cout << "\n";
        }

        // Next, the single rhs check.
        cout << "Complex single rhs check.\n";
        double tot_resid = 0.0;
        for (i = 0; i < 3; i++)
        {
            complex<double> elem = 0.0;
            for (j = 0; j < 3; j++)
            {
                elem += mat[i][j]*x[3][j];
            }
            tot_resid += real(conj(elem - b[3][i])*(elem - b[3][i]));
        }
        cout << "The residual of Ax = {{2},{4-3i},{5}} solve is " << sqrt(tot_resid) << ".\n";



        for (i = 0; i < 4; i++)
        {
            delete[] x[i];
            delete[] b[i];
            if (i < 3)
            {
                delete[] mat[i];
            }
        }
        delete[] x;
        delete[] b;
        delete[] mat; 
    }
    
    cout << "\nFinished gelim tests.\n\n"; 
    
    // Timing.
    timespec time1, time2, timediff;

    // Set parameters. 
    
    // L_x = L_y = Dimension for a square lattice.
    int square_size = 64; 
    
    // What masses should we use?
    double mass = 1e-2;
    
    // Solver precision.
    double outer_precision = 1e-10; 
    
    cout << "[OP]: Operator " << "Square Laplace" << " Mass " << mass << "\n";
    
    // Describe the fine lattice. 
    int x_fine = square_size;
    int y_fine = square_size;
    int fine_size = x_fine*y_fine;
    
    cout << "[VOL]: X " << x_fine << " Y " << y_fine << " Volume " << x_fine*y_fine;
    cout << "\n";
    
    // Do some allocation.
    // Initialize the lattice. Indexing: index = y*N + x.
    lhs = new double*[128]; // Try up to 128 lhs, rhs.
    rhs = new double*[128];
    for (i = 0; i < 128; i++)
    {
        lhs[i] = new double[fine_size]; zero<double>(lhs[i], fine_size);
        rhs[i] = new double[fine_size]; zero<double>(rhs[i], fine_size);
    }
    
    check = new double[fine_size];   
    tmp2 = new double[fine_size];
    // Zero it out.
    zero<double>(check, fine_size);
    zero<double>(tmp2, fine_size);
    //
    
    // Fill stagif.
    stagif.lattice = 0; // no lattice for now.
    stagif.mass = mass;  
    stagif.x_fine = x_fine;
    stagif.y_fine = y_fine; 
    stagif.Nc = 1; // Only relevant for laplace test only.
    
    // Create the verbosity structure.
    inversion_verbose_struct verb;
    verb.verbosity = VERB_DETAIL;
    verb.verb_prefix = "[BLOCKCG]: ";
    verb.precond_verbosity = VERB_DETAIL;
    verb.precond_verb_prefix = "Prec ";
    
    // Set a random source for all rhs. Should be the same every time b/c we hard code the seed above.
    for (i = 0; i < 128; i++)
    {
        gaussian<double>(rhs[i], fine_size, generator); 
    }
    
    // Orthonormalize vectors...
    //orthogonal<double>(rhs[1], rhs[0], fine_size);
    //normalize<double>(rhs[0], fine_size);
    //normalize<double>(rhs[1], fine_size);
    
    zero<double>(rhs[0], fine_size);
    rhs[0][3] = 1.0;
    zero<double>(rhs[1], fine_size);
    rhs[1][40] = 1.0;
    
    // Test an inversion.
    for (i = 0; i < 128; i++)
    {
        zero<double>(lhs[i], fine_size);
    }
    clock_gettime(CLOCK_PROCESS_CPUTIME_ID, &time1);
    invif = minv_vector_block_cg_mrhs(lhs, rhs, 2, fine_size, 100000, outer_precision, square_laplace, (void*)&stagif, &verb);
    clock_gettime(CLOCK_PROCESS_CPUTIME_ID, &time2); timediff = diff(time1, time2); cout << "Time " << ((double)(long int)(1000000000*timediff.tv_sec + timediff.tv_nsec))*1e-9 << "\n";
    
    // Check the inversion.
    verb.verb_prefix = "[CG]: ";
    for (i = 0; i < 128; i++)
    {
        zero<double>(lhs[i], fine_size);
    }
    clock_gettime(CLOCK_PROCESS_CPUTIME_ID, &time1);
    invif = minv_vector_cg(lhs[0], rhs[0], fine_size, 100000, outer_precision, square_laplace, (void*)&stagif, &verb);
    clock_gettime(CLOCK_PROCESS_CPUTIME_ID, &time2); timediff = diff(time1, time2); cout << "Time " << ((double)(long int)(1000000000*timediff.tv_sec + timediff.tv_nsec))*1e-9 << "\n";
    
    
    /*
    // Test BiCGStab-l inversions in powers of 2.
    for (i = 1; i <= 16; i*=2)
    {
        sstream.str(string()); // Clear the string stream.
        sstream << "[BICGSTAB-" << i << "]: "; // Set the string stream.
        verb.verb_prefix = sstream.str(); // Update the verbosity string.
        
        zero<double>(lhs, fine_size);
        clock_gettime(CLOCK_PROCESS_CPUTIME_ID, &time1);
        invif = minv_vector_bicgstab_l(lhs, rhs, fine_size, 100000, outer_precision, i, op, (void*)&stagif, &verb);
        clock_gettime(CLOCK_PROCESS_CPUTIME_ID, &time2); timediff = diff(time1, time2); cout << "Time " << ((double)(long int)(1000000000*timediff.tv_sec + timediff.tv_nsec))*1e-9 << "\n";
    }*/
    
    
    
    // Free the lattice.
    
    for (i = 0; i < 128; i++)
    {
        delete[] lhs[i];
        delete[] rhs[i];
    }
    
    delete[] lhs;
    delete[] rhs;
    delete[] check;
    delete[] tmp2; 
    
    
    return 0; 
}


// Custom routine to load gauge field.
void internal_load_gauge_u1(complex<double>* lattice, int x_fine, int y_fine, double BETA)
{
    if (x_fine == 32 && y_fine == 32)
    {
        if (abs(BETA - 3.0) < 1e-8)
        {
            read_gauge_u1(lattice, x_fine, y_fine, "../../multigrid/aa_mg/cfg/l32t32b30_heatbath.dat");
            cout << "[GAUGE]: Loaded a U(1) gauge field with non-compact beta = " << 1.0/sqrt(BETA) << "\n";
        }
        else if (abs(BETA - 6.0) < 1e-8)
        {
            read_gauge_u1(lattice, x_fine, y_fine, "../../multigrid/aa_mg/cfg/l32t32b60_heatbath.dat");
            cout << "[GAUGE]: Loaded a U(1) gauge field with non-compact beta = " << 1.0/sqrt(BETA) << "\n";
        }
        else if (abs(BETA - 10.0) < 1e-8)
        {
            read_gauge_u1(lattice, x_fine, y_fine, "../../multigrid/aa_mg/cfg/l32t32b100_heatbath.dat");
            cout << "[GAUGE]: Loaded a U(1) gauge field with non-compact beta = " << 1.0/sqrt(BETA) << "\n";
        }
        else if (abs(BETA - 10000.0) < 1e-8)
        {
            read_gauge_u1(lattice, x_fine, y_fine, "../../multigrid/aa_mg/cfg/l32t32bperturb_heatbath.dat");
            cout << "[GAUGE]: Loaded a U(1) gauge field with non-compact beta = " << 1.0/sqrt(BETA) << "\n";
        }
        else
        {
            cout << "[GAUGE]: Saved U(1) gauge field with correct beta, volume does not exist.\n";
        }
    }
    else if (x_fine == 64 && y_fine == 64)
    {
        if (abs(BETA - 3.0) < 1e-8)
        {
            read_gauge_u1(lattice, x_fine, y_fine, "../../multigrid/aa_mg/cfg/l64t64b30_heatbath.dat");
            cout << "[GAUGE]: Loaded a U(1) gauge field with non-compact beta = " << 1.0/sqrt(BETA) << "\n";
        }
        else if (abs(BETA - 6.0) < 1e-8)
        {
            read_gauge_u1(lattice, x_fine, y_fine, "../../multigrid/aa_mg/cfg/l64t64b60_heatbath.dat");
            cout << "[GAUGE]: Loaded a U(1) gauge field with non-compact beta = " << 1.0/sqrt(BETA) << "\n";
        }
        else if (abs(BETA - 10.0) < 1e-8)
        {
            read_gauge_u1(lattice, x_fine, y_fine, "../../multigrid/aa_mg/cfg/l64t64b100_heatbath.dat");
            cout << "[GAUGE]: Loaded a U(1) gauge field with non-compact beta = " << 1.0/sqrt(BETA) << "\n";
        }
        else if (abs(BETA - 10000.0) < 1e-8)
        {
            read_gauge_u1(lattice, x_fine, y_fine, "../../multigrid/aa_mg/cfg/l64t64bperturb_heatbath.dat");
            cout << "[GAUGE]: Loaded a U(1) gauge field with non-compact beta = " << 1.0/sqrt(BETA) << "\n";
        }
        else
        {
            cout << "[GAUGE]: Saved U(1) gauge field with correct beta, volume does not exist.\n";
        }
    }
    else if (x_fine == 128 && y_fine == 128)
    {
        if (abs(BETA - 6.0) < 1e-8)
        {
            read_gauge_u1(lattice, x_fine, y_fine, "../../multigrid/aa_mg/cfg/l128t128b60_heatbath.dat");
            cout << "[GAUGE]: Loaded a U(1) gauge field with non-compact beta = " << 1.0/sqrt(BETA) << "\n";
        }
        else if (abs(BETA - 10.0) < 1e-8)
        {
            read_gauge_u1(lattice, x_fine, y_fine, "../../multigrid/aa_mg/cfg/l128t128b100_heatbath.dat");
            cout << "[GAUGE]: Loaded a U(1) gauge field with non-compact beta = " << 1.0/sqrt(BETA) << "\n";
        }
        else if (abs(BETA - 10000.0) < 1e-8)
        {
            read_gauge_u1(lattice, x_fine, y_fine, "../../multigrid/aa_mg/cfg/l128t128bperturb_heatbath.dat");
            cout << "[GAUGE]: Loaded a U(1) gauge field with non-compact beta = " << 1.0/sqrt(BETA) << "\n";
        }
        else
        {
            cout << "[GAUGE]: Saved U(1) gauge field with correct beta, volume does not exist.\n";
        }
    }
    else
    {
        cout << "[GAUGE]: Saved U(1) gauge field with correct beta, volume does not exist.\n";
    }

}

// Square lattice.
// Kinetic term for a 2D laplace w/ period bc. Applies lhs = A*rhs.
// The unit vectors are e_1 = xhat, e_2 = yhat.
void square_laplace(double* lhs, double* rhs, void* extra_data)
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



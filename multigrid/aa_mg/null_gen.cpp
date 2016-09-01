
#include "null_gen.h"
#include "mg_complex.h"
#include "generic_vector.h"
#include "operators.h"

// Function to partition null vectors on the top level.
void null_partition_staggered(mg_operator_struct_complex* mgstruct, int num_null_vec, blocking_strategy bstrat, Lattice* Lat)
{
	int j; 
	int coord[2]; 
	int nd; 
	
	complex<double> *tmp, *tmp2; 
	
	switch (bstrat)
	{
		case BLOCK_EO:
			for (j = 0; j < mgstruct->curr_fine_size; j++)
			{
				if (Lat->index_is_even(j))
				{
					mgstruct->null_vectors[0][2*num_null_vec+1][j] = mgstruct->null_vectors[0][2*num_null_vec][j];
					mgstruct->null_vectors[0][2*num_null_vec][j] = 0.0;
				}
			}
			break;
		case BLOCK_TOPO:
			tmp = new complex<double>[mgstruct->curr_fine_size];
			tmp2 = new complex<double>[mgstruct->curr_fine_size];
			// Form vectors from (1 \pm \Gamma_5)/2.

			// Form \Gamma_5 null.

			// i/2 \Gamma_1 \Gamma_2
			zero<double>(tmp, mgstruct->curr_fine_size);
			zero<double>(tmp2, mgstruct->curr_fine_size);
			staggered_symmshift_y(tmp, mgstruct->null_vectors[0][2*num_null_vec], mgstruct->matrix_extra_data /* stagif */);
			staggered_symmshift_x(tmp2, tmp, mgstruct->matrix_extra_data);
			for (j = 0; j < mgstruct->curr_fine_size; j++)
			{
				mgstruct->null_vectors[0][2*num_null_vec+1][j] += complex<double>(0.0,0.5)*tmp2[j];
			}

			// -i/2 \Gamma_2 \Gamma_1
			zero<double>(tmp, mgstruct->curr_fine_size);
			zero<double>(tmp2, mgstruct->curr_fine_size);
			staggered_symmshift_x(tmp, mgstruct->null_vectors[0][2*num_null_vec], mgstruct->matrix_extra_data);
			staggered_symmshift_y(tmp2, tmp, mgstruct->matrix_extra_data);
			for (j = 0; j < mgstruct->curr_fine_size; j++)
			{
				mgstruct->null_vectors[0][2*num_null_vec+1][j] -= complex<double>(0.0,0.5)*tmp2[j];
			}

			// Form the two projectors.
			for (j = 0; j < mgstruct->curr_fine_size; j++)
			{
				mgstruct->null_vectors[0][2*num_null_vec][j] = 0.5*(mgstruct->null_vectors[0][2*num_null_vec][j]+mgstruct->null_vectors[0][2*num_null_vec+1][j]);
				mgstruct->null_vectors[0][2*num_null_vec+1][j] = mgstruct->null_vectors[0][2*num_null_vec][j]-mgstruct->null_vectors[0][2*num_null_vec+1][j];
			}
			delete[] tmp;
			delete[] tmp2; 

			break;
		case BLOCK_CORNER:

			for (j = 0; j < mgstruct->curr_fine_size; j++)
			{
				// Find x and y component. 
				Lat->index_to_coord(j, coord, nd);
				if (coord[0]%2 == 1 && coord[1]%2 == 0)
				{
					mgstruct->null_vectors[0][4*num_null_vec+1][j] = mgstruct->null_vectors[0][4*num_null_vec][j];
					mgstruct->null_vectors[0][4*num_null_vec][j] = 0.0;
				}
				else if (coord[0]%2 == 0 && coord[1]%2 == 1)
				{
					mgstruct->null_vectors[0][4*num_null_vec+2][j] = mgstruct->null_vectors[0][4*num_null_vec][j];
					mgstruct->null_vectors[0][4*num_null_vec][j] = 0.0;
				}
				else if (coord[0]%2 == 1 && coord[1]%2 == 1)
				{
					mgstruct->null_vectors[0][4*num_null_vec+3][j] = mgstruct->null_vectors[0][4*num_null_vec][j];
					mgstruct->null_vectors[0][4*num_null_vec][j] = 0.0;
				}
			}
			break;
		case BLOCK_NONE:
			// Nothing special to do.
			break;
	}
}

// Function to partition null vectors on the coarse level. Should update "Lattice" object for all levels.
void null_partition_coarse(complex<double>** null_vector, mg_operator_struct_complex* mgstruct, blocking_strategy bstrat)
{
	
}
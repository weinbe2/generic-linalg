
#include <random>
#include <iostream>

#include "mg.h"
#include "null_gen.h"
#include "mg_complex.h"
#include "generic_vector.h"
#include "operators.h"
#include "verbosity.h"

#include "generic_bicgstab.h"
#include "generic_cg.h"
#include "generic_cr.h"
#include "generic_gcr.h"
#include "generic_minres.h"
#include "generic_bicgstab_l.h"

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
				if (!Lat->index_is_even(j))
				{
					mgstruct->null_vectors[0][num_null_vec+mgstruct->n_vector/2][j] = mgstruct->null_vectors[0][num_null_vec][j];
					mgstruct->null_vectors[0][num_null_vec][j] = 0.0;
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
			staggered_symmshift_y(tmp, mgstruct->null_vectors[0][num_null_vec], mgstruct->matrix_extra_data /* stagif */);
			staggered_symmshift_x(tmp2, tmp, mgstruct->matrix_extra_data);
			for (j = 0; j < mgstruct->curr_fine_size; j++)
			{
				mgstruct->null_vectors[0][num_null_vec+mgstruct->n_vector/2][j] += complex<double>(0.0,0.5)*tmp2[j];
			}

			// -i/2 \Gamma_2 \Gamma_1
			zero<double>(tmp, mgstruct->curr_fine_size);
			zero<double>(tmp2, mgstruct->curr_fine_size);
			staggered_symmshift_x(tmp, mgstruct->null_vectors[0][num_null_vec], mgstruct->matrix_extra_data);
			staggered_symmshift_y(tmp2, tmp, mgstruct->matrix_extra_data);
			for (j = 0; j < mgstruct->curr_fine_size; j++)
			{
				mgstruct->null_vectors[0][num_null_vec+mgstruct->n_vector/2][j] -= complex<double>(0.0,0.5)*tmp2[j];
			}

			// Form the two projectors.
			for (j = 0; j < mgstruct->curr_fine_size; j++)
			{
				mgstruct->null_vectors[0][num_null_vec][j] = 0.5*(mgstruct->null_vectors[0][num_null_vec][j]+mgstruct->null_vectors[0][num_null_vec+mgstruct->n_vector/2][j]);
				mgstruct->null_vectors[0][num_null_vec+mgstruct->n_vector/2][j] = mgstruct->null_vectors[0][num_null_vec][j]-mgstruct->null_vectors[0][num_null_vec+mgstruct->n_vector/2][j];
			}
			delete[] tmp;
			delete[] tmp2; 

			break;
		case BLOCK_CORNER:

			for (j = 0; j < mgstruct->curr_fine_size; j++)
			{
				// Find x and y component. 
				Lat->index_to_coord(j, coord, nd);
				if (coord[0]%2 == 1 && coord[1]%2 == 0) // odd corner
				{
					mgstruct->null_vectors[0][num_null_vec+2*mgstruct->n_vector/4][j] = mgstruct->null_vectors[0][num_null_vec][j];
					mgstruct->null_vectors[0][num_null_vec][j] = 0.0;
				}
				else if (coord[0]%2 == 0 && coord[1]%2 == 1) // odd corner
				{
					mgstruct->null_vectors[0][num_null_vec+3*mgstruct->n_vector/4][j] = mgstruct->null_vectors[0][num_null_vec][j];
					mgstruct->null_vectors[0][num_null_vec][j] = 0.0;
				}
				else if (coord[0]%2 == 1 && coord[1]%2 == 1) // even corner
				{
					mgstruct->null_vectors[0][num_null_vec+mgstruct->n_vector/4][j] = mgstruct->null_vectors[0][num_null_vec][j];
					mgstruct->null_vectors[0][num_null_vec][j] = 0.0;
				}
			}
			break;
		case BLOCK_NONE:
			// Nothing special to do.
			break;
	}
}

// Function to partition null vectors on the coarse level. Should update "Lattice" object for all levels.
void null_partition_coarse(mg_operator_struct_complex* mgstruct, int num_null_vec, blocking_strategy bstrat)
{
	int j, c;
	
	switch (bstrat)
	{
		case BLOCK_EO:
		case BLOCK_TOPO:
			for (j = 0; j < mgstruct->curr_fine_size; j++)
			{
				c = j % mgstruct->n_vector; // What color index do we have?
											   // 0 to mgstruct->n_vector/2-1 is even, else is odd.
				//int x_coord = (i - c)/mgstruct->n_vector % mgstruct->curr_x_fine;
				//int y_coord = ((i - c)/mgstruct->n_vector - x_coord)/mgstruct->curr_x_fine;

				// If c is >= mgstruct->n_vector/2, it's odd!

				if (c >= mgstruct->n_vector/2)
				{
					mgstruct->null_vectors[mgstruct->curr_level][num_null_vec+mgstruct->n_vector/2][j] = mgstruct->null_vectors[mgstruct->curr_level][num_null_vec][j];
					mgstruct->null_vectors[mgstruct->curr_level][num_null_vec][j] = 0.0;
				}
			}
			break;
		case BLOCK_CORNER:
			for (j = 0; j < mgstruct->curr_fine_size; j++)
			{
				c = j % mgstruct->n_vector; // What color index do we have?
											   // 0 to mgstruct->n_vector/2-1 is even, else is odd.
				//int x_coord = (i - c)/mgstruct->n_vector % mgstruct->curr_x_fine;
				//int y_coord = ((i - c)/mgstruct->n_vector - x_coord)/mgstruct->curr_x_fine;


				if (c >= mgstruct->n_vector/4 && c < 2*mgstruct->n_vector/4)
				{
					mgstruct->null_vectors[mgstruct->curr_level][num_null_vec+mgstruct->n_vector/4][j] = mgstruct->null_vectors[mgstruct->curr_level][num_null_vec][j];
					mgstruct->null_vectors[mgstruct->curr_level][num_null_vec][j] = 0.0;
				}
				else if (c >= 2*mgstruct->n_vector/4 && c < 3*mgstruct->n_vector/4)
				{
					mgstruct->null_vectors[mgstruct->curr_level][num_null_vec+2*mgstruct->n_vector/4][j] = mgstruct->null_vectors[mgstruct->curr_level][num_null_vec][j];
					mgstruct->null_vectors[mgstruct->curr_level][num_null_vec][j] = 0.0;
				}
				else if (c >= 3*mgstruct->n_vector/4)
				{
					mgstruct->null_vectors[mgstruct->curr_level][num_null_vec+3*mgstruct->n_vector/4][j] = mgstruct->null_vectors[mgstruct->curr_level][num_null_vec][j];
					mgstruct->null_vectors[mgstruct->curr_level][num_null_vec][j] = 0.0;
				}
			}

			break;
		case BLOCK_NONE:
			// Nothing to do here...
			break;
	}
}

// Function to generate free field null vectors.
void null_generate_random_smooth(mg_operator_struct_complex* mgstruct, null_vector_params* nvec_params, bool do_gauge_transform, std::complex<double>* gauge_trans)
{
	int i,k;
	
	// Construct a single vector, partition as needed.
	for (i = 0; i < mgstruct->curr_fine_size; i++)
	{
		mgstruct->null_vectors[mgstruct->curr_level][0][i] = 1;
		if (do_gauge_transform && mgstruct->curr_level == 0) 
		{
			mgstruct->null_vectors[0][0][i] *= gauge_trans[i];
		}
	}

	// Aggregate in chirality (or corners) as needed.  
	// This is handled differently if we're on the top level or further down. 
	if (mgstruct->curr_level == 0) // top level routines
	{
		null_partition_staggered(mgstruct, 0, nvec_params->bstrat, mgstruct->latt[0]);
	}
	else // not on the top level
	{
		null_partition_coarse(mgstruct, 0, nvec_params->bstrat);
	}

	for (k = 0; k < nvec_params->null_partitions; k++)
	{
		normalize(mgstruct->null_vectors[mgstruct->curr_level][k*nvec_params->n_null_vector], mgstruct->curr_fine_size);
	}
}

// Function to generate null vectors by throwing a random source and smoothing.
void null_generate_random_smooth(mg_operator_struct_complex* mgstruct, null_vector_params* nvec_params, inversion_verbose_struct* verb, std::mt19937* generator)
{
	int i,j,k;
	inversion_info invif; 
	
	// Populate an inverter struct.
	solver_params solve;
	solve.tol = nvec_params->null_precisions[mgstruct->curr_level];
	solve.max_iters = nvec_params->null_max_iters[mgstruct->curr_level];
	solve.restart = nvec_params->null_restart;
	solve.restart_freq = nvec_params->null_restart_freq;
	solve.minres_omega = 0.67; // should expose.
	solve.bicgstabl_l = nvec_params->null_bicgstab_l;
	
	// We generate null vectors by solving Ax = 0, with a
	// gaussian initial guess.
	// For sanity with the residual, we really solve Ax = -Ax_0,
	// where x has a zero initial guess, x_0 is a random vector.
	complex<double>* rand_guess = new complex<double>[mgstruct->curr_fine_size];
	complex<double>* Arand_guess = new complex<double>[mgstruct->curr_fine_size];

	for (i = 0; i < mgstruct->n_vector/nvec_params->null_partitions; i++)
	{
		// Create a gaussian random source. 
		gaussian<double>(rand_guess, mgstruct->curr_fine_size, *generator);

		// Make orthogonal to previous solutions. May not be that necessary.
		if (i > 0) // If there are vectors to orthogonalize against...
		{
			for (j = 0; j < i; j++) // Iterate over all of them...
			{
				for (k = 0; k < (nvec_params->do_ortho_eo ? nvec_params->null_partitions : 1); k++) // And then iterate over even/odd or corners!
				{
					orthogonal<double>(rand_guess, mgstruct->null_vectors[mgstruct->curr_level][j+k*nvec_params->n_null_vector], mgstruct->curr_fine_size);
					if (nvec_params->do_global_ortho_conj)
					{
						conj<double>(mgstruct->null_vectors[mgstruct->curr_level][j+k*nvec_params->n_null_vector], mgstruct->curr_fine_size);
						orthogonal<double>(rand_guess, mgstruct->null_vectors[mgstruct->curr_level][j+k*nvec_params->n_null_vector], mgstruct->curr_fine_size);
						conj<double>(mgstruct->null_vectors[mgstruct->curr_level][j+k*nvec_params->n_null_vector], mgstruct->curr_fine_size);
					}
				}
			}
		}

		// Solve the residual equation. 
		zero<double>(Arand_guess, mgstruct->curr_fine_size);

		fine_square_staggered(Arand_guess, rand_guess, (void*)mgstruct); mgstruct->dslash_count->nullvectors[mgstruct->curr_level]++;

		for (j = 0; j < mgstruct->curr_fine_size; j++)
		{
		   Arand_guess[j] = -Arand_guess[j]; 
		}
		zero<double>(mgstruct->null_vectors[mgstruct->curr_level][i], mgstruct->curr_fine_size);

		// Encapsulating function for all null space solvers.
		invif = minv_unpreconditioned(mgstruct->null_vectors[mgstruct->curr_level][i], Arand_guess, mgstruct->curr_fine_size, nvec_params->null_gen, solve, fine_square_staggered, (void*)mgstruct, verb);

		mgstruct->dslash_count->nullvectors[mgstruct->curr_level] += invif.ops_count; 

		for (j = 0; j < mgstruct->curr_fine_size; j++)
		{
			mgstruct->null_vectors[mgstruct->curr_level][i][j] += rand_guess[j];
		}

		// Split into eo now if need be, otherwise we do it later.
		if (nvec_params->do_ortho_eo)
		{
			// Aggregate in chirality (or corners) as needed.  
			// This is handled differently if we're on the top level or further down. 
			
			if (mgstruct->curr_level == 0) // top level routines
			{
				null_partition_staggered(mgstruct, i, nvec_params->bstrat, mgstruct->latt[0]);
			}
			else // not on the top level
			{
				null_partition_coarse(mgstruct, i, nvec_params->bstrat);
			}
		}


		// Normalize new vectors.
		for (k = 0; k < (nvec_params->do_ortho_eo ? nvec_params->null_partitions : 1); k++)
		{
			normalize(mgstruct->null_vectors[mgstruct->curr_level][i+k*nvec_params->n_null_vector], mgstruct->curr_fine_size);
		}

		// Orthogonalize against previous vectors. 
		if (i > 0)
		{
			for (j = 0; j < i; j++)
			{
				cout << "[L" << mgstruct->curr_level+1 << "_NULLVEC]: Pre-orthog cosines of " << j << "," << i << " are: "; 
				for (k = 0; k < (nvec_params->do_ortho_eo ? nvec_params->null_partitions : 1); k++)
				{
					// Check dot product before normalization.
					cout << abs(dot<double>(mgstruct->null_vectors[mgstruct->curr_level][i+k*nvec_params->n_null_vector], mgstruct->null_vectors[mgstruct->curr_level][j+k*nvec_params->n_null_vector], mgstruct->curr_fine_size)/sqrt(norm2sq<double>(mgstruct->null_vectors[mgstruct->curr_level][i+k*nvec_params->n_null_vector],mgstruct->curr_fine_size)*norm2sq<double>(mgstruct->null_vectors[mgstruct->curr_level][j+k*nvec_params->n_null_vector],mgstruct->curr_fine_size))) << " ";

					orthogonal<double>(mgstruct->null_vectors[mgstruct->curr_level][i+k*nvec_params->n_null_vector], mgstruct->null_vectors[0][j+k*nvec_params->n_null_vector], mgstruct->curr_fine_size); 
					if (nvec_params->do_global_ortho_conj)
					{
						conj<double>(mgstruct->null_vectors[mgstruct->curr_level][j+k*nvec_params->n_null_vector], mgstruct->curr_fine_size);
						orthogonal<double>(mgstruct->null_vectors[mgstruct->curr_level][i+k*nvec_params->n_null_vector], mgstruct->null_vectors[0][j+k*nvec_params->n_null_vector], mgstruct->curr_fine_size); 
						conj<double>(mgstruct->null_vectors[mgstruct->curr_level][j+k*nvec_params->n_null_vector], mgstruct->curr_fine_size);
					}
				}
				cout << "\n";
			}
		}

		// Normalize again.
		for (k = 0; k < (nvec_params->do_ortho_eo ? nvec_params->null_partitions : 1); k++)
		{
			normalize(mgstruct->null_vectors[mgstruct->curr_level][i+k*nvec_params->n_null_vector], mgstruct->curr_fine_size);
		}

	}

	// If we didn't split null vectors before, we do it now.
	if (!nvec_params->do_ortho_eo)
	{
		for (i = 0; i < mgstruct->n_vector/nvec_params->null_partitions; i++)
		{
			// Aggregate in chirality (or corners) as needed.  
			// This is handled differently if we're on the top level or further down. 
			if (mgstruct->curr_level == 0) // top level routines
			{
				null_partition_staggered(mgstruct, i, nvec_params->bstrat, mgstruct->latt[0]);
			}
			else // not on the top level
			{
				null_partition_coarse(mgstruct, i, nvec_params->bstrat);
			}

			for (k = 0; k < nvec_params->null_partitions; k++)
			{
				normalize(mgstruct->null_vectors[mgstruct->curr_level][i+k*nvec_params->n_null_vector], mgstruct->curr_fine_size);
			}
		}
	}

	delete[] rand_guess; 
	delete[] Arand_guess; 
}

inversion_info minv_unpreconditioned(complex<double>* lhs, complex<double>* rhs, int size, mg_null_gen_type type, solver_params& params, void (*matrix_vector)(complex<double>*,complex<double>*,void*), void* extra_info, inversion_verbose_struct* verb)
{
	switch (type)
	{
		case NULL_GCR:
			if (params.restart)
			{
				return minv_vector_gcr_restart(lhs, rhs, size, params.max_iters, params.tol, params.restart_freq, matrix_vector, extra_info, verb);
			}
			else
			{
				return minv_vector_gcr(lhs, rhs, size, params.max_iters, params.tol, matrix_vector, extra_info, verb);
			}
			break;
		case NULL_BICGSTAB:
			if (params.restart)
			{
				return minv_vector_bicgstab_restart(lhs, rhs, size, params.max_iters, params.tol, params.restart_freq, matrix_vector, extra_info, verb);
			}
			else
			{
				return minv_vector_bicgstab(lhs, rhs, size, params.max_iters, params.tol, matrix_vector, extra_info, verb);
			}
			break;
		case NULL_CG:
			if (params.restart) // why would you do this I don't know it's CG come on
			{
				return minv_vector_cg_restart(lhs, rhs, size, params.max_iters, params.tol, params.restart_freq, matrix_vector, extra_info, verb);
			}
			else
			{
				return minv_vector_cg(lhs, rhs, size, params.max_iters, params.tol, matrix_vector, extra_info, verb);
			}
			break;
		case NULL_MINRES:
			// Restarting doesn't make sense for MinRes. 
			return minv_vector_minres(lhs, rhs, size, params.max_iters, params.tol, params.minres_omega, matrix_vector, extra_info, verb);
			break;
		case NULL_ARPACK: // it can't get here. 
			break;
		case NULL_BICGSTAB_L:
			if (params.restart)
			{
				return minv_vector_bicgstab_l_restart(lhs, rhs, size, params.max_iters, params.tol, params.restart_freq, params.bicgstabl_l, matrix_vector, extra_info, verb);
			}
			else
			{
				return minv_vector_bicgstab_l(lhs, rhs, size, params.max_iters, params.tol, params.bicgstabl_l, matrix_vector, extra_info, verb);
			}
			break;
		default:
			return inversion_info();
	}

	return inversion_info();
}



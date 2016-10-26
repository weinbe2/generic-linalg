
#include "string.h"

#include "input_params.h"


// Print usage. 
void display_usage()
{
    cout << "--help\n";
    cout << "--lattice-size [32, 64, 128] {##}                 (default 32x32)\n";
    cout << "--operator [laplace, laplace2, staggered\n";
    cout << "       g5_staggered, normal_staggered, index]     (default staggered)\n";
    cout << "--outer-solver [gcr, bicgstab, cg]                (default gcr)\n";
    cout << "--outer-precision [outer prec]                    (default 5e-7)\n";
    cout << "--outer-max-iter [outer maximum iterations]       (default 100000)\n";
    cout << "--outer-restart [yes, no]                         (default yes)\n";
    cout << "--outer-restart-freq [#]                          (default 64)\n";
    cout << "--null-operator [laplace, laplace2, staggered\n";
    cout << "       g5_staggered, normal_staggered, index]     (default staggered)\n";
    cout << "--null-solver [gcr, bicgstab, cg, cr, minres,\n";
    cout << "       arpack, bicgstab-l]                        (default bicgstab)\n";
    cout << "--null-precision [null prec] {#, #...}            (default 5e-5)\n";
    cout << "--null-max-iter [null nax iter] {#, #...}         (default 500)\n";
    cout << "--null-restart [yes, no]                          (default no)\n";
    cout << "--null-restart-freq [#]                           (default 8 if restarting is enabled)\n";
    cout << "--null-bicgstab-l [#]                             (default 4 if using bicgstab-l)\n"; 
    //cout << "--null-mass [mass] {#, #...}                      (default 1e-4)\n";
    cout << "--null-mass [mass]                                (default 1e-4)\n";
    cout << "--null-eo [corner, yes, no, topo]                 (default yes)\n";
    cout << "--null-global-ortho-conj [yes, no]                (default no, it only helps in some weird fluke cases)\n";
    cout << "--null-ortho-eo [yes, no]                         (default no)\n"; 
    cout << "--mass [mass]                                     (default 1e-2)\n";
    cout << "--blocksize [blocksize] {#, #...}                 (default 4, same for all levels)\n";
    cout << "--nvec [nvec]                                     (default 4)\n";
    cout << "--nrefine [number coarse]                         (default 1)\n";
    cout << "--cycle-type [v, k]                               (default v)\n";
    cout << "--mg-type [self, normal_eqn]                      (default self)\n";
    cout << "--smoother-solver [gcr, bicgstab, cg, minres,\n";
    cout << "         cr, bicgstab-l, none]                    (default gcr)\n"; 
    cout << "--smoother-type [self, normal_eqn]                (default self)\n";
    cout << "--npre-smooth [presmooth steps] {#, #, ...}       (default 6, same for all levels)\n";
    cout << "--npost-smooth [postsmooth steps] {#, #, ...}     (default 6, same for all levels)\n";
    cout << "--coarse-solver [gcr, bicgstab, cg, minres,\n";
    cout << "         cr, none]                                (default gcr)\n";
    cout << "--coarse-precision [coarse precision] {#, #, ...} (default 1e-2, same for all levels)\n";
    cout << "--coarse-max-iter [coarse maximum iterations]     (default 1024)\n";
    cout << "--gauge [unit, load, random]                      (default load)\n";
    cout << "--gauge-transform [yes, no]                       (default no)\n";
    cout << "--beta [3.0, 6.0, 10.0, 10000.0]                  (default 6.0)\n";
    cout << "--load-cfg [path]                                 (default do not load, overrides beta)\n";
    cout << "--do-freetest [yes, no]                           (default no)\n"; 
    cout << "--do-eigentest [yes, no]                          (default no)\n";
    cout << "--n-ev [all, # smallest]                          (default all)\n";
    cout << "--n-cv [#]                                        (default min(all, 2.5*n-ev))\n";
}

// Parse inputs, put into mg_input_params.
int parse_inputs(int argc, char** argv, mg_input_params *params)
{
    for (int i = 1; i < argc; i++)
    {
        if (strcmp(argv[i], "--help") == 0)
        {
            display_usage();
            return 0;
        }
        if (i+1 != argc)
        {
            if (strcmp(argv[i], "--mass") == 0)
            {
                params->mass = atof(argv[i+1]);
                i++;
            }
            else if (strcmp(argv[i], "--operator") == 0)
            {
                if (strcmp(argv[i+1], "laplace") == 0)
                {
                    params->opt = LAPLACE; 
                }
                else if (strcmp(argv[i+1], "laplace2") == 0)
                {
                    params->opt = LAPLACE_NC2; 
                }
                else if (strcmp(argv[i+1], "staggered") == 0)
                {
                    params->opt = STAGGERED; 
                }
                else if (strcmp(argv[i+1], "g5_staggered") == 0)
                {
                    params->opt = G5_STAGGERED; 
                }
                else if (strcmp(argv[i+1], "normal_staggered") == 0)
                {
                    params->opt = STAGGERED_NORMAL;
                }
                else if (strcmp(argv[i+1], "index") == 0)
                {
                    params->opt = STAGGERED_INDEX;
                }
                i++;
            }
            else if (strcmp(argv[i], "--outer-solver") == 0)
            {
                if (strcmp(argv[i+1], "gcr") == 0)
                {
                    params->out_solve = OUTER_GCR;
                }
                else if (strcmp(argv[i+1], "bicgstab") == 0)
                {
                    params->out_solve = OUTER_BICGSTAB;
                }
                else if (strcmp(argv[i+1], "cg") == 0)
                {
                    params->out_solve = OUTER_CG;
                }
                i++;
            }
            else if (strcmp(argv[i], "--outer-precision") == 0)
            {
                params->outer_precision = atof(argv[i+1]);
                i++;
            }
            else if (strcmp(argv[i], "--outer-max-iter") == 0)
            {
                params->outer_max_iter = atoi(argv[i+1]);
                i++;
            }
            else if (strcmp(argv[i], "--outer-restart") == 0)
            {
                if (strcmp(argv[i+1], "yes") == 0)
                {
                    params->outer_restart = true;
                }
                else if (strcmp(argv[i+1], "no") == 0)
                {
                    params->outer_restart = false;
                }
                i++;
            }
            else if (strcmp(argv[i], "--outer-restart-freq") == 0)
            {
                params->outer_restart_freq = atoi(argv[i+1]);
                i++;
            }
            else if (strcmp(argv[i], "--null-operator") == 0)
            {
                if (strcmp(argv[i+1], "laplace") == 0)
                {
                    params->nvec_params.opt_null = LAPLACE; 
                }
                else if (strcmp(argv[i+1], "laplace2") == 0)
                {
                    params->nvec_params.opt_null = LAPLACE_NC2; 
                }
                else if (strcmp(argv[i+1], "staggered") == 0)
                {
                    params->nvec_params.opt_null = STAGGERED; 
                }
                else if (strcmp(argv[i+1], "g5_staggered") == 0)
                {
                    params->nvec_params.opt_null = G5_STAGGERED; 
                }
                else if (strcmp(argv[i+1], "normal_staggered") == 0)
                {
                    params->nvec_params.opt_null = STAGGERED_NORMAL;
                }
                else if (strcmp(argv[i+1], "index") == 0)
                {
                    params->nvec_params.opt_null = STAGGERED_INDEX;
                }
                i++;
            }
            else if (strcmp(argv[i], "--blocksize") == 0)
            {
                params->blocksizes[0] = atoi(argv[i+1]);
                i++;
                while (i+1 != argc && argv[i+1][0] != '-')
                {
                    params->blocksizes.push_back(atoi(argv[i+1]));
                    i++;
                }
            }
            else if (strcmp(argv[i], "--nvec") == 0)
            {
                params->nvec_params.n_null_vector = atoi(argv[i+1]);
                i++;
            }
            else if (strcmp(argv[i], "--null-precision") == 0)
            {
                params->nvec_params.null_precisions[0] = atof(argv[i+1]);
                i++;
                while (i+1 != argc && argv[i+1][0] != '-')
                {
                    params->nvec_params.null_precisions.push_back(atof(argv[i+1]));
                    i++;
                }
            }
            else if (strcmp(argv[i], "--null-max-iter") == 0)
            {
                params->nvec_params.null_max_iters[0] = atoi(argv[i+1]);
                i++;
                while (i+1 != argc && argv[i+1][0] != '-')
                {
                    params->nvec_params.null_max_iters.push_back(atoi(argv[i+1]));
                    i++;
                }
            }
            else if (strcmp(argv[i], "--null-restart") == 0)
            {
                if (strcmp(argv[i+1], "yes") == 0)
                {
                    params->nvec_params.null_restart = true;
                }
                else if (strcmp(argv[i+1], "no") == 0)
                {
                    params->nvec_params.null_restart = false;
                }
                i++;
            }
            else if (strcmp(argv[i], "--null-restart-freq") == 0)
            {
                params->nvec_params.null_restart_freq = atoi(argv[i+1]);
                i++;
            }
            else if (strcmp(argv[i], "--null-bicgstab-l") == 0)
            {
                params->nvec_params.null_bicgstab_l = atoi(argv[i+1]);
                i++;
            }
            else if (strcmp(argv[i], "--null-mass") == 0)
            {
                params->nvec_params.null_mass = atof(argv[i+1]);
                i++;
            }
            else if (strcmp(argv[i], "--null-solver") == 0)
            {
                if (strcmp(argv[i+1], "gcr") == 0)
                {
                    params->nvec_params.null_gen = MINV_GCR;
                }
                else if (strcmp(argv[i+1], "bicgstab") == 0)
                {
                    params->nvec_params.null_gen = MINV_BICGSTAB;
                }
                else if (strcmp(argv[i+1], "cg") == 0)
                {
                    params->nvec_params.null_gen = MINV_CG;
                }
                else if (strcmp(argv[i+1], "cr") == 0)
                {
                    params->nvec_params.null_gen = MINV_CR;
                }
                else if (strcmp(argv[i+1], "minres") == 0)
                {
                    params->nvec_params.null_gen = MINV_MINRES;
                }
                else if (strcmp(argv[i+1], "arpack") == 0)
                {
                    params->nvec_params.null_gen = MINV_INVALID;
                    params->nvec_use_eigen = true; 
                    //cout << "[ERROR]: Cannot use eigenvectors as null vectors without arpack bindings.\n";
                    return 0;
                }
                else if (strcmp(argv[i+1], "bicgstab-l") == 0)
                {
                    params->nvec_params.null_gen = MINV_BICGSTAB_L;
                }
                i++;
            }
            else if (strcmp(argv[i], "--null-eo") == 0)
            {
                if (strcmp(argv[i+1], "yes") == 0)
                {
                    params->nvec_params.bstrat = BLOCK_EO; // even/odd
                }
                else if (strcmp(argv[i+1], "corner") == 0)
                {
                    params->nvec_params.bstrat = BLOCK_CORNER; // corners
                }
                else if (strcmp(argv[i+1], "topo") == 0)
                {
                    params->nvec_params.bstrat = BLOCK_TOPO; // chirality as defined by taste singlet. 
                }
                else // none.
                {
                    params->nvec_params.bstrat = BLOCK_NONE; // fully reduce. 
                }
                i++;
            }
            else if (strcmp(argv[i], "--null-global-ortho-conj") == 0)
            {
                if (strcmp(argv[i+1], "yes") == 0)
                {
                    params->nvec_params.do_global_ortho_conj = true; // yes, globally orthogonalize null vectors against previous and conj.
                }
                else if (strcmp(argv[i+1], "no") == 0)
                {
                    params->nvec_params.do_global_ortho_conj = false;
                }
                i++;
            }
            else if (strcmp(argv[i], "--null-ortho-eo") == 0)
            {
                if (strcmp(argv[i+1], "yes") == 0)
                {
                    params->nvec_params.do_ortho_eo = true; // yes, split null vectors into eo before orthogonalizing.
                }
                else if (strcmp(argv[i+1], "no") == 0)
                {
                    params->nvec_params.do_ortho_eo = false; 
                }
                i++;
            }
            else if (strcmp(argv[i], "--nrefine") == 0)
            {
                params->n_refine = atoi(argv[i+1]);
                i++;
            }
            else if (strcmp(argv[i], "--multi-strategy") == 0)
            {
                if (strcmp(argv[i+1], "smooth") == 0)
                {
                    params->mlevel_type = MLEVEL_SMOOTH;
                }
                else if (strcmp(argv[i+1], "recursive") == 0)
                {
                    params->mlevel_type = MLEVEL_RECURSIVE; 
                }
                i++;
            }
            else if (strcmp(argv[i], "--cycle-type") == 0) // aliased to multi-strategy
            {
                if (strcmp(argv[i+1], "v") == 0)
                {
                    params->mlevel_type = MLEVEL_SMOOTH;
                }
                else if (strcmp(argv[i+1], "k") == 0)
                {
                    params->mlevel_type = MLEVEL_RECURSIVE; 
                }
                i++;
            }
            else if (strcmp(argv[i], "--coarse-solver") == 0)
            {
                if (strcmp(argv[i+1], "cr") == 0)
                {
                    params->in_solve = CR;
                }
                else if (strcmp(argv[i+1], "bicgstab") == 0)
                {
                    params->in_solve = BICGSTAB;
                }
                else if (strcmp(argv[i+1], "cg") == 0)
                {
                    params->in_solve = CG;
                }
                else if (strcmp(argv[i+1], "gcr") == 0)
                {
                    params->in_solve = GCR;
                }
                else if (strcmp(argv[i+1], "minres") == 0)
                {
                    params->in_solve = MINRES;
                }
                else if (strcmp(argv[i+1], "none") == 0)
                {
                    params->in_solve = NONE;
                }
                i++;
            }
            else if (strcmp(argv[i], "--coarse-precision") == 0)
            {
                params->inner_precisions[0] = atof(argv[i+1]);
                i++;
                while (i+1 != argc && argv[i+1][0] != '-')
                {
                    params->inner_precisions.push_back(atof(argv[i+1]));
                    i++;
                }
            }
            else if (strcmp(argv[i], "--coarse-max-iter") == 0)
            {
                params->inner_max = atoi(argv[i+1]);
                i++;
            }
            else if (strcmp(argv[i], "--gauge") == 0)
            {
                if (strcmp(argv[i+1], "unit") == 0)
                {
                    params->gauge_load = GAUGE_UNIT;
                }
                else if (strcmp(argv[i+1], "random") == 0)
                {
                    params->gauge_load = GAUGE_RANDOM;
                }
                else
                {
                    params->gauge_load = GAUGE_LOAD;
                }
                i++;
            }
            else if (strcmp(argv[i], "--gauge-transform") == 0)
            {
                if (strcmp(argv[i+1], "yes") == 0)
                {
                    params->do_gauge_transform = true;
                }
                else
                {
                    params->do_gauge_transform = false;
                }
                i++;
            }
            else if (strcmp(argv[i], "--beta") == 0)
            {
                params->beta = atof(argv[i+1]);
                i++;
            }
            else if (strcmp(argv[i], "--lattice-size") == 0)
            {
                params->lattice_size_x = atoi(argv[i+1]);
                
                if (i+2 != argc)
                {
                    if (argv[i+2][0] == '-' && argv[i+2][1] == '-') // look for --
                    {
                        params->lattice_size_y = params->lattice_size_x;
                    }
                    else // assume number
                    {
                        params->lattice_size_y = atoi(argv[i+2]);
                        i++;
                    }
                }
                else
                {
                    params->lattice_size_y = params->lattice_size_x; // At the end, don't try to grab the next element!
                }
                i++;
            }
            else if (strcmp(argv[i], "--mg-type") == 0)
            {
                if (strcmp(argv[i+1], "self") == 0)
                {
                    params->normal_eqn_mg = false;
                }
                else if (strcmp(argv[i+1], "normal_eqn") == 0)
                {
                    params->normal_eqn_mg = true;
                }
                i++;
            }
            else if (strcmp(argv[i], "--smoother-solver") == 0)
            {
                if (strcmp(argv[i+1], "gcr") == 0)
                {
                    params->in_smooth = MINV_GCR;
                }
                else if (strcmp(argv[i+1], "bicgstab") == 0)
                {
                    params->in_smooth = MINV_BICGSTAB;
                }
                else if (strcmp(argv[i+1], "cg") == 0)
                {
                    params->in_smooth = MINV_CG;
                }
                else if (strcmp(argv[i+1], "cr") == 0)
                {
                    params->in_smooth = MINV_CR;
                }
                else if (strcmp(argv[i+1], "minres") == 0)
                {
                    params->in_smooth = MINV_MINRES;
                }
                else if (strcmp(argv[i+1], "bicgstab-l") == 0)
                {
                    params->in_smooth = MINV_BICGSTAB_L;
                }
                else if (strcmp(argv[i+1], "none") == 0)
                {
                    params->in_smooth = MINV_INVALID;
                }
                i++;
            }
            else if (strcmp(argv[i], "--smoother-type") == 0)
            {
                if (strcmp(argv[i+1], "self") == 0)
                {
                    params->normal_eqn_smooth = false;
                }
                else if (strcmp(argv[i+1], "normal_eqn") == 0)
                {
                    params->normal_eqn_smooth = true;
                }
                i++;
            }
            else if (strcmp(argv[i], "--npre-smooth") == 0)
            {
                params->pre_smooths[0] = atoi(argv[i+1]);
                i++;
                while (i+1 != argc && argv[i+1][0] != '-')
                {
                    params->pre_smooths.push_back(atoi(argv[i+1]));
                    i++;
                }
            }
            else if (strcmp(argv[i], "--npost-smooth") == 0)
            {
                params->post_smooths[0] = atoi(argv[i+1]);
                i++;
                while (i+1 != argc && argv[i+1][0] != '-')
                {
                    params->post_smooths.push_back(atoi(argv[i+1]));
                    i++;
                }
            }
            else if (strcmp(argv[i], "--load-cfg") == 0)
            {
                params->load_cfg = argv[i+1];
                params->do_load = true;
                i++;
            }
            else if (strcmp(argv[i], "--do-freetest") == 0)
            {
                if (strcmp(argv[i+1], "yes") == 0)
                {
                    params->do_free = true;
                }
                else
                {
                    params->do_free = false;
                }
                i++;
            }
//#ifdef EIGEN_TEST
            else if (strcmp(argv[i], "--do-eigentest") == 0)
            {
                if (strcmp(argv[i+1], "yes") == 0)
                {
                    params->do_eigentest = true;
                }
                else
                {
                    params->do_eigentest = false;
                }
                i++;
            }
            else if (strcmp(argv[i], "--n-ev") == 0)
            {
                if (strcmp(argv[i+1], "all") == 0)
                {
                    params->set_eigen = -1;
                    params->set_cv = -1;
                }
                else
                {
                    params->set_eigen = atoi(argv[i+1]);
                }
                i++;
            }
            else if (strcmp(argv[i], "--n-cv") == 0)
            {
                params->set_cv = atoi(argv[i+1]);
                i++;
            }
//#endif // EIGENTEST
            else
            {
                display_usage();
                return 0;
            }
        }
    }
    
    return 1; 
}

// Set defaults.
mg_input_params::mg_input_params()
{
    ////// TOP LEVEL PROPERTIES ////////
    
    
    // Describe the staggered fermions
    // double
    mass = 0.01; // Can be overridden on command line with --mass 
    
    // What operator are we using for the solve? (Laplace is free only.) Can set on command line.
    // op_type
    opt = STAGGERED; // STAGGERED, LAPLACE, LAPLACE_NC2, G5_STAGGERED
    
    // L_x = L_y = Dimension for a lattice.
    // int
    lattice_size_x = 32; // Can be set on command line with --lattice-size. 
    lattice_size_y = 32; 

    
    /////// OUTER SOLVER ////////
    
    
    // Solver algorithm to use.
    //outer_solver
    out_solve = OUTER_GCR; 
    
    // Target relative residual for outer solver.
    //double
    outer_precision = 5e-7; 
    
    // Are we using restarts for the otuer solver?
    //bool
    outer_restart = true;
    
    // How frequently are we restarting, if we are?
    //int
    outer_restart_freq = 64; 
    
    // What's the maximum number of outer iterations before we give up?
    //int
    outer_max_iter = 100000;
    
    
    /////// NULL VECTOR GENERATION ///////////////
    
    
    // What operator are we using for null vector generation? Can set on command line.
    //op_type
    nvec_params.opt_null = STAGGERED;
    
    // Null vector generation
    
    // If GEN_NULL_VECTOR isn't defined:
    //   1 for just const vector, 2 for const + even/odd vector, 4 for each corner
    //    of the hypercube.
    // If GEN_NULL_VECTOR is defined and eo = 0:
    //   Generate "n_null_vector" null vectors which are block orthogonalized.
    // IF GEN_NULL_VECTOR is defined and eo = 1:
    //   Generate "n_null_vector" null vectors, partition into even and odd.
    //    Total number of null vectors is 2*VECTOR_COUNT. 
    // IF GEN_NULL_VECTOR is defined and eo = 3:
    //   Generate "n_null_vector" null vectors, partition into four corners.
    //    Total number of null vectors is 4*VECTOR_COUNT. 
    //int
    nvec_params.n_null_vector = 4; // Note: Gets multiplied by 2 for LAPLACE_NC2 test.
                           // Can be overriden on command line with --nvec
    
    // How are we generating null vectors?
    //mg_null_gen_type
    nvec_params.null_gen = MINV_BICGSTAB; // MINV_GCR, MINV_CG, MINV_MINRES, MINV_BICGSTAB_L
    
    // Relative precision we solve the residual equation to, per level.
    //vector<double>
    nvec_params.null_precisions.push_back(5e-5); // Can be overriden on command line with --null-precision
    
    // Maximum number of iterations for null vector generation, per level.
    //vector<int>
    nvec_params.null_max_iters.push_back(500); // Can be overriden on command line with --null-max-iter
    
    // Do we restart our null space solver?
    //bool
    nvec_params.null_restart = false; 
    
    // if we are restarting, what's the frequency?
    //int
    nvec_params.null_restart_freq = 8; 
    
    // What l to use if we're using BiCGStab-l for null vector generation.
    //int
    nvec_params.null_bicgstab_l = 4; 
    
    // What mass do we use for null vector generation?
    //double
    nvec_params.null_mass = 1e-4; 
    
    // How many ways do we split the null vectors up?
    //int
    nvec_params.null_partitions = 2; 
    
    // blocking strategy
    //blocking_strategy
    nvec_params.bstrat = BLOCK_EO; // BLOCK_NONE, BLOCK_EO, BLOCK_CORNER
    
    // Do we globally orthogonalize null vectors both against previous null vectors and their conjugate?
    //bool
    nvec_params.do_global_ortho_conj = false;
    
    // Do we split null vectors into even/odd, then orthogonalize, or do we orthogonalize first?
    //bool
    nvec_params.do_ortho_eo = false; 
    
    // Do we just generate null vectors with arpack instead?
    nvec_use_eigen = false; 
    
    
    //////// MG PROPERTIES ////////////
    
    // vector of block sizes.
    //vector<int>
    blocksizes.push_back(4); 
    
    // Multigrid information. 
    //int
    n_refine = 1; // 1 = two level V cycle, 2 = three level V cycle, etc. 
                      // Can be set on command line with --nrefine
    
    // Inner solver.
    //mg_multilevel_type
    mlevel_type = MLEVEL_SMOOTH; // V cycle = MLEVEL_SMOOTH, K cycle = MLEVEL_RECURSIVE   
    
    // MG type---do I use D_coarse everywhere, or D^\dag_coarse D_coarse?
    //bool
    normal_eqn_mg = false; 
    
    
    //////////// MG SMOOTHER //////////
    
    
    // What solver do we use to smooth?
    //inner_solver
    in_smooth = MINV_GCR; // MINV_BICGSTAB, MINV_BICGSTAB_L, MINV_CG, MINV_CR, MINV_MINRES, MINV_INVALID for none. 
    
    // What relaxation parameter do we use (for MINRES)
    //double
    omega_smooth = 0.67; // for MINRES only. 
    
    // How many pre-smoothing iterations do we perform?
    //vector<int>
    pre_smooths.push_back(6);
    
    // How many post-smoothing iterations do we perform?
    //vector<int> 
    post_smooths.push_back(6);
    
    // Do we smooth with the normal equations?
    //bool
    normal_eqn_smooth = false; // do we smooth with the normal equations?
    
    
    ////////// COARSE SOLVE ///////////
    
    
    // What solver do we use for coarse solvers?
    //inner_solver
    in_solve = GCR; //CR; //GCR; 
    
    // What inner precision do we solve to?
    //vector<double>
    inner_precisions.push_back(1e-2);
    
    // What is the restart frequency of our inner solve?
    //int
    inner_restart = 64;
    
    // What is te maximum number of iterations for our inner solve?
    //int
    inner_max = 1024;
    
    
    /////////// GAUGE FIELD /////////////
    
    
    // How are we creating the gauge field? Load it, random, unit? Can set on command line.
    //gauge_create_type
    gauge_load = GAUGE_LOAD; // GAUGE_LOAD, GAUGE_UNIT, GAUGE_RANDOM
    
    // Should we do a random gauge rotation? Can set on command line.
    //bool
    do_gauge_transform = false; 
    
    // Gauge field information.
    //double
    beta = 6.0; // For random gauge field, phase angles have std.dev. 1/sqrt(beta).
                // For heatbath gauge field, corresponds to non-compact beta.
                // Can be set on command line with --beta.
    
    // Load an external cfg?
    //char*
    load_cfg = NULL;
    //bool
    do_load = false; 
    
    
    ///////////// TEST TYPES //////////////
    
    
    // Do the free field test? Overrides # null vectors, partitioning scheme, etc. 
    //bool
    do_free = false; 
    
//#ifdef EIGEN_TEST
    //bool
    do_eigentest = false;
    
    //int
    set_eigen = -1; // default for generating all eigenvalues.
    //int
    set_cv = -1; // default for generating all eigenvalues, or 
//#endif // EIGEN_TEST 
    

}


// Verbosity header.

#ifndef VERBOSITY
#define VERBOSITY


enum inversion_verbose_level
{
    VERB_PASS_THROUGH = -1, // For preconditioned verbosity only. This means to take whatever's used in verbosity. 
    VERB_NONE = 0, // Inverter prints nothing.
    VERB_SUMMARY = 1, // Inverter prints info about self, total number of iterations at end.
    VERB_RESTART_DETAIL = 2, // SUMMARY + prints residual at each restart (for restarted inverters). Equivalent to SUMMARY for non-restarted version.
    VERB_DETAIL = 3 // RESTART_DETAIL + Prints relative residual at every step. 
};

struct inversion_verbose_struct
{
    inversion_verbose_level verbosity; // How verbose to be.
    std::string verb_prefix; // A string to prefix all verbosity printouts with.
    inversion_verbose_level precond_verbosity; // How verbose the preconditioner should be. 
    std::string precond_verb_prefix; // Prefix for preconditioner. 
};

// Shuffles preconditioned values into current values. 
void shuffle_verbosity_precond(inversion_verbose_struct* verb_new, inversion_verbose_struct* verb);

// Shuffles values for restarted solves.
void shuffle_verbosity_restart(inversion_verbose_struct* verb_new, inversion_verbose_struct* verb);

// Properly prints relative residual.
void print_verbosity_resid(inversion_verbose_struct* verb, std::string alg, int iter, double relres); 

// Properly prints summary at end of inversion.
void print_verbosity_summary(inversion_verbose_struct* verb, std::string alg, bool success, int iter, double relres); 

// Properly prints summary at restart.
void print_verbosity_restart(inversion_verbose_struct* verb, std::string alg, int iter, double relres); 

#endif // VERBOSITY
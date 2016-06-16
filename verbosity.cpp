
#include <iostream>
#include <string>
#include "verbosity.h"

// Properly prints relative residual.
void print_verbosity_resid(inversion_verbose_struct* verb, std::string alg, int iter, double relres)
{
    if (verb != 0)
    {
        if (verb->verbosity == VERB_DETAIL)
        {
            std::cout << verb->verb_prefix << alg << " Iter " << iter << " RelRes " << relres << "\n";
        }
    }
    
    return; 
}

// Properly prints summary at end of inversion.
void print_verbosity_summary(inversion_verbose_struct* verb, std::string alg, bool success, int iter, double relres)
{
    if (verb != 0)
    {
        if (verb->verbosity == VERB_SUMMARY || verb->verbosity == VERB_RESTART_DETAIL || verb->verbosity == VERB_DETAIL)
        {
            std::cout << verb->verb_prefix << alg << " Success " << (success ? "Y" : "N") << " Iter " << iter << " RelRes " << relres << "\n";
        }
    }
    
    return; 
}

// Properly prints summary at restart.
void print_verbosity_restart(inversion_verbose_struct* verb, std::string alg, int iter, double relres)
{
    if (verb != 0)
    {
        if (verb->verbosity == VERB_RESTART_DETAIL || verb->verbosity == VERB_DETAIL)
        {
            std::cout << verb->verb_prefix << alg << " Iter " << iter << " RelRes " << relres << "\n";
        }
    }
}
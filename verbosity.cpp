
#include <iostream>
#include <string>
#include "verbosity.h"

// Shuffles preconditioned values into current values. 
void shuffle_verbosity_precond(inversion_verbose_struct* verb_new, inversion_verbose_struct* verb)
{
    verb_new->verbosity = (verb->verbosity == VERB_RESTART_DETAIL || verb->verbosity == VERB_SUMMARY) ? VERB_NONE : verb->verbosity; 
    verb_new->verb_prefix = verb->verb_prefix;
    verb_new->precond_verbosity = verb->precond_verbosity;
    verb_new->precond_verb_prefix = verb->precond_verb_prefix; 
}

// Shuffles values for restarted solves.
void shuffle_verbosity_restart(inversion_verbose_struct* verb_new, inversion_verbose_struct* verb)
{
    if (verb_new->precond_verbosity == VERB_PASS_THROUGH)
    {
        verb_new->verbosity = verb->verbosity;
    }
    else
    {
        verb_new->verbosity = verb->precond_verbosity;
    }
    verb_new->verb_prefix = verb->precond_verb_prefix;
    verb_new->precond_verbosity = VERB_NONE;
    
}

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
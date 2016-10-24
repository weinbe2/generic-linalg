
#ifndef MG
#define MG

enum inner_solver
{
    NONE = 0,
    MINRES = 1,
    CG = 2,
    GCR = 3,
    BICGSTAB = 4,
    CR = 5,
    BICGSTAB_L = 6,
};

enum outer_solver
{
    OUTER_GCR = 0, // VPGCR
    OUTER_CG = 1, // FPCG
    OUTER_BICGSTAB = 2 //PBICGSTAB
}; 


#endif // MG

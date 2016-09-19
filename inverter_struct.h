// Wed Aug 24 15:46:50 EDT 2016
// Evan S Weinberg
// Include file for the inverter info struct.

#ifndef ESW_INVERTER_STRUCT
#define ESW_INVERTER_STRUCT

// Struct that contains information that
// all matrix functions return.
// Return for various matrix functions.
struct inversion_info
{
  double resSq; // squared residual.
  int iter; // number of iterations.
  bool success; // did we reach residual?
  std::string name; // name of algorithm.
  int ops_count; // how many times was the matrix op called?
    
  double* resSqmrhs; // squared residual. Only assigned for multirhs solves. 
  int n_rhs; 
  
  // Constructor
  inversion_info() : resSq(0.0), iter(0), success(false), name(""), ops_count(0), resSqmrhs(0), n_rhs(0)
  { }
  
  // Destructor.
  ~inversion_info()
  {
    if (resSqmrhs != 0)
    {
      delete[] resSqmrhs;
      resSqmrhs = 0;
    }
  }
  
  // Copy constructor.
  inversion_info(const inversion_info &obj)
  {
    resSq = obj.resSq;
    iter = obj.iter;
    success = obj.success;
    name = obj.name;
    ops_count = obj.ops_count;
    
    n_rhs = obj.n_rhs;
    if (resSqmrhs != 0)
    {
      resSqmrhs = new double[n_rhs];
      for (int i = 0; i < n_rhs; i++)
      {
        resSqmrhs[i] = obj.resSqmrhs[i];
      }
    }
    else
    {
      resSqmrhs = 0;
    }
  }
};

#endif
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
};

#endif
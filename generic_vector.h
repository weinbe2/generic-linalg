// Mon Jan 11 16:53:20 EST 2016
// Evan S Weinberg
// Header file for templated vector operations.

#include <complex>

// Zeroes a vector.
template<typename T> inline void zero(T* v1, int size)
{
   int i;
   for (i = 0; i < size; i++)
   {
      v1[i] = 0.0;
   }
}

// Zeroes a complex vector. 
template <typename T> inline void zero(std::complex<T>* v1, int size)
{
  // Initialize.
  int i;
  
  for (i = 0; i < size; i++)
  {
    v1[i] = 0.0;
  }
  
}


// Computes v1 dot v2.
template<typename T> inline T dot(T* v1, T* v2, int size)
{
  // Initialize.
  int i;
  T res = (T)0.0;
  
  for (i = 0; i < size; i++)
  {
    res = res + v1[i]*v2[i];
  }
  
  return res;
}

// computes conj(v1) dot v2.
template <typename T> inline std::complex<T> dot(std::complex<T>* v1, std::complex<T>* v2, int size)
{
  // Initialize.
  int i;
  std::complex<T> res = (T)0.0;
  
  for (i = 0; i < size; i++)
  {
    res = res + conj(v1[i])*v2[i];
  }
  
  return res;
}

template <typename T>
inline T norm2sq(T* v1,  int size)
{
  // Initialize.
  int i;
  T res = (T)0.0;
  
  for (i = 0; i < size; i++)
  {
    res = res + v1[i]*v1[i];
  }
  
  return res;
}

template <typename T>
inline T norm2sq(std::complex<T>* v1, int size)
{
  // Initialize.
  int i;
  T res = (T)0.0;
  
  for (i = 0; i < size; i++)
  {
    res = res + real(conj(v1[i])*v1[i]);
  }
  
  return res;
}


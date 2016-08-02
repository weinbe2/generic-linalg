// This file defines a trait to homogenize functions that could take real or complex numbers. 
// For ex, for a double or complex<double>, it will return double for both as the internal type.
// It will also correctly return the real part. 
// Ex: 
//
// double v1 = 1.0;
// std::complex<double> v2 = std::complex<double>(2.0,1.0);
// RealType<std::complex<double> >::Type is double.
// RealType<typeof(v2)>::Type is double.
// double v3 = RealType<typeof(v2)>::Real(v2) would return 2.

#ifndef GENERIC_TRAITS
#define GENERIC_TRAITS

#include <complex>

// Traits demo based on https://erdani.com/publications/traits.html

// Most general case not implemented.
template <typename T> struct RealType;

// For real types.
template <>
struct RealType<float>
{
   typedef float Type;
   static Type Real(Type in)
   {
      return in;
   }
};

template <>
struct RealType<double>
{
   typedef double Type;
   static Type Real(Type in)
   {
      return in;
   }
};


// For complex types.
template <>
struct RealType <std::complex<float> >
{
   typedef float Type;
   static Type Real(std::complex<Type> in)
   {
      return in.real();
   }
};

template <>
struct RealType <std::complex<double> >
{
   typedef double Type;
   static Type Real(std::complex<Type> in)
   {
      return in.real();
   }
};

#endif

#ifndef PSCF_CUDA_COMPLEX_H
#define PSCF_CUDA_COMPLEX_H

/*
* PSCF - Polymer Self-Consistent Field
*
* Copyright 2015 - 2026, The Regents of the University of Minnesota
* Distributed under the terms of the GNU General Public License.
*/

#include <pscf/backend/cuda/cudaTypes.h>
#include <pscf/math/arithmetic.h>
#include <complex>

namespace Pscf {

   /**
   * \defgroup Pscf_Cuda_Complex_Module Complex Arithmetic (GPU)
   *
   * Complex arithmetic functions using the complex type cudaComplex
   * used in GPU code that interfaces with cufft.
   *
   * \ingroup Pscf_Cuda_Module
   */

   /*
   * Types cudaComplex and cudaReal are defined in pscf/cuda/cudaTypes.h
   * (in the Prdc::Cuda namespace) as aliases for cufft complex and 
   * real types. They may be either single or double precision.
   */

   // Real and imaginary components

   /**
   * Return the real part of a complex number.
   *
   * \ingroup Pscf_Cuda_Complex_Module
   *
   * \param a complex argument (input)
   */
   inline 
   cudaReal real(cudaComplex const& a)
   {  return a.x; }

   /**
   * Return the imaginary part of a complex number.
   *
   * \ingroup Pscf_Cuda_Complex_Module
   *
   * \param a complex argument (input)
   */
   inline 
   cudaReal imag(cudaComplex const& a)
   {  return a.y; }

   // Absolute magnitude

   /**
   * Return absolute magnitude of a complex number.
   *
   * \ingroup Pscf_Cuda_Complex_Module
   *
   * \param a complex argument (in)
   */
   inline 
   cudaReal abs(cudaComplex const& a)
   {  return sqrt(a.x * a.x + a.y * a.y); }

   /**
   * Return square of absolute magnitude of a complex number.
   *
   * \ingroup Pscf_Cuda_Complex_Module
   *
   * \param a complex argument (in)
   */
   inline 
   cudaReal absSq(cudaComplex const& a)
   {  return (a.x * a.x + a.y * a.y); }

   // Complex Conjugation

   /**
   * Compute complex conjugate, z = a^*.
   *
   * \ingroup Pscf_Cuda_Complex_Module
   *
   * \param z complex conjugate of argument (out)
   * \param a complex argument (in)
   */
   inline
   void conj(cudaComplex& z, cudaComplex const& a)
   {
      z.x = a.x;
      z.y = -a.y;
   }

   /**
   * In-place complex conjugation of a complex number, a = a^* .
   *
   * \ingroup Pscf_Cuda_Complex_Module
   *
   * \param a argument (in) and complex conjugate (out)
   */
   inline
   void conj(cudaComplex& a)
   {
      a.x = a.x;
      a.y = -a.y;
   }

   // Assignment

   /**
   * Create a complex number from real and imaginary parts, z = a + ib.
   *
   * \ingroup Pscf_Cuda_Complex_Module
   *
   * \param z complex (out)
   * \param a real part (in)
   * \param b imaginary part (in)
   */
   inline
   void assign(cudaComplex& z, cudaReal const& a, cudaReal const& b)
   {
      z.x = a;
      z.y = b;
   }

   /**
   * Assign a real input to a complex variable.
   *
   * \ingroup Pscf_Cuda_Complex_Module
   *
   * \param z complex (out)
   * \param a real (in)
   */
   inline
   void assign(cudaComplex& z, cudaReal const& a)
   {
      z.x = a;
      z.y = 0;
   }

   /**
   * Assign a complex input to a complex variable.
   *
   * \ingroup Pscf_Cuda_Complex_Module
   *
   * \param z complex (out)
   * \param a complex (in)
   */
   inline
   void assign(cudaComplex& z, cudaComplex const& a)
   {
      z.x = a.x;
      z.y = a.y;
   }

   /**
   * Assign a std::complex input to a complex variable.
   *
   * \ingroup Pscf_Cuda_Complex_Module
   *
   * \param z complex (out)
   * \param a std::complex (in)
   */
   inline
   void assign(cudaComplex & z, 
               std::complex<cudaReal> const& a)
   {
      z.x = a.real();
      z.y = a.imag();
   }

   /**
   * Assign a complex input to a std::complex variable.
   *
   * \ingroup Pscf_Cuda_Complex_Module
   *
   * \param z std::complex (out)
   * \param a complex (in)
   */
   inline
   void assign(std::complex<cudaReal> & z, 
               cudaComplex const& a)
   {  z = std::complex<cudaReal>(a.x, a.y); }

   // Addition

   /**
   * Addition of two complex numbers, z = a + b.
   *
   * \ingroup Pscf_Cuda_Complex_Module
   *
   * \param z complex sum (out)
   * \param a complex summand (in)
   * \param b complex summand (in)
   */
   inline
   void add(cudaComplex& z, cudaComplex const& a, 
            cudaComplex const& b)
   {
      z.x = a.x + b.x;
      z.y = a.y + b.y;
   }

   /**
   * Addition of a complex and real number, z = a + b.
   *
   * \ingroup Pscf_Cuda_Complex_Module
   *
   * \param z complex sum (out)
   * \param a complex summand (in)
   * \param b real summand (in)
   */
   inline
   void add(cudaComplex& z, cudaComplex const& a, 
            cudaReal const& b)
   {
      z.x = a.x + b;
      z.y = a.y;
   }

   /**
   * In place addition of complex numbers, a += b.
   *
   * \ingroup Pscf_Cuda_Complex_Module
   *
   * \param a complex summand (in) and sum (out)
   * \param b complex summand (in)
   */
   inline
   void addEq(cudaComplex& a, cudaComplex const& b)
   {
      a.x += b.x;
      a.y += b.y;
   }

   /**
   * In place addition of a complex and real number, a += b.
   *
   * \ingroup Pscf_Cuda_Complex_Module
   *
   * \param a complex summand (in) and sum (out)
   * \param b real summand (in)
   */
   inline
   void addEq(cudaComplex& a, cudaReal const& b)
   {
      a.x += b;
   }

   // Subtraction

   /**
   * Subtraction of two complex numbers, z = a - b.
   *
   * \ingroup Pscf_Cuda_Complex_Module
   *
   * \param z complex difference (out)
   * \param a complex 1st argument (in)
   * \param b complex 2nd argument (in)
   */
   inline
   void sub(cudaComplex& z, cudaComplex const& a, 
            cudaComplex const& b)
   {
      z.x = a.x - b.x;
      z.y = a.y - b.y;
   }

   /**
   * Subtraction of a real number from a complex number, z = a - b.
   *
   * \ingroup Pscf_Cuda_Complex_Module
   *
   * \param z complex difference (out)
   * \param a complex 1st argument (in)
   * \param b real 2nd argument (in)
   */
   inline
   void sub(cudaComplex& z, cudaComplex const& a, 
            cudaReal const& b)
   {
      z.x = a.x - b;
      z.y = a.y;
   }

   /**
   * In place subtraction of two complex numbers, a -= b.
   *
   * \ingroup Pscf_Cuda_Complex_Module
   *
   * \param a complex argument (in) and difference (out)
   * \param b complex argument (in)
   */
   inline
   void subEq(cudaComplex & a, cudaComplex const& b)
   {
      a.x -= b.x;
      a.y -= b.y;
   }

   /**
   * In place subtraction of real number from a complex number, a -= b.
   *
   * \ingroup Pscf_Cuda_Complex_Module
   *
   * \param a complex argument (in) and difference (out)
   * \param b real argument (in)
   */
   inline
   void subEq(cudaComplex & a, cudaReal const& b)
   {
      a.x -= b;
   }

   /**
   * Return square of the absolute magnitude of a complex difference.
   *
   * This function returns |a-b|^2 for complex a and b.
   *
   * \ingroup Pscf_Cuda_Complex_Module
   *
   * \param a complex 1st argument (in)
   * \param b complex 2nd argument (in)
   */
   inline 
   cudaReal absSqDiff(cudaComplex const& a, 
                      cudaComplex const& b)
   {
      cudaComplex z;
      sub(z, a, b);
      return absSq(z);
   }

   // Multiplication

   /**
   * Multiplication of two complex numbers, z = a * b.
   *
   * \ingroup Pscf_Cuda_Complex_Module
   *
   * \param z complex product (out)
   * \param a complex factor (in)
   * \param b complex factor (in)
   */
   inline
   void mul(cudaComplex& z, cudaComplex const& a, 
            cudaComplex const& b)
   {
      z.x = a.x * b.x - a.y * b.y;
      z.y = a.y * b.x + a.x * b.y;
   }

   /**
   * Multiplication of complex and real numbers, z = a * b.
   *
   * \ingroup Pscf_Cuda_Complex_Module
   *
   * \param z complex product (out)
   * \param a complex factor (in)
   * \param b real factor (in)
   */
   inline
   void mul(cudaComplex& z, cudaComplex const& a, 
            cudaReal const& b)
   {
      z.x = a.x*b;
      z.y = a.y*b;
   }

   /**
   * In place multiplication of two complex number, a *= b.
   *
   * \ingroup Pscf_Cuda_Complex_Module
   *
   * \param a complex factor (in) and product (out)
   * \param b complex factor (in)
   */
   inline
   void mulEq(cudaComplex & a, cudaComplex const& b)
   {
      cudaReal a0;
      a0   = a.x * b.x - a.y * b.y;
      a.y = a.y * b.x + a.x * b.y;
      a.x = a0;
   }

   /**
   * In place multiplication of a complex and real number, a *= b.
   *
   * \ingroup Pscf_Cuda_Complex_Module
   *
   * \param a complex factor (in) and product (out)
   * \param b real factor (in)
   */
   inline
   void mulEq(cudaComplex & a, cudaReal const& b)
   {
      a.x *= b;
      a.y *= b;
   }

   /**
   * Compute complex square of a complex number, z = a * a.
   *
   * \ingroup Pscf_Cuda_Complex_Module
   *
   * \param z complex product (out)
   * \param a complex factor (in)
   */
   inline
   void square(cudaComplex& z, cudaComplex const& a)
   {
      z.x = a.x * a.x - a.y * a.y;
      z.y = 2.0 * a.y * a.x;
   }

   // Division

   /**
   * Division of two complex numbers, z = a / b .
   *
   * \ingroup Pscf_Cuda_Complex_Module
   *
   * \param z complex ratio (out)
   * \param a complex numerator (in)
   * \param b complex denominator (in)
   */
   inline
   void div(cudaComplex& z, cudaComplex const& a, 
            cudaComplex const& b)
   {
      cudaReal bSq = b.x * b.x + b.y * b.y;
      z.x = (a.x * b.x + a.y * b.y) / bSq;
      z.y = (a.y * b.x - a.x * b.y) / bSq;
   }

   /**
   * Division of a complex number by a real number, z = a / b .
   *
   * \ingroup Pscf_Cuda_Complex_Module
   *
   * \param z complex ratio (out)
   * \param a complex numerator (in)
   * \param b real denominator (in)
   */
   inline
   void div(cudaComplex& z, cudaComplex const& a, 
            cudaReal const& b)
   {
      z.x = a.x / b;
      z.y = a.y / b;
   }

   /**
   * In place division of two complex number, a /= b.
   *
   * \ingroup Pscf_Cuda_Complex_Module
   *
   * \param a complex numerator (in) and ratio (out)
   * \param b complex denominator (in)
   */
   inline
   void divEq(cudaComplex & a, cudaComplex const & b)
   {
      cudaReal bSq = b.x * b.x + b.y * b.y;
      cudaReal a0 = (a.x * b.x + a.y * b.y)/bSq;
      a.y = (a.y * b.x - a.x * b.y)/bSq;
      a.x = a0;
   }

   /**
   * In place division of a complex number by a real number, a /= b.
   *
   * \ingroup Pscf_Cuda_Complex_Module
   *
   * \param a complex numerator (in) and ratio (out)
   * \param b real denominator (in)
   */
   inline
   void divEq(cudaComplex & a, cudaReal const& b)
   {
      a.x /= b;
      a.y /= b;
   }

   // Inversion

   /**
   * Inversion of a complex number, z = 1 / a .
   *
   * \ingroup Pscf_Cuda_Complex_Module
   *
   * \param z inverse (out)
   * \param a argument (in)
   */
   inline
   void inverse(cudaComplex& z, 
                cudaComplex const & a)
   {
      cudaReal aSq = a.x * a.x + a.y * a.y;
      z.x =   a.x/aSq;
      z.y = - a.y/aSq;
   }

   // Exponentiation and logarithm

   /**
   * Exponentation of a ffts_complex variable, z = exp(a).
   *
   * \ingroup Pscf_Cuda_Complex_Module
   *
   * \param z exponent (out)
   * \param a argument (in)
   */
   inline
   void assignExp(cudaComplex & z, 
                  cudaComplex const & a)
   {  
      std::complex<cudaReal> arg 
                  = std::complex<cudaReal>(a.x, a.y); 
      std::complex<cudaReal> result = std::exp(arg);
      z.x = result.real();
      z.y = result.imag();
   }

   /**
   * Logarithm of a cudaComplex variable, z = log(a).
   *
   * \ingroup Pscf_Cuda_Complex_Module
   *
   * \param z logarithm (out)
   * \param a argument (in)
   */
   inline
   void assignLog(cudaComplex & z, 
                  cudaComplex const & a)
   {  
      std::complex<cudaReal> arg 
                 = std::complex<cudaReal>(a.x, a.y); 
      std::complex<cudaReal> result = std::log(arg);
      z.x = result.real();
      z.y = result.imag();
   }

   // Pseudo-constructor

   /*
   * Pseudo-constructor function for cudaComplex.
   *
   * \ingroup Pscf_Cuda_Complex_Module
   *
   * \param x  real part (in)
   * \param y  imaginary part (in)
   */
   __host__ __device__ static inline   
   cudaComplex makeComplex(cudaReal x, cudaReal y)
   {
      cudaComplex result;
      result.x = x;
      result.y = y;
      return result;
   }

} // namespace Pscf
#endif

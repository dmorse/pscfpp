#ifndef PSCF_CUDA_COMPLEX_H
#define PSCF_CUDA_COMPLEX_H

/*
* PSCF Package - Polymer Self-Consistent Field 
*
* Copyright 2016 - 2023, The Regents of the University of Minnesota
* Distributed under the terms of the GNU General Public License.
*/

#include <pscf/cuda/GpuTypes.h>
#include <complex>

namespace Pscf {
namespace Cuda {

   /*
   * Types Cuda::Complex and Cuda::Real are defined in pscf/cpu/GpuTypes.h
   * as aliases for cufft complex and real types.
   */

   // Real and imaginary components

   /**
   * Return the real part of a complex number.
   *
   * \ingroup Pscf_Cuda_Complex_Module
   *
   * \param a complex argument (input)
   */
   inline Real real(Complex const& a)
   {  return a.x; }

   /**
   * Return the imaginary part of a complex number.
   *
   * \ingroup Pscf_Cuda_Complex_Module
   *
   * \param a complex argument (input)
   */
   inline Real imag(Complex const& a)
   {  return a.y; }

   // Absolute magnitude

   /**
   * Return absolute magnitude of a complex number.
   *
   * \ingroup Pscf_Cuda_Complex_Module
   *
   * \param a complex argument (in)
   */
   inline Real abs(Complex const& a)
   {  return sqrt(a.x * a.x + a.y * a.y); }

   /**
   * Return square of absolute magnitude of a complex number.
   *
   * \ingroup Pscf_Cuda_Complex_Module
   *
   * \param a complex argument (in)
   */
   inline Real absSq(Complex const& a)
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
   void conj(Complex& z, Complex const& a)
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
   void conj(Complex& a)
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
   void assign(Complex& z, Real const& a, Real const& b)
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
   void assign(Complex& z, Real const& a)
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
   void assign(Complex& z, Complex const& a)
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
   void assign(Complex & z, std::complex<Real> const& a)
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
   void assign(std::complex<Real> & z, Complex const& a)
   {  z = std::complex<Real>(a.x, a.y); }

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
   void add(Complex& z, Complex const& a, Complex const& b)
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
   void add(Complex& z, Complex const& a, Real const& b)
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
   void addEq(Complex& a, Complex const& b)
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
   void addEq(Complex& a, Real const& b)
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
   void sub(Complex& z, Complex const& a, Complex const& b)
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
   void sub(Complex& z, Complex const& a, Real const& b)
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
   void subEq(Complex & a, Complex const& b)
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
   void subEq(Complex & a, Real const& b)
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
   inline Real absSqDiff(Complex const& a, Complex const& b)
   {
      Complex z;
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
   void mul(Complex& z, Complex const& a, Complex const& b)
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
   void mul(Complex& z, Complex const& a, Real const& b)
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
   void mulEq(Complex & a, Complex const& b)
   {
      Real a0;
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
   void mulEq(Complex & a, Real const& b)
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
   void square(Complex& z, Complex const& a)
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
   void div(Complex& z, Complex const& a, Complex const& b)
   {
      Real bSq = b.x * b.x + b.y * b.y;
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
   void div(Complex& z, Complex const& a, Real const& b)
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
   void divEq(Complex & a, Complex const & b)
   {
      Real bSq = b.x * b.x + b.y * b.y;
      Real a0 = (a.x * b.x + a.y * b.y)/bSq;
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
   void divEq(Complex & a, Real const& b)
   {
      a.x /= b;
      a.y /= b;
   }

} // Pscf::Cuda
} // Pscf
#endif

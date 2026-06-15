#ifndef PSCF_CPU_COMPLEX_H
#define PSCF_CPU_COMPLEX_H

/*
* PSCF - Polymer Self-Consistent Field
*
* Copyright 2015 - 2025, The Regents of the University of Minnesota
* Distributed under the terms of the GNU General Public License.
*/

#include <pscf/math/arithmetic.h>  // built-in and standard lib types

#include <fftw3.h>
#include <complex>
#include <iostream>

namespace Pscf {

   /**
   * \defgroup Pscf_Cpu_Complex_Module Complex Arithmetic (CPU)
   *
   * Complex arithmetic functions using the complex type fftw_complex 
   * defined by the FFTW fast Fourier transform library.
   *
   * \ingroup Pscf_Cpu_Module
   */

   // Real and imaginary components of fftw_complex numbers

   /**
   * Return the real part of an fftw_complex number.
   *
   * \ingroup Pscf_Cpu_Complex_Module
   *
   * \param a complex argument (in)
   */
   inline 
   double real(fftw_complex const & a)
   {  return a[0]; }

   /**
   * Return the imaginary part of an fftw_complex number.
   *
   * \ingroup Pscf_Cpu_Complex_Module
   *
   * \param a complex argument (in)
   */
   inline 
   double imag(fftw_complex const & a)
   {  return a[1]; }

   // Absolute magnitude

   /**
   * Return absolute magnitude of an fftw_complex number.
   *
   * \ingroup Pscf_Cpu_Complex_Module
   *
   * \param a complex argument (in)
   */
   inline 
   double abs(fftw_complex const & a)
   {  return sqrt(a[0] * a[0] + a[1] * a[1]); }

   /**
   * Return square of absolute magnitude of an fftw_complex number.
   *
   * \ingroup Pscf_Cpu_Complex_Module
   *
   * \param a complex argument (in)
   */
   inline 
   double absSq(fftw_complex const & a)
   {  return (a[0] * a[0] + a[1] * a[1]); }

   // Complex Conjugation

   /**
   * Complex conjugate of an fftw_complex, z = a^*.
   *
   * \ingroup Pscf_Cpu_Complex_Module
   *
   * \param z complex conjugate of argument (out)
   * \param a complex argument (in)
   */
   inline
   void conj(fftw_complex& z, fftw_complex const & a)
   {
      z[0] = a[0];
      z[1] = -a[1];
   }

   /**
   * In-place complex conjugation of an fftw_complex number, a = a^* .
   *
   * \ingroup Pscf_Cpu_Complex_Module
   *
   * \param a argument and output (in/out)
   */
   inline
   void conj(fftw_complex& a)
   {
      a[0] = a[0];
      a[1] = -a[1];
   }

   // Assignment

   /**
   * Create an fftw_complex from real and imaginary parts, z = a + ib.
   *
   * \ingroup Pscf_Cpu_Complex_Module
   *
   * \param z complex (out)
   * \param a real part (in)
   * \param b imaginary part (in)
   */
   inline
   void assign(fftw_complex& z, double const & a, double const & b)
   {
      z[0] = a;
      z[1] = b;
   }

   /**
   * Assign a real input to an fftw_complex variable.
   *
   * \ingroup Pscf_Cpu_Complex_Module
   *
   * \param z complex (out)
   * \param a real (in)
   */
   inline
   void assign(fftw_complex& z, double const & a)
   {
      z[0] = a;
      z[1] = 0.0;
   }

   /**
   * Assign an fftw_complex input to an fftw_complex variable.
   *
   * \ingroup Pscf_Cpu_Complex_Module
   *
   * \param z complex (out)
   * \param a complex (in)
   */
   inline
   void assign(fftw_complex& z, fftw_complex const & a)
   {
      z[0] = a[0];
      z[1] = a[1];
   }

   /**
   * Assign a std::complex input to an fftw_complex variable.
   *
   * \ingroup Pscf_Cpu_Complex_Module
   *
   * \param z complex (out)
   * \param a std::complex (in)
   */
   inline
   void assign(fftw_complex & z, std::complex<double> const& a)
   {
      z[0] = a.real();
      z[1] = a.imag();
   }

   /**
   * Assign an fftw_complex input to a std::complex variable.
   *
   * \ingroup Pscf_Cpu_Complex_Module
   *
   * \param z std::complex (out)
   * \param a complex (in)
   */
   inline
   void assign(std::complex<double> & z, fftw_complex const& a)
   {  z = std::complex<double>(a[0], a[1]); }

   // Addition

   /**
   * Addition of fftw_complex numbers, z = a + b.
   *
   * \ingroup Pscf_Cpu_Complex_Module
   *
   * \param z complex sum (out)
   * \param a complex summand (in)
   * \param b complex summand (in)
   */
   inline
   void add(fftw_complex& z, fftw_complex const& a, fftw_complex const& b)
   {
      z[0] = a[0] + b[0];
      z[1] = a[1] + b[1];
   }

   /**
   * Addition of a complex and real number, z = a + b.
   *
   * \ingroup Pscf_Cpu_Complex_Module
   *
   * \param z complex sum (out)
   * \param a complex summand (in)
   * \param b real summand (in)
   */
   inline
   void add(fftw_complex& z, fftw_complex const& a, double const& b)
   {
      z[0] = a[0] + b;
      z[1] = a[1];
   }

   /**
   * In-place addition of fftw_complex numbers, a += b.
   *
   * \ingroup Pscf_Cpu_Complex_Module
   *
   * \param a complex summand (in) and sum (out)
   * \param b complex summand (in)
   */
   inline
   void addEq(fftw_complex& a, fftw_complex const& b)
   {
      a[0] += b[0];
      a[1] += b[1];
   }

   /**
   * In-place addition of a complex and real number, a += b.
   *
   * \ingroup Pscf_Cpu_Complex_Module
   *
   * \param a complex summand (in) and sum (out)
   * \param b real summand (in)
   */
   inline
   void addEq(fftw_complex& a, double const& b)
   {
      a[0] += b;
   }

   // Subtraction

   /**
   * Subtraction of fftw_complex numbers, z = a - b.
   *
   * \ingroup Pscf_Cpu_Complex_Module
   *
   * \param z complex difference (out)
   * \param a complex 1st argument (in)
   * \param b complex 2nd argument (in)
   */
   inline
   void sub(fftw_complex& z, fftw_complex const& a, fftw_complex const& b)
   {
      z[0] = a[0] - b[0];
      z[1] = a[1] - b[1];
   }

   /**
   * Subtraction of a real number from an fftw_complex number, z = a - b.
   *
   * \ingroup Pscf_Cpu_Complex_Module
   *
   * \param z complex difference (out)
   * \param a complex 1st argument (in)
   * \param b real 2nd argument (in)
   */
   inline
   void sub(fftw_complex& z, fftw_complex const& a, double const& b)
   {
      z[0] = a[0] - b;
      z[1] = a[1];
   }

   /**
   * In-place subtraction of fftw_complex numbers, a -= b.
   *
   * \ingroup Pscf_Cpu_Complex_Module
   *
   * \param a complex argument (in) and difference (out)
   * \param b complex argument (in)
   */
   inline
   void subEq(fftw_complex & a, fftw_complex const& b)
   {
      a[0] -= b[0];
      a[1] -= b[1];
   }

   /**
   * In-place subtraction of real from an fftw_complex number, a -= b.
   *
   * \ingroup Pscf_Cpu_Complex_Module
   *
   * \param a complex argument (in) and difference (out)
   * \param b real argument (in)
   */
   inline
   void subEq(fftw_complex & a, double const& b)
   {
      a[0] -= b;
   }

   /**
   * Return square of the absolute magnitude of a complex difference.
   *
   * This function returns |a-b|^2 for complex a and b.
   *
   * \ingroup Pscf_Cpu_Complex_Module
   *
   * \param a complex 1st argument (in)
   * \param b complex 2nd argument (in)
   */
   inline 
   double absSqDiff(fftw_complex const& a, fftw_complex const& b)
   {
      fftw_complex z;
      sub(z, a, b);
      return absSq(z);
   }

   // Multiplication

   /**
   * Multiplication of fftw_complex numbers, z = a * b.
   *
   * \ingroup Pscf_Cpu_Complex_Module
   *
   * \param z complex product (out)
   * \param a complex factor (in)
   * \param b complex factor (in)
   */
   inline
   void mul(fftw_complex& z, fftw_complex const& a, fftw_complex const& b)
   {
      z[0] = a[0] * b[0] - a[1] * b[1];
      z[1] = a[1] * b[0] + a[0] * b[1];
   }

   /**
   * Multiplication of an fftw_complex and real number, z = a * b.
   *
   * \ingroup Pscf_Cpu_Complex_Module
   *
   * \param z complex product (out)
   * \param a complex factor (in)
   * \param b real factor (in)
   */
   inline
   void mul(fftw_complex& z, fftw_complex const& a, double const& b)
   {
      z[0] = a[0]*b;
      z[1] = a[1]*b;
   }

   /**
   * In-place multiplication of two complex number, a *= b.
   *
   * \ingroup Pscf_Cpu_Complex_Module
   *
   * \param a complex factor (in) and product (out)
   * \param b complex factor (in)
   */
   inline
   void mulEq(fftw_complex & a, fftw_complex const& b)
   {
      double a0;
      a0   = a[0] * b[0] - a[1] * b[1];
      a[1] = a[1] * b[0] + a[0] * b[1];
      a[0] = a0;
   }

   /**
   * In-place multiplication of a complex and real number, a *= b.
   *
   * \ingroup Pscf_Cpu_Complex_Module
   *
   * \param a complex factor (in) and product (out)
   * \param b real factor (in)
   */
   inline
   void mulEq(fftw_complex & a, double const& b)
   {
      a[0] *= b;
      a[1] *= b;
   }

   /**
   * Complex square of an fftw_complex number, z = a * a.
   *
   * \ingroup Pscf_Cpu_Complex_Module
   *
   * \param z complex product (out)
   * \param a complex factor (in)
   */
   inline
   void square(fftw_complex& z, fftw_complex const& a)
   {
      z[0] = a[0] * a[0] - a[1] * a[1];
      z[1] = 2.0 * a[1] * a[0];
   }

   // Division

   /**
   * Division of fftw_complex numbers, z = a / b .
   *
   * \ingroup Pscf_Cpu_Complex_Module
   *
   * \param z complex ratio (out)
   * \param a complex numerator (in)
   * \param b complex denominator (in)
   */
   inline
   void div(fftw_complex& z, fftw_complex const& a, fftw_complex const& b)
   {
      double bSq = b[0] * b[0] + b[1] * b[1];
      z[0] = (a[0] * b[0] + a[1] * b[1])/bSq;
      z[1] = (a[1] * b[0] - a[0] * b[1])/bSq;
   }

   /**
   * Division of an fftw_complex number by a real number, z = a / b .
   *
   * \ingroup Pscf_Cpu_Complex_Module
   *
   * \param z complex ratio (out)
   * \param a complex numerator (in)
   * \param b real denominator (in)
   */
   inline
   void div(fftw_complex& z, fftw_complex const& a, double const& b)
   {
      z[0] = a[0]/b;
      z[1] = a[1]/b;
   }

   /**
   * In-place division of fftw_complex numbers, a /= b.
   *
   * \ingroup Pscf_Cpu_Complex_Module
   *
   * \param a complex numerator (in) and ratio (out)
   * \param b complex denominator (in)
   */
   inline
   void divEq(fftw_complex & a, fftw_complex const & b)
   {
      double bSq = b[0] * b[0] + b[1] * b[1];
      double a0 = (a[0] * b[0] + a[1] * b[1])/bSq;
      a[1] = (a[1] * b[0] - a[0] * b[1])/bSq;
      a[0] = a0;
   }

   /**
   * In-place division of an fftw_complex number by a real number, a /= b.
   *
   * \ingroup Pscf_Cpu_Complex_Module
   *
   * \param a complex numerator (in) and ratio (out)
   * \param b real denominator (in)
   */
   inline
   void divEq(fftw_complex & a, double const & b)
   {
      a[0] /= b;
      a[1] /= b;
   }

   // Inversion

   /**
   * Inversion of an fftw_complex number, z = 1 / a .
   *
   * \ingroup Pscf_Cpu_Complex_Module
   *
   * \param z inverse (out)
   * \param a argument (in)
   */
   inline
   void inverse(fftw_complex& z, fftw_complex const & a)
   {
      double aSq = a[0] * a[0] + a[1] * a[1];
      z[0] =   a[0]/aSq;
      z[1] = - a[1]/aSq;
   }

   // Exponentiation and logarithm

   /**
   * Exponentation of a ffts_complex variable, z = exp(a).
   *
   * \ingroup Pscf_Cpu_Complex_Module
   *
   * \param z exponent (out)
   * \param a argument (in)
   */
   void assignExp(fftw_complex & z, fftw_complex const & a);

   #if 0
   {
      std::complex<double> arg = std::complex<double>(a[0], a[1]); 
      std::complex<double> result = std::exp(arg);
      z[0] = result.real();
      z[1] = result.imag();
   }
   #endif

   /**
   * Logarithm of an fftw_complex variable, z = log(a).
   *
   * \ingroup Pscf_Cpu_Complex_Module
   *
   * \param z logarithm (out)
   * \param a argument (in)
   */
   void assignLog(fftw_complex & z, fftw_complex const & a);

   #if 0
   {  
      std::complex<double> arg = std::complex<double>(a[0], a[1]); 
      std::complex<double> result = std::log(arg);
      z[0] = result.real();
      z[1] = result.imag();
   }
   #endif

   // Stream IO operators

   /**
   * Stream extraction operator for fftw_complex
   *
   * \ingroup Pscf_Cpu_Complex_Module
   *
   * \param is  input stream
   * \param z   complex number
   */
   std::istream& operator >> (std::istream& is, fftw_complex & z);

   /**
   * Stream insertion operator for fftw_complex
   *
   * \ingroup Pscf_Cpu_Complex_Module
   *
   * \param os  output stream
   * \param z  complex number
   * \return modified output stream
   */
   std::ostream& operator << (std::ostream& os, fftw_complex const & z);

} // namespace Pscf
#endif

#ifndef PRDC_CPU_COMPLEX_H
#define PRDC_CPU_COMPLEX_H

/*
* PSCF Package - Polymer Self-Consistent Field Theory
*
* Copyright 2016 - 2022, The Regents of the University of Minnesota
* Distributed under the terms of the GNU General Public License.
*/

#include <pscf/math/complex.h>

#include <fftw3.h>
#include <complex>
#include <iostream>

namespace Pscf {

   /*
   * Types Cpu::Complex and Cpu::Real are defined in prdc/cpu/types.h
   * as aliases for fftw_complex and double, respectively.
   */

   // Real and imaginary components

   /**
   * Return the real part of a complex number.
   *
   * \ingroup Prdc_Cpu_Complex_Module
   *
   * \param a complex argument (input)
   */
   template <> inline 
   double real<fftw_complex, double>(fftw_complex const& a)
   {  return a[0]; }

   /**
   * Return the imaginary part of a complex number.
   *
   * \ingroup Prdc_Cpu_Complex_Module
   *
   * \param a complex argument (input)
   */
   template <> inline 
   double imag<fftw_complex, double>(fftw_complex const& a)
   {  return a[1]; }

   // Absolute magnitude

   /**
   * Return absolute magnitude of a complex number.
   *
   * \ingroup Prdc_Cpu_Complex_Module
   *
   * \param a complex argument (in)
   */
   template <> inline 
   double abs<fftw_complex, double>(fftw_complex const& a)
   {  return sqrt(a[0] * a[0] + a[1] * a[1]); }

   /**
   * Return square of absolute magnitude of a complex number.
   *
   * \ingroup Prdc_Cpu_Complex_Module
   *
   * \param a complex argument (in)
   */
   template <> inline 
   double absSq<fftw_complex, double>(fftw_complex const& a)
   {  return (a[0] * a[0] + a[1] * a[1]); }

   // fftw_complex Conjugation

   /**
   * Compute complex conjugate, z = a^*.
   *
   * \ingroup Prdc_Cpu_Complex_Module
   *
   * \param z complex conjugate of argument (out)
   * \param a complex argument (in)
   */
   template <> inline
   void conj(fftw_complex& z, fftw_complex const& a)
   {
      z[0] = a[0];
      z[1] = -a[1];
   }

   /**
   * In-place complex conjugation of a complex number, a = a^* .
   *
   * \ingroup Prdc_Cpu_Complex_Module
   *
   * \param a argument (in) and complex conjugate (out)
   */
   template <> inline
   void conj(fftw_complex& a)
   {
      a[0] = a[0];
      a[1] = -a[1];
   }

   // Assignment

   /**
   * Create a complex number from real and imaginary parts, z = a + ib.
   *
   * \ingroup Prdc_Cpu_Complex_Module
   *
   * \param z complex (out)
   * \param a real part (in)
   * \param b imaginary part (in)
   */
   template <> inline
   void assign(fftw_complex& z, double const& a, double const& b)
   {
      z[0] = a;
      z[1] = b;
   }

   /**
   * Assign a real input to a complex variable.
   *
   * \ingroup Prdc_Cpu_Complex_Module
   *
   * \param z complex (out)
   * \param a real (in)
   */
   template <> inline
   void assign(fftw_complex& z, double const& a)
   {
      z[0] = a;
      z[1] = 0;
   }

   /**
   * Assign a complex input to a complex variable.
   *
   * \ingroup Prdc_Cpu_Complex_Module
   *
   * \param z complex (out)
   * \param a complex (in)
   */
   template <> inline
   void assign(fftw_complex& z, fftw_complex const& a)
   {
      z[0] = a[0];
      z[1] = a[1];
   }

   /**
   * Assign a std::complex input to a complex variable.
   *
   * \ingroup Prdc_Cpu_Complex_Module
   *
   * \param z complex (out)
   * \param a std::complex (in)
   */
   template <> inline
   void assign(fftw_complex & z, std::complex<double> const& a)
   {
      z[0] = a.real();
      z[1] = a.imag();
   }

   /**
   * Assign a complex input to a std::complex variable.
   *
   * \ingroup Prdc_Cpu_Complex_Module
   *
   * \param z std::complex (out)
   * \param a complex (in)
   */
   template <> inline
   void assign(std::complex<double> & z, fftw_complex const& a)
   {  z = std::complex<double>(a[0], a[1]); }

   // Addition

   /**
   * Addition of two complex numbers, z = a + b.
   *
   * \ingroup Prdc_Cpu_Complex_Module
   *
   * \param z complex sum (out)
   * \param a complex summand (in)
   * \param b complex summand (in)
   */
   template <> inline
   void add(fftw_complex& z, fftw_complex const& a, fftw_complex const& b)
   {
      z[0] = a[0] + b[0];
      z[1] = a[1] + b[1];
   }

   /**
   * Addition of a complex and real number, z = a + b.
   *
   * \ingroup Prdc_Cpu_Complex_Module
   *
   * \param z complex sum (out)
   * \param a complex summand (in)
   * \param b real summand (in)
   */
   template <> inline
   void add(fftw_complex& z, fftw_complex const& a, double const& b)
   {
      z[0] = a[0] + b;
      z[1] = a[1];
   }

   /**
   * In place addition of complex numbers, a += b.
   *
   * \ingroup Prdc_Cpu_Complex_Module
   *
   * \param a complex summand (in) and sum (out)
   * \param b complex summand (in)
   */
   template <> inline
   void addEq(fftw_complex& a, fftw_complex const& b)
   {
      a[0] += b[0];
      a[1] += b[1];
   }

   /**
   * In place addition of a complex and real number, a += b.
   *
   * \ingroup Prdc_Cpu_Complex_Module
   *
   * \param a complex summand (in) and sum (out)
   * \param b real summand (in)
   */
   template <> inline
   void addEq(fftw_complex& a, double const& b)
   {
      a[0] += b;
   }

   // Subtraction

   /**
   * Subtraction of two complex numbers, z = a - b.
   *
   * \ingroup Prdc_Cpu_Complex_Module
   *
   * \param z complex difference (out)
   * \param a complex 1st argument (in)
   * \param b complex 2nd argument (in)
   */
   template <> inline
   void sub(fftw_complex& z, fftw_complex const& a, fftw_complex const& b)
   {
      z[0] = a[0] - b[0];
      z[1] = a[1] - b[1];
   }

   /**
   * Subtraction of a real number from a complex number, z = a - b.
   *
   * \ingroup Prdc_Cpu_Complex_Module
   *
   * \param z complex difference (out)
   * \param a complex 1st argument (in)
   * \param b real 2nd argument (in)
   */
   template <> inline
   void sub(fftw_complex& z, fftw_complex const& a, double const& b)
   {
      z[0] = a[0] - b;
      z[1] = a[1];
   }

   /**
   * In place subtraction of two complex numbers, a -= b.
   *
   * \ingroup Prdc_Cpu_Complex_Module
   *
   * \param a complex argument (in) and difference (out)
   * \param b complex argument (in)
   */
   template <> inline
   void subEq(fftw_complex & a, fftw_complex const& b)
   {
      a[0] -= b[0];
      a[1] -= b[1];
   }

   /**
   * In place subtraction of real number from a complex number, a -= b.
   *
   * \ingroup Prdc_Cpu_Complex_Module
   *
   * \param a complex argument (in) and difference (out)
   * \param b real argument (in)
   */
   template <> inline
   void subEq(fftw_complex & a, double const& b)
   {
      a[0] -= b;
   }

   /**
   * Return square of the absolute magnitude of a complex difference.
   *
   * This function returns |a-b|^2 for complex a and b.
   *
   * \ingroup Prdc_Cpu_Complex_Module
   *
   * \param a complex 1st argument (in)
   * \param b complex 2nd argument (in)
   */
   template <> inline 
   double absSqDiff(fftw_complex const& a, fftw_complex const& b)
   {
      fftw_complex z;
      sub(z, a, b);
      return absSq<fftw_complex, double>(z);
   }

   // Multiplication

   /**
   * Multiplication of two complex numbers, z = a * b.
   *
   * \ingroup Prdc_Cpu_Complex_Module
   *
   * \param z complex product (out)
   * \param a complex factor (in)
   * \param b complex factor (in)
   */
   template <> inline
   void mul(fftw_complex& z, fftw_complex const& a, fftw_complex const& b)
   {
      z[0] = a[0] * b[0] - a[1] * b[1];
      z[1] = a[1] * b[0] + a[0] * b[1];
   }

   /**
   * Multiplication of complex and real numbers, z = a * b.
   *
   * \ingroup Prdc_Cpu_Complex_Module
   *
   * \param z complex product (out)
   * \param a complex factor (in)
   * \param b real factor (in)
   */
   template <> inline
   void mul(fftw_complex& z, fftw_complex const& a, double const& b)
   {
      z[0] = a[0]*b;
      z[1] = a[1]*b;
   }

   /**
   * In place multiplication of two complex number, a *= b.
   *
   * \ingroup Prdc_Cpu_Complex_Module
   *
   * \param a complex factor (in) and product (out)
   * \param b complex factor (in)
   */
   template <> inline
   void mulEq(fftw_complex & a, fftw_complex const& b)
   {
      double a0;
      a0   = a[0] * b[0] - a[1] * b[1];
      a[1] = a[1] * b[0] + a[0] * b[1];
      a[0] = a0;
   }

   /**
   * In place multiplication of a complex and real number, a *= b.
   *
   * \ingroup Prdc_Cpu_Complex_Module
   *
   * \param a complex factor (in) and product (out)
   * \param b real factor (in)
   */
   template <> inline
   void mulEq(fftw_complex & a, double const& b)
   {
      a[0] *= b;
      a[1] *= b;
   }

   /**
   * Compute complex square of a complex number, z = a * a.
   *
   * \ingroup Prdc_Cpu_Complex_Module
   *
   * \param z complex product (out)
   * \param a complex factor (in)
   */
   template <> inline
   void square(fftw_complex& z, fftw_complex const& a)
   {
      z[0] = a[0] * a[0] - a[1] * a[1];
      z[1] = 2.0 * a[1] * a[0];
   }

   // Division

   /**
   * Division of two complex numbers, z = a / b .
   *
   * \ingroup Prdc_Cpu_Complex_Module
   *
   * \param z complex ratio (out)
   * \param a complex numerator (in)
   * \param b complex denominator (in)
   */
   template <> inline
   void div(fftw_complex& z, fftw_complex const& a, fftw_complex const& b)
   {
      double bSq = b[0] * b[0] + b[1] * b[1];
      z[0] = (a[0] * b[0] + a[1] * b[1])/bSq;
      z[1] = (a[1] * b[0] - a[0] * b[1])/bSq;
   }

   /**
   * Division of a complex number by a real number, z = a / b .
   *
   * \ingroup Prdc_Cpu_Complex_Module
   *
   * \param z complex ratio (out)
   * \param a complex numerator (in)
   * \param b real denominator (in)
   */
   template <> inline
   void div(fftw_complex& z, fftw_complex const& a, double const& b)
   {
      z[0] = a[0]/b;
      z[1] = a[1]/b;
   }

   /**
   * In place division of two complex number, a /= b.
   *
   * \ingroup Prdc_Cpu_Complex_Module
   *
   * \param a complex numerator (in) and ratio (out)
   * \param b complex denominator (in)
   */
   template <> inline
   void divEq(fftw_complex & a, fftw_complex const & b)
   {
      double bSq = b[0] * b[0] + b[1] * b[1];
      double a0 = (a[0] * b[0] + a[1] * b[1])/bSq;
      a[1] = (a[1] * b[0] - a[0] * b[1])/bSq;
      a[0] = a0;
   }

   /**
   * In place division of a complex number by a real number, a /= b.
   *
   * \ingroup Prdc_Cpu_Complex_Module
   *
   * \param a complex numerator (in) and ratio (out)
   * \param b real denominator (in)
   */
   template <> inline
   void divEq(fftw_complex & a, double const& b)
   {
      a[0] /= b;
      a[1] /= b;
   }

   /**
   * Stream extraction operator for fftw_complex
   *
   * \param is  input stream
   * \param z   complex number
   */
   std::istream& operator >> (std::istream& is, fftw_complex & z);

   /**
   * Stream insertion operator for fftw_complex
   *
   * \param os  output stream
   * \param z  complex number
   * \return modified output stream
   */
   std::ostream& operator << (std::ostream& os, fftw_complex const & z);

   #if 0
   /**
   * Serialization function template for fftw_complex number.
   *
   * Implementation serializes real part, then imaginary part.
   *
   * \param ar Archive object
   * \param z complex number
   * \param version version id
   */
   template <typename Archive>
   void serialize(Archive& ar, double z[2], const unsigned int version = 0)
   {
      serialize(ar, z[0], version);
      serialize(ar, z[1], version);
   }
   #endif

} // namespace Pscf
#endif

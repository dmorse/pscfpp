#ifndef PSCF_CPU_COMPLEX_H
#define PSCF_CPU_COMPLEX_H

/*
* PSCF Package - Polymer Self-Consistent Field 
*
* Copyright 2016 - 2023, The Regents of the University of Minnesota
* Distributed under the terms of the GNU General Public License.
*/

#include <fftw3.h>
#include <complex>

namespace Pscf {
namespace Cpu {

   // Real and imaginary components

   /**
   * Return the real part of a complex number.
   *
   * \ingroup Pscf_Cpu_Module
   *
   * \param a complex argument (input)
   */
   inline double real(fftw_complex const& a) 
   {  return a[0]; }

   /**
   * Return the imaginary part of a complex number.
   *
   * \ingroup Pscf_Cpu_Module
   *
   * \param a complex argument (input)
   */
   inline double imag(fftw_complex const& a) 
   {  return a[1]; }

   // Absolute magnitude

   /**
   * Return absolute magnitude of a complex number.
   *
   * \ingroup Pscf_Cpu_Module
   *
   * \param a complex argument (in)
   */
   inline double abs(fftw_complex const& a) 
   {  return sqrt(a[0] * a[0] + a[1] * a[1]); }

   /**
   * Return square of absolute magnitude of a complex number.
   *
   * \ingroup Pscf_Cpu_Module
   *
   * \param a complex argument (in)
   */
   inline double absSq(fftw_complex const& a) 
   {  return (a[0] * a[0] + a[1] * a[1]); }

   // Complex Conjugation

   /**
   * Compute complex conjugate, z = a^*.
   *
   * \ingroup Pscf_Cpu_Module
   *
   * \param z complex conjugate of argument (out)
   * \param a complex argument (in)
   */
   inline 
   void conj(fftw_complex& z, fftw_complex const& a)
   {
      z[0] = a[0];
      z[1] = -a[1];
   }

   /**
   * In-place complex conjugation of a complex number, a = a^* .
   *
   * \ingroup Pscf_Cpu_Module
   *
   * \param a argument (in) and complex conjugate (out)
   */
   inline 
   void conj(fftw_complex& a)
   {
      a[0] = a[0];
      a[1] = -a[1];
   }

   // Assignment 

   /**
   * Create a complex number from real and imaginary parts, z = a + ib.
   *
   * \ingroup Pscf_Cpu_Module
   *
   * \param z complex (out)
   * \param a real part (in)
   * \param b imaginary part (in)
   */
   inline 
   void assign(fftw_complex& z, double const& a, double const& b)
   {
      z[0] = a;
      z[1] = b;
   }

   /**
   * Assign a real input to a complex variable. 
   *
   * \ingroup Pscf_Cpu_Module
   *
   * \param z complex (out)
   * \param a real (in)
   */
   inline 
   void assign(fftw_complex& z, double const& a)
   {
      z[0] = a;
      z[1] = 0;
   }

   /**
   * Assign a complex input to a complex variable. 
   *
   * \ingroup Pscf_Cpu_Module
   *
   * \param z complex (out)
   * \param a complex (in)
   */
   inline 
   void assign(fftw_complex& z, fftw_complex const& a)
   {
      z[0] = a[0];
      z[1] = a[1];
   }

   /**
   * Assign a std::complex input to a complex variable.
   *
   * \ingroup Pscf_Cpu_Module
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
   * Assign a complex input to a std::complex variable.
   *
   * \ingroup Pscf_Cpu_Module
   *
   * \param z std::complex (out)
   * \param a complex (in)
   */
   inline 
   void assign(std::complex<double> & z, fftw_complex const& a)
   {  z = std::complex<double>(a[0], a[1]); }

   // Addition

   /**
   * Addition of two complex numbers, z = a + b.
   *
   * \ingroup Pscf_Cpu_Module
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
   * \ingroup Pscf_Cpu_Module
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
   * In place addition of complex numbers, a += b.
   *
   * \ingroup Pscf_Cpu_Module
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
   * In place addition of a complex and real number, a += b.
   *
   * \ingroup Pscf_Cpu_Module
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
   * Subtraction of two complex numbers, z = a - b.
   *
   * \ingroup Pscf_Cpu_Module
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
   * Subtraction of a real number from a complex number, z = a - b.
   *
   * \ingroup Pscf_Cpu_Module
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
   * In place subtraction of two complex numbers, a -= b.
   *
   * \ingroup Pscf_Cpu_Module
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
   * In place subtraction of real number from a complex number, a -= b.
   *
   * \ingroup Pscf_Cpu_Module
   *
   * \param a complex argument (in) and difference (out)
   * \param b real argument (in)
   */
   inline 
   void subEq(fftw_complex & a, double const& b)
   {   
      a[0] -= b; 
   }

   // Multiplication

   /**
   * Multiplication of two complex numbers, z = a * b.
   *
   * \ingroup Pscf_Cpu_Module
   *
   * \param z complex product (out)
   * \param a complex factor (in)
   * \param b complex factor (in)
   */
   inline 
   void mul(fftw_complex& z, fftw_complex const& a, fftw_complex const& b)
   {   
      z[0] = a[0]*b[0] - a[1]*b[1];
      z[1] = a[1]*b[0] + a[0]*b[1];
   }

   /**
   * Multiplication of complex and real numbers, z = a * b.
   *
   * \ingroup Pscf_Cpu_Module
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
   * In place multiplication of two complex number, a *= b.
   *
   * \ingroup Pscf_Cpu_Module
   *
   * \param a complex factor (in) and product (out)
   * \param b complex factor (in)
   */
   inline 
   void mulEq(fftw_complex & a, fftw_complex const& b)
   {  
      double a0;
      a0   = a[0]*b[0] - a[1]*b[1];
      a[1] = a[1]*b[0] + a[0]*b[1];
      a[0] = a0;
   }

   /**
   * In place multiplication of a complex and real number, a *= b.
   *
   * \ingroup Pscf_Cpu_Module
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

   // Division

   /**
   * Division of two complex numbers, z = a / b .
   *
   * \ingroup Pscf_Cpu_Module
   *
   * \param z complex ratio (out)
   * \param a complex numerator (in)
   * \param b complex denominator (in)
   */
   inline 
   void div(fftw_complex& z, fftw_complex const& a, fftw_complex const& b)
   {
      double bSq = b[0]*b[0] + b[1]*b[1];
      z[0] = (a[0]*b[0] + a[1]*b[1])/bSq;
      z[1] = (a[1]*b[0] - a[0]*b[1])/bSq;
   }

   /**
   * Division of a complex number by a real number, z = a / b .
   *
   * \ingroup Pscf_Cpu_Module
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
   * In place division of two complex number, a /= b.
   *
   * \ingroup Pscf_Cpu_Module
   *
   * \param a complex numerator (in) and ratio (out)
   * \param b complex denominator (in)
   */
   inline 
   void divEq(fftw_complex & a, fftw_complex const & b)
   {
      double bSq = b[0]*b[0] + b[1]*b[1];
      double a0 = (a[0]*b[0] + a[1]*b[1])/bSq;
      a[1] = (a[1]*b[0] - a[0]*b[1])/bSq;
      a[0] = a0;
   }

   /**
   * In place division of a complex number by a real number, a /= b.
   *
   * \ingroup Pscf_Cpu_Module
   *
   * \param a complex numerator (in) and ratio (out)
   * \param b real denominator (in)
   */
   inline 
   void divEq(fftw_complex & a, double const& b)
   {   
      a[0] /= b;
      a[1] /= b; 
   }

} // Pscf::Cpu
} // Pscf
#endif

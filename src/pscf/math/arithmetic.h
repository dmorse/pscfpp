#ifndef PSCF_MATH_ARITHMETIC_H
#define PSCF_MATH_ARITHMETIC_H

/*
* PSCF - Polymer Self-Consistent Field
*
* Copyright 2015 - 2026, The Regents of the University of Minnesota
* Distributed under the terms of the GNU General Public License.
*/

#include <cmath>
#include <complex>
//#include <iostream>

namespace Pscf {

   /**
   * \defgroup Pscf_Math_Arithmetic_Module Arithmetic Functions
   *
   * Overloaded functions and function templates for arithmetic operations 
   * using built-in real types and/or C++ standard library complex types.
   * Corresponding functions that use complex types used by the FFTW and
   * cufft libraries are defined in the files pscf/cpu/complex.h and 
   * pscf/cuda/complex.h
   *
   * Convention: Functions for which the result or output could be a 
   * complex number provide this as a modified value of the first 
   * parameter of the function, which must be passed as a non-const 
   * reference.  Functions for which the output is always real (such as 
   * an absolute value function) instead provide this as the function 
   * return value. This convention creates interfaces that allow the use
   * of data types to represent complex number for which no assignment 
   * operator is defined, which cannot be returned as function return 
   * values.
   *
   * \ingroup Pscf_Math_Module
   */

   /**
   * Assign a double value.
   *
   * \ingroup Pscf_Math_Arithmetic_Module
   *
   * \param z value (out)
   * \param a value (in)
   */
   inline
   void assign(double & z, double const & a)
   {  z = a; }

   /**
   * Assign a real float value.
   *
   * \ingroup Pscf_Math_Arithmetic_Module
   *
   * \param z value (out)
   * \param a value (in)
   */
   inline
   void assign(float & z, float const & a)
   {  z = a; }

   /**
   * Assign one std::complex<RT> variable to another.
   *
   * \ingroup Pscf_Math_Arithmetic_Module
   *
   * \param z value (out)
   * \param a value (in)
   */
   template <typename RT>  inline
   void assign(std::complex<RT> & z, std::complex<RT> const & a)
   {  z = a; }

   /**
   * Create std::complex<RT> from real and imaginary parts, z = a + ib.
   *
   * \ingroup Pscf_Math_Arithmetic_Module
   *
   * \param z complex (out)
   * \param a real part (in)
   * \param b imaginary part (in)
   */
   template <typename RT> inline
   void assign(std::complex<RT> & z, RT const & a, RT const & b)
   {  z = std::complex<RT>(a, b); }

   /**
   * Assign a real input to a std::complex variable.
   *
   * \ingroup Pscf_Math_Arithmetic_Module
   *
   * \param z complex (out)
   * \param a real (in)
   */
   template <typename RT> inline
   void assign(std::complex<RT> & z, RT const & a)
   {  z = a; }

   // Addition

   /**
   * Addition of two double precision numbers, z = a + b.
   *
   * \ingroup Pscf_Math_Arithmetic_Module
   *
   * \param z sum (out)
   * \param a summand (in)
   * \param b summand (in)
   */
   inline
   void add(double & z, double const & a, double const & b)
   {  z = a + b; }

   /**
   * Addition of two real float numbers, z = a + b.
   *
   * \ingroup Pscf_Math_Arithmetic_Module
   *
   * \param z sum (out)
   * \param a summand (in)
   * \param b summand (in)
   */
   inline
   void add(float & z, float const & a, float const & b)
   {  z = a + b; }

   /**
   * Addition of two std::complex variables, z = a + b.
   *
   * \ingroup Pscf_Math_Arithmetic_Module
   *
   * \param z sum (out)
   * \param a summand (in)
   * \param b summand (in)
   */
   template <typename RT> inline
   void add(std::complex<RT> & z, 
            std::complex<RT> const & a, std::complex<RT> const & b)
   {  z = a + b; }

   /**
   * Addition of a std::complex and real number, z = a + b.
   *
   * \ingroup Pscf_Math_Arithmetic_Module
   *
   * \param z complex sum (out)
   * \param a complex summand (in)
   * \param b real summand (in)
   */
   template <typename RT> inline
   void add(std::complex<RT> & z, 
            std::complex<RT> const & a, RT const & b)
   {  z = a + b; }

   /**
   * In place addition a += b (double).
   *
   * \ingroup Pscf_Math_Arithmetic_Module
   *
   * \param a summand (in) and sum (out)
   * \param b summand (in)
   */
   inline
   void addEq(double & a, double const & b)
   {  a += b; }

   /**
   * In place addition a += b (double).
   *
   * \ingroup Pscf_Math_Arithmetic_Module
   *
   * \param a summand (in) and sum (out)
   * \param b summand (in)
   */
   inline
   void addEq(float & a, float const & b)
   {  a += b; }

   /**
   * In place addition a += b (std::complex<RT>)
   *
   * \ingroup Pscf_Math_Arithmetic_Module
   *
   * \param a summand (in) and sum (out)
   * \param b summand (in)
   */
   template <typename RT> inline
   void addEq(std::complex<RT> & a, std::complex<RT> const & b)
   {  a += b; }

   /**
   * In place addition of std::complex and real numbers, a += b.
   *
   * \ingroup Pscf_Math_Arithmetic_Module
   *
   * \param a summand (in) and sum (out)
   * \param b summand (in)
   */
   template <typename RT> inline
   void addEq(std::complex<RT> & a, RT const & b)
   {  a += b; }

   // Subtraction

   /**
   * Subtraction of double precision numbers, z = a - b.
   *
   * \ingroup Pscf_Math_Arithmetic_Module
   *
   * \param z difference (out)
   * \param a 1st argument (in)
   * \param b 2nd argument (in)
   */
   inline
   void sub(double & z, double const & a, double const & b)
   {  z = a - b; }

   /**
   * Subtraction z = a - b (real float).
   *
   * \ingroup Pscf_Math_Arithmetic_Module
   *
   * \param z difference (out)
   * \param a 1st argument (in)
   * \param b 2nd argument (in)
   */
   inline
   void sub(float & z, float const & a, float const & b)
   {  z = a - b; }

   /**
   * Subtraction z = a - b (std::complex<RT>)
   *
   * \ingroup Pscf_Math_Arithmetic_Module
   *
   * \param z difference (out)
   * \param a 1st argument (in)
   * \param b 2nd argument (in)
   */
   template <typename RT> inline
   void sub(std::complex<RT> & z, 
            std::complex<RT> const & a, std::complex<RT> const & b)
   {  z = a - b; }

   /**
   * Subtraction of real from complex, z = a - b.
   *
   * \ingroup Pscf_Math_Arithmetic_Module
   *
   * \param z complex difference (out)
   * \param a complex 1st argument (in)
   * \param b real 2nd argument (in)
   */
   template <typename RT> inline
   void sub(std::complex<RT> & z, 
            std::complex<RT> const & a, RT const & b)
   {  z = a - b; }

   /**
   * In place subtraction a -= b (double).
   *
   * \ingroup Pscf_Math_Arithmetic_Module
   *
   * \param a argument (in) and difference (out)
   * \param b argument (in)
   */
   inline
   void subEq(double & a, double const & b)
   {  a -= b; }

   /**
   * In place subtraction, a -= b (float).
   *
   * \ingroup Pscf_Math_Arithmetic_Module
   *
   * \param a argument (in) and difference (out)
   * \param b argument (in)
   */
   inline
   void subEq(float & a, float const & b)
   {  a -= b; }

   /**
   * In place subtraction a -= b (std::complex<RT>).
   *
   * \ingroup Pscf_Math_Arithmetic_Module
   *
   * \param a argument (in) and difference (out)
   * \param b argument (in)
   */
   template <typename RT> inline
   void subEq(std::complex<RT> & a, std::complex<RT> const & b)
   {  a -= b; }

   /**
   * In place subtraction of real from std::complex, a -= b.
   *
   * \ingroup Pscf_Math_Arithmetic_Module
   *
   * \param a argument (in) and difference (out)
   * \param b argument (in)
   */
   template <typename RT> inline
   void subEq(std::complex<RT> & a, RT const & b)
   {  a -= b; }

   // Multiplication

   /**
   * Multiplication of double real numbers, z = a * b.
   *
   * \ingroup Pscf_Math_Arithmetic_Module
   *
   * \param z product (out)
   * \param a factor (in)
   * \param b factor (in)
   */
   inline
   void mul(double & z, double const & a, double const & b)
   {  z = a * b; }

   /**
   * Multiplication z = a * b (float).
   *
   * \ingroup Pscf_Math_Arithmetic_Module
   *
   * \param z product (out)
   * \param a factor (in)
   * \param b factor (in)
   */
   inline
   void mul(float & z, float const & a, float const & b)
   {  z = a * b; }

   /**
   * Multiplication z = a * b (std::complex).
   *
   * \ingroup Pscf_Math_Arithmetic_Module
   *
   * \param z product (out)
   * \param a factor (in)
   * \param b factor (in)
   */
   template <typename RT> inline
   void mul(std::complex<RT> & z, 
            std::complex<RT> const & a, std::complex<RT> const & b)
   {  z = a * b; }

   /**
   * Multiplication of std::complex and real, z = a * b.
   *
   * \ingroup Pscf_Math_Arithmetic_Module
   *
   * \param z complex product (out)
   * \param a complex factor (in)
   * \param b real factor (in)
   */
   template <typename RT> inline
   void mul(std::complex<RT> & z, 
            std::complex<RT> const & a, 
            RT const & b)
   {  z = a * b; }

   /**
   * In place multiplication a *= b (double).
   *
   * \ingroup Pscf_Math_Arithmetic_Module
   *
   * \param a factor (in) and product (out)
   * \param b factor (in)
   */
   inline
   void mulEq(double & a, double const & b)
   {  a *= b; }

   /**
   * In place multiplication a *= b (float).
   *
   * \ingroup Pscf_Math_Arithmetic_Module
   *
   * \param a factor (in) and product (out)
   * \param b factor (in)
   */
   inline
   void mulEq(float & a, float const & b)
   {  a *= b; }

   /**
   * Compute complex square of a std::complex, z = a * a.
   *
   * \ingroup Pscf_Math_Arithmetic_Module
   *
   * \param z complex product (out)
   * \param a complex factor (in)
   */
   template <typename RT> inline
   void square(std::complex<RT> & z, std::complex<RT> const & a)
   {  z = a * a; }

   // Division

   /**
   * Division z = a / b (double).
   *
   * \ingroup Pscf_Math_Arithmetic_Module
   *
   * \param z ratio (out)
   * \param a numerator (in)
   * \param b denominator (in)
   */
   inline
   void div(double & z, double const & a, double const & b)
   {  z = a/b; }

   /**
   * Division z = a / b (float).
   *
   * \ingroup Pscf_Math_Arithmetic_Module
   *
   * \param z ratio (out)
   * \param a numerator (in)
   * \param b denominator (in)
   */
   inline
   void div(float & z, float const & a, float const & b)
   {  z = a / b; }

   /**
   * Division z = a / b (std::complex)
   *
   * \ingroup Pscf_Math_Arithmetic_Module
   *
   * \param z ratio (out)
   * \param a numerator (in)
   * \param b denominator (in)
   */
   template <typename RT> inline
   void div(std::complex<RT> & z, 
            std::complex<RT> const & a, std::complex<RT> const & b)
   {  z = a / b; }

   /**
   * Division of a std::complex number by real, z = a / b .
   *
   * \ingroup Pscf_Math_Arithmetic_Module
   *
   * \param z complex ratio (out)
   * \param a complex numerator (in)
   * \param b real denominator (in)
   */
   template <typename RT> inline
   void div(std::complex<RT> & z, 
            std::complex<RT> const & a, RT const & b)
   {  z = a / b; }

   /**
   * In place division a /= b (double).
   *
   * \ingroup Pscf_Math_Arithmetic_Module
   *
   * \param a complex numerator (in) and ratio (out)
   * \param b complex denominator (in)
   */
   inline
   void divEq(double & a, double const & b)
   {  a /= b; }

   /**
   * In place division a /= b (float).
   *
   * \ingroup Pscf_Math_Arithmetic_Module
   *
   * \param a complex numerator (in) and ratio (out)
   * \param b complex denominator (in)
   */
   inline
   void divEq(float & a, float const & b)
   {  a /= b; }

   /**
   * In place division a /= b (std::complex).
   *
   * \ingroup Pscf_Math_Arithmetic_Module
   *
   * \param a complex numerator (in) and ratio (out)
   * \param b complex denominator (in)
   */
   template <typename RT> inline
   void divEq(std::complex<RT> & a, std::complex<RT> const & b)
   {  a /= b; }

   /**
   * In place division of std::complex number by real, a /= b.
   *
   * \ingroup Pscf_Math_Arithmetic_Module
   *
   * \param a complex numerator (in) and ratio (out)
   * \param b real denominator (in)
   */
   template <typename RT> inline
   void divEq(std::complex<RT> & a, RT const & b)
   {  a /= b; }

   // Inversion

   /**
   * Inverse, z = 1/a (double).
   *
   * \ingroup Pscf_Math_Arithmetic_Module
   *
   * \param z inverse (out)
   * \param a argument (in)
   */
   inline
   void inverse(double & z, double const & a)
   {  z = 1.0/a; }

   /**
   * Inverse, z = 1 / a (float).
   *
   * \ingroup Pscf_Math_Arithmetic_Module
   *
   * \param z inverse (out)
   * \param a argument (in)
   */
   inline
   void inverse(float & z, float const & a)
   {  z = 1.0/a; }

   // Exponential function

   /**
   * Exponentiation, z = exp(a) (double).
   *
   * \ingroup Pscf_Math_Arithmetic_Module
   *
   * \param z exponent (out)
   * \param a argument (in)
   */
   inline
   void assignExp(double & z, double const & a)
   {  z = std::exp(a); }

   /**
   * Exponent, z = exp(a) (float).
   *
   * \ingroup Pscf_Math_Arithmetic_Module
   *
   * \param z exponent (out)
   * \param a argument (in)
   */
   inline
   void assignExp(float & z, float const & a)
   {  z = std::exp(a); }

   // Natural logarithm function

   /**
   * Logarithm, z = exp(a) (double).
   *
   * \ingroup Pscf_Math_Arithmetic_Module
   *
   * \param z logarithm (out)
   * \param a argument (in)
   */
   inline
   void assignLog(double & z, double const & a)
   {  z = std::log(a); }

   /**
   * Logarithm, z = exp(a) (float).
   *
   * \ingroup Pscf_Math_Arithmetic_Module
   *
   * \param z logarithm (out)
   * \param a argument (in)
   */
   inline
   void assignLog(float & z, float const & a)
   {  z = std::log(a); }

} // namespace Pscf
#endif

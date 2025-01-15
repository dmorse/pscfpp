#ifndef PSCF_MATH_COMPLEX_H
#define PSCF_MATH_COMPLEX_H

/*
* PSCF Package  - Polymer Self-Consistent Field Theory
*
* Copyright 2016 - 2022, The Regents of the University of Minnesota
* Distributed under the terms of the GNU General Public License.
*/

#include <complex>
#include <iostream>

namespace Pscf {

   // Real and imaginary components

   /**
   * Return the real part of a complex number.
   *
   * \ingroup Pscf_Cpu_Complex_Module
   *
   * \param a complex argument (input)
   */
   template <typename CT, typename RT> 
   RT real(CT const& a);

   /**
   * Return the imaginary part of a complex number.
   *
   * \ingroup Pscf_Math_Complex_Module
   *
   * \param a complex argument (input)
   */
   template <typename CT, typename RT> 
   RT imag(CT const& a);

   // Absolute magnitude

   /**
   * Return absolute magnitude of a complex number.
   *
   * \ingroup Pscf_Math_Complex_Module
   *
   * \param a complex argument (in)
   */
   template <typename CT, typename RT> 
   RT abs(CT const& a);

   /**
   * Return square of absolute magnitude of a complex number.
   *
   * \ingroup Pscf_Math_Complex_Module
   *
   * \param a complex argument (in)
   */
   template <typename CT, typename RT> 
   RT absSq(CT const& a);

   // CT Conjugation

   /**
   * Compute complex conjugate, z = a^*.
   *
   * \ingroup Pscf_Math_Complex_Module
   *
   * \param z complex conjugate of argument (out)
   * \param a complex argument (in)
   */
   template <typename CT>
   void conj(CT& z, CT const& a);

   /**
   * In-place complex conjugation of a complex number, a = a^* .
   *
   * \ingroup Pscf_Math_Complex_Module
   *
   * \param a argument (in) and complex conjugate (out)
   */
   template <typename CT>
   void conj(CT& a);

   // Assignment

   /**
   * Create a complex number from real and imaginary parts, z = a + ib.
   *
   * \ingroup Pscf_Math_Complex_Module
   *
   * \param z complex (out)
   * \param a real part (in)
   * \param b imaginary part (in)
   */
   template <typename CT, typename RT>
   void assign(CT& z, RT const& a, RT const& b);

   /**
   * Assign a real input to a complex variable.
   *
   * \ingroup Pscf_Math_Complex_Module
   *
   * \param z complex (out)
   * \param a real (in)
   */
   template <typename CT, typename RT>
   void assign(CT& z, RT const& a);

   /**
   * Assign a complex input to a complex variable, z=a.
   *
   * \ingroup Pscf_Math_Complex_Module
   *
   * \param z complex (out)
   * \param a complex (in)
   */
   template <typename CT>
   void assign(CT& z, CT const& a);

   /**
   * Assign a std::complex input to a complex variable, z=a.
   *
   * \ingroup Pscf_Math_Complex_Module
   *
   * \param z complex (out)
   * \param a std::complex (in)
   */
   template <typename CT, typename RT>
   void assign(CT & z, std::complex<RT> const& a);

   /**
   * Assign a complex input to a std::complex variable, z=a.
   *
   * \ingroup Pscf_Math_Complex_Module
   *
   * \param z std::complex (out)
   * \param a complex (in)
   */
   template <typename CT, typename RT>
   void assign(std::complex<RT> & z, CT const& a);

   // Addition

   /**
   * Addition of two complex numbers, z = a + b.
   *
   * \ingroup Pscf_Math_Complex_Module
   *
   * \param z complex sum (out)
   * \param a complex summand (in)
   * \param b complex summand (in)
   */
   template <typename CT>
   void add(CT& z, CT const& a, CT const& b);

   /**
   * Addition of a complex and real number, z = a + b.
   *
   * \ingroup Pscf_Math_Complex_Module
   *
   * \param z complex sum (out)
   * \param a complex summand (in)
   * \param b real summand (in)
   */
   template <typename CT, typename RT>
   void add(CT& z, CT const& a, RT const& b);

   /**
   * In place addition of complex numbers, a += b.
   *
   * \ingroup Pscf_Math_Complex_Module
   *
   * \param a complex summand (in) and sum (out)
   * \param b complex summand (in)
   */
   template <typename CT>
   void addEq(CT& a, CT const& b);

   /**
   * In place addition of a complex and real number, a += b.
   *
   * \ingroup Pscf_Math_Complex_Module
   *
   * \param a complex summand (in) and sum (out)
   * \param b real summand (in)
   */
   template <typename CT, typename RT>
   void addEq(CT& a, RT const& b);

   // Subtraction

   /**
   * Subtraction of two complex numbers, z = a - b.
   *
   * \ingroup Pscf_Math_Complex_Module
   *
   * \param z complex difference (out)
   * \param a complex 1st argument (in)
   * \param b complex 2nd argument (in)
   */
   template <typename CT>
   void sub(CT& z, CT const& a, CT const& b);

   /**
   * Subtraction of a real number from a complex number, z = a - b.
   *
   * \ingroup Pscf_Math_Complex_Module
   *
   * \param z complex difference (out)
   * \param a complex 1st argument (in)
   * \param b real 2nd argument (in)
   */
   template <typename CT, typename RT>
   void sub(CT& z, CT const& a, RT const& b);

   /**
   * In place subtraction of two complex numbers, a -= b.
   *
   * \ingroup Pscf_Math_Complex_Module
   *
   * \param a complex argument (in) and difference (out)
   * \param b complex argument (in)
   */
   template <typename CT>
   void subEq(CT & a, CT const& b);

   /**
   * In place subtraction of real number from a complex number, a -= b.
   *
   * \ingroup Pscf_Math_Complex_Module
   *
   * \param a complex argument (in) and difference (out)
   * \param b real argument (in)
   */
   template <typename CT, typename RT>
   void subEq(CT & a, RT const& b);

   /**
   * Return square of the absolute magnitude of a complex difference.
   *
   * This function returns |a-b|^2 for complex a and b.
   *
   * \ingroup Pscf_Math_Complex_Module
   *
   * \param a complex 1st argument (in)
   * \param b complex 2nd argument (in)
   */
   template <typename CT, typename RT> 
   RT absSqDiff(CT const& a, CT const& b);

   // Multiplication

   /**
   * Multiplication of two complex numbers, z = a * b.
   *
   * \ingroup Pscf_Math_Complex_Module
   *
   * \param z complex product (out)
   * \param a complex factor (in)
   * \param b complex factor (in)
   */
   template <typename CT>
   void mul(CT& z, CT const& a, CT const& b);

   /**
   * Multiplication of complex and real numbers, z = a * b.
   *
   * \ingroup Pscf_Math_Complex_Module
   *
   * \param z complex product (out)
   * \param a complex factor (in)
   * \param b real factor (in)
   */
   template <typename CT, typename RT>
   void mul(CT& z, CT const& a, RT const& b);

   /**
   * In place multiplication of two complex number, a *= b.
   *
   * \ingroup Pscf_Math_Complex_Module
   *
   * \param a complex factor (in) and product (out)
   * \param b complex factor (in)
   */
   template <typename CT>
   void mulEq(CT & a, CT const& b);

   /**
   * In place multiplication of a complex and real number, a *= b.
   *
   * \ingroup Pscf_Math_Complex_Module
   *
   * \param a complex factor (in) and product (out)
   * \param b real factor (in)
   */
   template <typename CT, typename RT>
   void mulEq(CT & a, RT const& b);

   /**
   * Compute complex square of a complex number, z = a * a.
   *
   * \ingroup Pscf_Math_Complex_Module
   *
   * \param z complex product (out)
   * \param a complex factor (in)
   */
   template <typename CT>
   void square(CT& z, CT const& a);

   // Division

   /**
   * Division of two complex numbers, z = a / b .
   *
   * \ingroup Pscf_Math_Complex_Module
   *
   * \param z complex ratio (out)
   * \param a complex numerator (in)
   * \param b complex denominator (in)
   */
   template <typename CT>
   void div(CT& z, CT const& a, CT const& b);

   /**
   * Division of a complex number by a real number, z = a / b .
   *
   * \ingroup Pscf_Math_Complex_Module
   *
   * \param z complex ratio (out)
   * \param a complex numerator (in)
   * \param b real denominator (in)
   */
   template <typename CT, typename RT>
   void div(CT& z, CT const& a, RT const& b);

   /**
   * In place division of two complex number, a /= b.
   *
   * \ingroup Pscf_Math_Complex_Module
   *
   * \param a complex numerator (in) and ratio (out)
   * \param b complex denominator (in)
   */
   template <typename CT>
   void divEq(CT & a, CT const & b);

   /**
   * In place division of a complex number by a real number, a /= b.
   *
   * \ingroup Pscf_Math_Complex_Module
   *
   * \param a complex numerator (in) and ratio (out)
   * \param b real denominator (in)
   */
   template <typename CT, typename RT>
   void divEq(CT & a, RT const& b);

   #if 0
   /**
   * Stream extraction operator for CT
   *
   * \param is  input stream
   * \param z   complex number
   */
   std::istream& operator >> (std::istream& is, CT & z);

   /**
   * Stream insertion operator for CT
   *
   * \param os  output stream
   * \param z  complex number
   * \return modified output stream
   */
   std::ostream& operator << (std::ostream& os, CT const & z);
   #endif

} // namespace Pscf
#endif

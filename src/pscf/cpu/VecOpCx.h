#ifndef PSCF_CPU_VEC_OP_CX_H
#define PSCF_CPU_VEC_OP_CX_H

/*
* PSCF - Polymer Self-Consistent Field
*
* Copyright 2015 - 2026, The Regents of the University of Minnesota
* Distributed under the terms of the GNU General Public License.
*/

#include "VecOp.h"
#include <fftw3.h>

// Forward declaration
namespace Util {
   template <typename T> class Array;
}

namespace Pscf {

   using namespace Util;

   /*
   * Functions for element-wise complex vector operations on a CPU.
   *
   * Naming conventions are the same as those described in the file 
   * cpu/VecOp.h, which contains declarations of analogous functions
   * for real-valued arrays. Many operations involving complex arrays 
   * are implemented in both purely complex form, in which all inputs
   * are complex numbers, and "mixed" form in which one or more of the 
   * inputs is a real array or a real scalar.
   *
   * Functions defined in this file use the complex type fftw_complex
   * that is defined for use with the FFTW fast Fourier transform
   * library.
   */

   namespace VecOp {

      // Real and imaginary parts

      /**
      * Copy real part of a complex array to a real array.
      *
      * \ingroup Pscf_Cpu_VecOp_Module
      *
      * \param a  real array (LHS)
      * \param b  complex array (RHS)
      */
      void real(Array<double>& a, Array<fftw_complex> const & b);

      /**
      * Copy imaginary part of a complex array to a real array.
      *
      * \ingroup Pscf_Cpu_VecOp_Module
      *
      * \param a  real array (LHS)
      * \param b  complex array (RHS)
      */
      void imag(Array<double>& a, Array<fftw_complex> const & b);

      // Assignment

      /**
      * Vector assignment, a[i] = b[i] (complex).
      *
      * \ingroup Pscf_Cpu_VecOp_Module
      *
      * \param a  complex array (LHS)
      * \param b  complex array (RHS)
      */
      void eqV(Array<fftw_complex>& a, Array<fftw_complex> const & b);

      /**
      * Vector assignment, a[i] = (b[i], c[i]) (complex, real & imaginary).
      *
      * \ingroup Pscf_Cpu_VecOp_Module
      *
      * \param a  complex array (LHS)
      * \param b  real array, real part (RHS)
      * \param c  real array, imaginary part (RHS)
      */
      void eqV(Array<fftw_complex>& a,
               Array<double> const & b,
               Array<double> const & c);

      /**
      * Vector assignment, a[i] = b[i] (mixed, real b).
      *
      * \ingroup Pscf_Cpu_VecOp_Module
      *
      * \param a  complex array (LHS)
      * \param b  real array (RHS)
      */
      void eqV(Array<fftw_complex>& a, Array<double> const & b);

      /**
      * Vector-scalar assignment, a[i] = b (complex).
      *
      * \ingroup Pscf_Cpu_VecOp_Module
      *
      * \param a  complex array (LHS)
      * \param b  real scalar (RHS)
      */
      void eqS(Array<fftw_complex>& a, fftw_complex b);

      /**
      * Vector-scalar assignment, a[i] = b (mixed, real scalar b).
      *
      * \ingroup Pscf_Cpu_VecOp_Module
      *
      * \param a  complex array (LHS)
      * \param b  real scalar (RHS)
      */
      void eqS(Array<fftw_complex>& a, double b);

      // Addition

      /**
      * Vector addition, a[i] = b[i] + c[i] (complex).
      *
      * \ingroup Pscf_Cpu_VecOp_Module
      *
      * \param a  complex array (LHS)
      * \param b  complex array (RHS)
      * \param c  complex array (RHS)
      */
      void addVV(Array<fftw_complex>& a,
                 Array<fftw_complex> const & b,
                 Array<fftw_complex> const & c);

      /**
      * Vector addition, a[i] = b[i] + c[i] (mixed, real array c).
      *
      * \ingroup Pscf_Cpu_VecOp_Module
      *
      * \param a  complex array (LHS)
      * \param b  complex array (RHS)
      * \param c  real array (RHS)
      */
      void addVV(Array<fftw_complex> & a,
                 Array<fftw_complex> const & b,
                 Array<double> const & c);

      /**
      * Vector-scalar addition, a[i] = b[i] + c (complex).
      *
      * \ingroup Pscf_Cpu_VecOp_Module
      *
      * \param a  complex array (LHS)
      * \param b  complex array (RHS)
      * \param c  complex scalar (RHS)
      */
      void addVS(Array<fftw_complex> & a,
                 Array<fftw_complex> const & b,
                 fftw_complex c);

      /**
      * Vector-scalar addition, a[i] = b[i] + c (mixed).
      *
      * \ingroup Pscf_Cpu_VecOp_Module
      *
      * \param a  complex array (LHS)
      * \param b  complex array (RHS)
      * \param c  real scalar (RHS)
      */
      void addVS(Array<fftw_complex> & a,
                 Array<fftw_complex> const & b,
                 double c);

      // Subtraction

      /**
      * Vector subtraction, a[i] = b[i] - c[i] (complex).
      *
      * \ingroup Pscf_Cpu_VecOp_Module
      *
      * \param a  complex array (LHS)
      * \param b  complex array (RHS)
      * \param c  complex array (RHS)
      */
      void subVV(Array<fftw_complex> & a,
                 Array<fftw_complex> const & b,
                 Array<fftw_complex> const & c);

      /**
      * Vector subtraction, a[i] = b[i] - c[i] (mixed, real array c).
      *
      * \ingroup Pscf_Cpu_VecOp_Module
      *
      * \param a  complex array (LHS)
      * \param b  complex array (RHS)
      * \param c  real array (RHS)
      */
      void subVV(Array<fftw_complex> & a,
                 Array<fftw_complex> const & b,
                 Array<double> const & c);

      /**
      * Vector-scalar subtraction, a[i] = b[i] - c (complex).
      *
      * \ingroup Pscf_Cpu_VecOp_Module
      *
      * \param a  complex array (LHS)
      * \param b  complex array (RHS)
      * \param c  complex scalar (RHS)
      */
      void subVS(Array<fftw_complex> & a,
                 Array<fftw_complex> const & b,
                 fftw_complex c);

      /**
      * Vector-scalar subtraction, a[i] = b[i] - c (mixed, real scalar c).
      *
      * \ingroup Pscf_Cpu_VecOp_Module
      *
      * \param a  complex array (LHS)
      * \param b  complex array (RHS)
      * \param c  real scalar (RHS)
      */
      void subVS(Array<fftw_complex> & a,
                 Array<fftw_complex> const & b,
                 double c);

      // Multiplication

      /**
      * Vector multiplication, a[i] = b[i] * c[i] (complex).
      *
      * \ingroup Pscf_Cpu_VecOp_Module
      *
      * \param a  complex array (LHS)
      * \param b  complex array (RHS)
      * \param c  complex array (RHS)
      */
      void mulVV(Array<fftw_complex> & a,
                 Array<fftw_complex> const & b,
                 Array<fftw_complex> const & c);

      /**
      * Vector multiplication, a[i] = b[i] * c[i] (mixed, real array c).
      *
      * \ingroup Pscf_Cpu_VecOp_Module
      *
      * \param a  complex array (LHS)
      * \param b  complex array (RHS)
      * \param c  real array (RHS)
      */
      void mulVV(Array<fftw_complex> & a,
                 Array<fftw_complex> const & b,
                 Array<double> const & c);

      /**
      * Vector-scalar multiplication, a[i] = b[i] * c (complex).
      *
      * \ingroup Pscf_Cpu_VecOp_Module
      *
      * \param a  complex array (LHS)
      * \param b  complex array (RHS)
      * \param c  complex scalar (RHS)
      */
      void mulVS(Array<fftw_complex> & a,
                 Array<fftw_complex> const & b,
                 fftw_complex c);

      /**
      * Vector-scalar multiplication, a[i] = b[i] * c (mixed, real c).
      *
      * \ingroup Pscf_Cpu_VecOp_Module
      *
      * \param a  complex array (LHS)
      * \param b  complex array (RHS)
      * \param c  real scalar (RHS)
      */
      void mulVS(Array<fftw_complex> & a,
                 Array<fftw_complex> const & b,
                 double c);

      // Division

      /**
      * Vector division, a[i] = b[i] / c[i] (complex).
      *
      * \ingroup Pscf_Cpu_VecOp_Module
      *
      * \param a  complex array (LHS)
      * \param b  complex array (RHS)
      * \param c  complex array (RHS)
      */
      void divVV(Array<fftw_complex> & a,
                 Array<fftw_complex> const & b,
                 Array<fftw_complex> const & c);

      /**
      * Vector division, a[i] = b[i] / c[i] (mixed, real array c).
      *
      * \ingroup Pscf_Cpu_VecOp_Module
      *
      * \param a  complex array (LHS)
      * \param b  complex array (RHS)
      * \param c  real array (RHS)
      */
      void divVV(Array<fftw_complex> & a,
                 Array<fftw_complex> const & b,
                 Array<double> const & c);

      /**
      * Vector-scalar division, a[i] = b[i] / c (complex).
      *
      * \ingroup Pscf_Cpu_VecOp_Module
      *
      * \param a  complex array (LHS)
      * \param b  complex array (RHS)
      * \param c  complex scalar (RHS)
      */
      void divVS(Array<fftw_complex> & a,
                 Array<fftw_complex> const & b,
                 fftw_complex c);

      /**
      * Vector-scalar division, a[i] = b[i] / c (mixed, real scalar c).
      *
      * \ingroup Pscf_Cpu_VecOp_Module
      *
      * \param a  complex array (LHS)
      * \param b  complex array (RHS)
      * \param c  real scalar (RHS)
      */
      void divVS(Array<fftw_complex> & a,
                       Array<fftw_complex> const & b,
                       double c);

      // In-place addition

      /*
      * Vector in-place addition, a[i] += b[i] (complex).
      *
      * \ingroup Pscf_Cpu_VecOp_Module
      *
      * \param a  complex array (LHS)
      * \param b  complex array (RHS)
      */
      void addEqV(Array<fftw_complex> & a, Array<fftw_complex> const & b);

      /*
      * Vector in-place addition, a[i] += (b[i], c[i]) (complex).
      *
      * This function add real array b to the real part of complex array a,
      * and adds real array c to the imaginary part of complex array a.
      *
      * \ingroup Pscf_Cpu_VecOp_Module
      *
      * \param a  complex array (LHS)
      * \param b  real array, increment of real part (RHS)
      * \param c  real array, increment of imaginary part (RHS)
      */
      void addEqV(Array<fftw_complex> & a, 
                  Array<double> const & b,
                  Array<double> const & c);

      /*
      * Vector in-place addition, a[i] += b[i] (mixed, reall array b).
      *
      * \ingroup Pscf_Cpu_VecOp_Module
      *
      * \param a  complex array (LHS)
      * \param b  real array (RHS)
      */
      void addEqV(Array<fftw_complex> & a, Array<double> const & b);

      /*
      * Vector-scalar in-place addition, a[i] += b (complex).
      *
      * \ingroup Pscf_Cpu_VecOp_Module
      *
      * \param a  complex array (LHS)
      * \param b  complex scalar (RHS)
      */
      void addEqS(Array<fftw_complex> & a, fftw_complex b);

      /**
      * Vector-scalar in-place addition, a[i] += b (mixed, real scalar b).
      *
      * \ingroup Pscf_Cpu_VecOp_Module
      *
      * \param a  complex array (LHS)
      * \param b  real scalar (RHS)
      */
      void addEqS(Array<fftw_complex> & a, double b);

      // In-place subtraction

      /**
      * Vector in-place subtraction, a[i] -= b[i] (complex).
      *
      * \ingroup Pscf_Cpu_VecOp_Module
      *
      * \param a  complex array (LHS)
      * \param b  complex array (RHS)
      */
      void subEqV(Array<fftw_complex>& a, Array<fftw_complex> const & b);

      /**
      * Vector in-place subtraction, a[i] -= b[i] (mixed, real array b).
      *
      * \ingroup Pscf_Cpu_VecOp_Module
      *
      * \param a  complex array (LHS)
      * \param b  real array (RHS)
      */
      void subEqV(Array<fftw_complex>& a, Array<double> const & b);

      /**
      * Vector-scalar in-place subtraction, a[i] -= b (complex).
      *
      * \ingroup Pscf_Cpu_VecOp_Module
      *
      * \param a  complex array (LHS)
      * \param b  complex scalar (RHS)
      */
      void subEqS(Array<fftw_complex>& a, fftw_complex b);

      /**
      * Vector-scalar in-place subtraction, a[i] -= b (mixed, real b).
      *
      * \ingroup Pscf_Cpu_VecOp_Module
      *
      * \param a  complex array (LHS)
      * \param b  real scalar (RHS)
      */
      void subEqS(Array<fftw_complex>& a, double b);

      // In-place multiplication

      /**
      * Vector in-place multiplication, a[i] *= b[i] (complex).
      *
      * \ingroup Pscf_Cpu_VecOp_Module
      *
      * \param a  complex array (LHS)
      * \param b  complex array (RHS)
      */
      void mulEqV(Array<fftw_complex>& a, Array<fftw_complex> const & b);

      /**
      * Vector in-place multiplication, a[i] *= b[i] (mixed, real array b).
      *
      * \ingroup Pscf_Cpu_VecOp_Module
      *
      * \param a  complex array (LHS)
      * \param b  real array (RHS)
      */
      void mulEqV(Array<fftw_complex>& a, Array<double> const & b);

      /**
      * Vector-scalar in-place multiplication, a[i] *= b[i] (complex).
      *
      * \ingroup Pscf_Cpu_VecOp_Module
      *
      * \param a  complex array (LHS)
      * \param b  complex scalar (RHS)
      */
      void mulEqS(Array<fftw_complex>& a, fftw_complex const & b);

      /**
      * Vector-scalar in-place multiplication, a[i] *= b (mixed real b).
      *
      * \ingroup Pscf_Cpu_VecOp_Module
      *
      * \param a  complex array (LHS)
      * \param b  real scalar (RHS)
      */
      void mulEqS(Array<fftw_complex>& a, double b);

      // In-place division

      /**
      * Vector in-place division, a[i] /= b[i] (complex).
      *
      * \ingroup Pscf_Cpu_VecOp_Module
      *
      * \param a  complex array (LHS)
      * \param b  complex array (RHS)
      */
      void divEqV(Array<fftw_complex>& a, Array<fftw_complex> const & b);

      /**
      * Vector in-place division, a[i] /= b[i] (mixed, real array b).
      *
      * \ingroup Pscf_Cpu_VecOp_Module
      *
      * \param a  complex array (LHS)
      * \param b  real array (RHS)
      */
      void divEqV(Array<fftw_complex>& a, Array<double> const & b);

      /**
      * Vector-scalar in-place division, a[i] /= b (complex).
      *
      * \ingroup Pscf_Cpu_VecOp_Module
      *
      * \param a  complex array (LHS)
      * \param b  complex scalar (RHS)
      */
      void divEqS(Array<fftw_complex>& a, fftw_complex b);

      /**
      * Vector-scalar in-place division, a[i] /= b (mixed, real scalar b).
      *
      * \ingroup Pscf_Cpu_VecOp_Module
      *
      * \param a  complex array (LHS)
      * \param b  real scalar (RHS)
      */
      void divEqS(Array<fftw_complex>& a, double b);

      // Exponentiation

      /**
      * Vector exponentiation, a[i] = exp(b[i]) (complex).
      *
      * \ingroup Pscf_Cpu_VecOp_Module
      *
      * \param a  complex array (LHS)
      * \param b  complex array (RHS)
      */
      void expV(Array<fftw_complex> & a, Array<fftw_complex> const & b);

      // Complex Square

      /**
      * Elementwise complex square, a[i] = b[i] * b[i] (complex).
      *
      * \ingroup Pscf_Cpu_VecOp_Module
      *
      * \param a  complex array (LHS)
      * \param b  complex array (RHS)
      */
      void sqV(Array<fftw_complex> & a, Array<fftw_complex> const & b);

      // Even powers of complex absolute magnitude 

      /**
      * Square of absolute magnitude, a[i] = |b[i]|^2 (complex).
      *
      * \ingroup Pscf_Cpu_VecOp_Module
      *
      * \param a  real array (LHS)
      * \param b  complex array (RHS)
      */
      void sqAbsV(Array<double> & a, Array<fftw_complex> const & b);

      /**
      * Fourth power of absolute magnitude, a[i] = |b[i]|^4 (complex).
      *
      * \ingroup Pscf_Cpu_VecOp_Module
      *
      * \param a  real array (LHS)
      * \param b  complex array (RHS)
      */
      void sqSqAbsV(Array<double> & a, Array<fftw_complex> const & b);

      // Combined Operations

      /**
      * Vector division in-place w/ coeff., a[i] /= (b[i] * c).
      *
      * \ingroup Pscf_Cpu_VecOp_Module
      *
      * \param a  complex array (LHS)
      * \param b  real array (RHS)
      * \param c  input scalar (RHS)
      */
      void divEqVc(Array<fftw_complex>& a,
                   Array<double> const & b,
                   double const c);

   } // namespace VecOp

} // namespace Pscf
#endif

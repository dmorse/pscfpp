#ifndef PSCF_CPU_COMPLEX_H
#define PSCF_CPU_COMPLEX_H

/*
* PSCF Package - Polymer Self-Consistent Field 
*
* Copyright 2016 - 2023, The Regents of the University of Minnesota
* Distributed under the terms of the GNU General Public License.
*/

#include <fftw3.h>

namespace Pscf {
namespace Cpu {

   // Addition

   inline 
   void add(fftw_complex& z, fftw_complex const& a, fftw_complex const& b)
   {   
      z[0] = a[0] + b[0]; 
      z[1] = a[1] + b[1]; 
   }

   inline 
   void add(fftw_complex& z, fftw_complex const& a, double const& b)
   {   
      z[0] = a[0] + b; 
      z[1] = a[1]; 
   }

   inline
   void addEq(fftw_complex& a, fftw_complex const& b)
   {   
      a[0] += b[0]; 
      a[1] += b[1]; 
   }

   inline 
   void addEq(fftw_complex& a, double const& b)
   {   
      a[0] += b; 
   }

   // Subtraction

   inline 
   void sub(fftw_complex& z, fftw_complex const& a, fftw_complex const& b)
   {   
      z[0] = a[0] - b[0]; 
      z[1] = a[1] - b[1]; 
   }

   inline
   void sub(fftw_complex& z, fftw_complex const& a, double const& b)
   {   
      z[0] = a[0] - b; 
      z[1] = a[1]; 
   }

   inline 
   void subEq(fftw_complex & a, fftw_complex const& b)
   {   
      a[0] -= b[0]; 
      a[1] -= b[1]; 
   }

   inline 
   void subEq(fftw_complex & a, double const& b)
   {   
      a[0] -= b; 
   }

   // Multiplication

   inline 
   void mul(fftw_complex& z, fftw_complex const& a, fftw_complex const& b)
   {   
      z[0] = a[0]*b[0] - a[1]*b[1];
      z[1] = a[1]*b[0] + a[0]*b[1];
   }

   inline 
   void mul(fftw_complex& z, fftw_complex const& a, double const& b)
   {   
      z[0] = a[0]*b;
      z[1] = a[1]*b; 
   }

   inline 
   void mulEq(fftw_complex & a, fftw_complex const& b)
   {  
      double a0;
      a0   = a[0]*b[0] - a[1]*b[1];
      a[1] = a[1]*b[0] + a[0]*b[1];
      a[0] = a0;
   }

   inline 
   void mulEq(fftw_complex & a, double const& b)
   {   
      a[0] *= b;
      a[1] *= b; 
   }

   // Division

   inline 
   void div(fftw_complex& z, fftw_complex const& a, fftw_complex const& b)
   {
      double bSq = b[0]*b[0] + b[1]*b[1];
      z[0] = (a[0]*b[0] + a[1]*b[1])/bSq;
      z[1] = (a[1]*b[0] - a[0]*b[1])/bSq;
   }

   inline 
   void div(fftw_complex& z, fftw_complex const& a, double const& b)
   {   
      z[0] = a[0]/b;
      z[1] = a[1]/b; 
   }

   inline 
   void divEq(fftw_complex & a, fftw_complex const & b)
   {
      double bSq = b[0]*b[0] + b[1]*b[1];
      double a0 = (a[0]*b[0] + a[1]*b[1])/bSq;
      a[1] = (a[1]*b[0] - a[0]*b[1])/bSq;
      a[0] = a0;
   }

   inline 
   void divEq(fftw_complex & a, double const& b)
   {   
      a[0] /= b;
      a[1] /= b; 
   }

} // Pscf::Cpu
} // Pscf
#endif

#ifndef PSCF_CUDA_COMPLEX_H
#define PSCF_CUDA_COMPLEX_H

/*
* PSCF Package - Polymer Self-Consistent Field 
*
* Copyright 2016 - 2023, The Regents of the University of Minnesota
* Distributed under the terms of the GNU General Public License.
*/

#include <pscf/cuda/GpuTypes.h>

namespace Pscf {
namespace Cuda {

   // Real and imaginary components

   inline double real(cudaComplex const& a) 
   {  return a.x; }

   inline double imag(cudaComplex const& a) 
   {  return a.y; }

   // Absolute magnitude

   inline double abs(cudaComplex const& a) 
   {  return sqrt(a.x*a.x + a.y*a.y); }

   inline double absSq(cudaComplex const& a) 
   {  return (a.x*a.x + a.y*a.y); }

   // Complex Conjugation

   /*
   * Compute complex conjugate, a = b^*
   *
   * \param a result (left hand side)
   * \param b right hand side
   */
   inline 
   void conj(cudaComplex& a, cudaComplex const& b)
   {
      a.x = b.x;
      a.y = -b.y;
   }

   /*
   * In place conjugation of argument, a => a^*.
   *
   */
   inline 
   void conj(cudaComplex& a)
   {
      a.x = a.x;
      a.y = -a.y;
   }

   // Assignment 

   inline 
   void assign(cudaComplex& a, cudaComplex const& b)
   {
      a.x = b.x;
      a.y = b.y;
   }

   inline 
   void assign(cudaComplex & a, std::complex<double> const& b) 
   {  
      a.x = b.real();
      a.y = b.imag();
   }

   inline 
   void assign(std::complex<double> & a, cudaComplex const& b)
   {  a = std::complex<double>(b.x, b.y); }

   // Addition

   inline 
   void add(cudaComplex& z, cudaComplex const& a, cudaComplex const& b)
   {   
      z.x = a.x + b.x; 
      z.y = a.y + b.y; 
   }

   inline 
   void add(cudaComplex& z, cudaComplex const& a, cudaReal const& b)
   {   
      z.x = a.x + b; 
      z.y = a.y; 
   }

   inline
   void addEq(cudaComplex& a, cudaComplex const& b)
   {   
      a.x += b.x; 
      a.y += b.y; 
   }

   inline 
   void addEq(cudaComplex& a, cudaReal const& b)
   {   
      a.x += b; 
   }

   // Subtraction

   inline 
   void sub(cudaComplex& z, cudaComplex const& a, cudaComplex const& b)
   {   
      z.x = a.x - b.x; 
      z.y = a.y - b.y; 
   }

   inline
   void sub(cudaComplex& z, cudaComplex const& a, cudaReal const& b)
   {   
      z.x = a.x - b; 
      z.y = a.y; 
   }

   inline 
   void subEq(cudaComplex & a, cudaComplex const& b)
   {   
      a.x -= b.x; 
      a.y -= b.y; 
   }

   inline 
   void subEq(cudaComplex & a, cudaReal const& b)
   {   
      a.x -= b; 
   }

   // Multiplication

   inline 
   void mul(cudaComplex& z, cudaComplex const& a, cudaComplex const& b)
   {
      z.x = a.x * b.x - a.y * b.y;
      z.y = a.y * b.x + a.x * b.y;
   }

   inline 
   void mul(cudaComplex& z, cudaComplex const& a, cudaReal const& b)
   {   
      z.x = a.x * b;
      z.y = a.y * b; 
   }

   inline 
   void mulEq(cudaComplex & a, cudaComplex const& b)
   {  
      cudaReal a0;
      a0   = a.x * b.x - a.y * b.y;
      a.y = a.y * b.x + a.x * b.y;
      a.x = a0;
   }

   inline 
   void mulEq(cudaComplex & a, cudaReal const& b)
   {   
      a.x *= b;
      a.y *= b; 
   }

   // Division

   inline 
   void div(cudaComplex& z, cudaComplex const& a, cudaComplex const& b)
   {
      cudaReal bSq = b.x * b.x + b.y * b.y;
      z.x = (a.x * b.x + a.y * b.y)/bSq;
      z.y = (a.y * b.x - a.x * b.y)/bSq;
   }

   inline 
   void div(cudaComplex& z, cudaComplex const& a, cudaReal const& b)
   {   
      z.x = a.x/b;
      z.y = a.y/b; 
   }

   inline 
   void divEq(cudaComplex & a, cudaComplex const & b)
   {
      cudaReal bSq = b.x * b.x + b.y * b.y;
      cudaReal a0 = (a.x * b.x + a.y * b.y)/bSq;
      a.y = (a.y * b.x - a.x * b.y)/bSq;
      a.x = a0;
   }

   inline 
   void divEq(cudaComplex & a, cudaReal const& b)
   {   
      a.x /= b;
      a.y /= b; 
   }

} // Pscf::Cuda
} // Pscf
#endif

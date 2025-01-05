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

   // Addition

   inline 
   void add(hostComplex& z, hostComplex const& a, hostComplex const& b)
   {   
      z.x = a.x + b.x; 
      z.y = a.y + b.y; 
   }

   inline 
   void add(hostComplex& z, hostComplex const& a, hostReal const& b)
   {   
      z.x = a.x + b; 
      z.y = a.y; 
   }

   inline
   void addEq(hostComplex& a, hostComplex const& b)
   {   
      a.x += b.x; 
      a.y += b.y; 
   }

   inline 
   void addEq(hostComplex& a, hostReal const& b)
   {   
      a.x += b; 
   }

   // Subtraction

   inline 
   void sub(hostComplex& z, hostComplex const& a, hostComplex const& b)
   {   
      z.x = a.x - b.x; 
      z.y = a.y - b.y; 
   }

   inline
   void sub(hostComplex& z, hostComplex const& a, hostReal const& b)
   {   
      z.x = a.x - b; 
      z.y = a.y; 
   }

   inline 
   void subEq(hostComplex & a, hostComplex const& b)
   {   
      a.x -= b.x; 
      a.y -= b.y; 
   }

   inline 
   void subEq(hostComplex & a, hostReal const& b)
   {   
      a.x -= b; 
   }

   // Multiplication

   inline 
   void mul(hostComplex& z, hostComplex const& a, hostComplex const& b)
   {
      z.x = a.x * b.x - a.y * b.y;
      z.y = a.y * b.x + a.x * b.y;
   }

   inline 
   void mul(hostComplex& z, hostComplex const& a, hostReal const& b)
   {   
      z.x = a.x * b;
      z.y = a.y * b; 
   }

   inline 
   void mulEq(hostComplex & a, hostComplex const& b)
   {  
      hostReal a0;
      a0   = a.x * b.x - a.y * b.y;
      a.y = a.y * b.x + a.x * b.y;
      a.x = a0;
   }

   inline 
   void mulEq(hostComplex & a, hostReal const& b)
   {   
      a.x *= b;
      a.y *= b; 
   }

   // Division

   inline 
   void div(hostComplex& z, hostComplex const& a, hostComplex const& b)
   {
      hostReal bSq = b.x * b.x + b.y * b.y;
      z.x = (a.x * b.x + a.y * b.y)/bSq;
      z.y = (a.y * b.x - a.x * b.y)/bSq;
   }

   inline 
   void div(hostComplex& z, hostComplex const& a, hostReal const& b)
   {   
      z.x = a.x/b;
      z.y = a.y/b; 
   }

   inline 
   void divEq(hostComplex & a, hostComplex const & b)
   {
      hostReal bSq = b.x * b.x + b.y * b.y;
      hostReal a0 = (a.x * b.x + a.y * b.y)/bSq;
      a.y = (a.y * b.x - a.x * b.y)/bSq;
      a.x = a0;
   }

   inline 
   void divEq(hostComplex & a, hostReal const& b)
   {   
      a.x /= b;
      a.y /= b; 
   }

} // Pscf::Cuda
} // Pscf
#endif

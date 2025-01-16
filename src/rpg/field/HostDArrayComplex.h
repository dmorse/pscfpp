#ifndef RPG_HOST_D_ARRAY_COMPLEX_H
#define RPG_HOST_D_ARRAY_COMPLEX_H

/*
* PSCF - Polymer Self-Consistent Field Theory
*
* Copyright 2016 - 2022, The Regents of the University of Minnesota
* Distributed under the terms of the GNU General Public License.
*/

#include <pscf/cuda/HostDArray.h>     // base class
#include <pscf/cuda/GpuTypes.h>       // template parameter

namespace Pscf {
namespace Rpg {

   using namespace Util;

   class HostDArrayComplex : public HostDArray<cudaComplex>
   {

   public:

      typedef cudaComplex Complex;

      typedef cudaReal    Real;

      typedef HostDArray<cudaComplex> Base;

      /**
      * Default constructor.
      */
      HostDArrayComplex();

      /**
      * Allocating constructor.
      */
      HostDArrayComplex(int capacity);

      /**
      * Copy constructor.
      */
      HostDArrayComplex(HostDArrayComplex const & other);

      /**
      * Destructor.
      */
      ~HostDArrayComplex();

      using Base::allocate;
      using Base::deallocate;
      using Base::operator =;
      using Base::operator [];
      using Base::capacity;
      using Base::cArray;
      using Base::isAllocated;

   };

} // namespace Rpg
} // namespace Pscf
#endif

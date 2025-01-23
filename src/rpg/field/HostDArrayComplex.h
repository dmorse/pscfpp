#ifndef RPG_HOST_D_ARRAY_COMPLEX_H
#define RPG_HOST_D_ARRAY_COMPLEX_H

/*
* PSCF - Polymer Self-Consistent Field Theory
*
* Copyright 2016 - 2022, The Regents of the University of Minnesota
* Distributed under the terms of the GNU General Public License.
*/

#include <prdc/cuda/types.h>          // template parameter
#include <pscf/cuda/HostDArray.h>     // base class

namespace Pscf {
namespace Rpg {

   using namespace Util;

   /**
   * HostDArray containing cudaComplex elements.
   *
   * Defined to add typedefs.
   */
   class HostDArrayComplex : public HostDArray<cudaComplex>
   {

   public:

      /**
      * Type of each element.
      */
      typedef cudaComplex ElementType;

      /**
      * Complex number type.
      */
      typedef cudaComplex Complex;

      /**
      * Type of real or imaginary part of a Complex number.
      */
      typedef cudaReal    Real;

      /**
      * Base class type.
      */
      typedef HostDArray<cudaComplex> Base;

      // Member functions

      /**
      * Default constructor.
      */
      HostDArrayComplex();

      /**
      * Allocating constructor.
      *
      * This constructor allocates memory for the array.
      *
      * \param capacity  desired capacity of array
      */
      HostDArrayComplex(int capacity);

      /**
      * Copy constructor.
      *
      * Perform a deep copy of all array elements.
      *
      * \param other  other object being copied to this one.
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

#ifndef RPG_HOST_D_ARRAY_COMPLEX_H
#define RPG_HOST_D_ARRAY_COMPLEX_H

/*
* PSCF - Polymer Self-Consistent Field
*
* Copyright 2015 - 2025, The Regents of the University of Minnesota
* Distributed under the terms of the GNU General Public License.
*/

#include <pscf/cuda/cudaTypes.h>     // template parameter
#include <pscf/cuda/HostDArray.h>    // base class

namespace Pscf {
namespace Prdc {
namespace Cuda {

   using namespace Util;

   /**
   * HostDArray containing cudaComplex elements.
   *
   * Defined to add typedefs.
   */
   class HostDArrayComplex : public HostDArray<cudaComplex>
   {

   public:

      // Public type name aliases

      /**
      * Type of each element.
      */
      using ValueType = cudaComplex;

      /**
      * Type of real or imaginary part of a complex number.
      */
      using RealType = cudaReal;

      /**
      * Base class type.
      */
      using Base = HostDArray<cudaComplex>;

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

      // Inherited member functions
      using Base::allocate;
      using Base::deallocate;
      using Base::operator =;
      using Base::operator [];
      using Base::capacity;
      using Base::cArray;
      using Base::isAllocated;

   };

} // namespace Cuda
} // namespace Prdc
} // namespace Pscf
#endif

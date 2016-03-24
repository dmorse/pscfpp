#ifndef PSSP_FIELD_H
#define PSSP_FIELD_H

/*
* PSCF - Polymer Self-Consistent Field Theory
*
* Copyright 2016, The Regents of the University of Minnesota
* Distributed under the terms of the GNU General Public License.
*/

#include <util/containers/DArray.h>

namespace Pscf { 
namespace Pssp
{ 

   using namespace Util;

   template <typename T = double>
   class Field : protected DArray<T>
   {
 
      public:

      Field();
 
      // Access as one-dimensional array 
      T DArray<T>::operator [](int i);
      const T& DArray<T>::operator [](int i) const;
   
      // Assignment operator
      Field<T>& operator = (Field<T>& other);
   
      // Assignment, assigns all elements to a common scalar.
      Field<T>& operator = (T& scalar);
   
      // Increment operator
      Field<T>& operator += (Field<T>& other);
   
      // Decrement operator
      Field<T>& operator -= (Field<T>& other);
   
      // Multipication by scalar operator
      Field<T>& operator *= (T scalar);
   
      /**
      * Pointwise multipication of field operator.
      */
      GridField<D, T>& operator *= (GridField<D, T>& other);
   
      /**
      * Set to zero.
      */
      void setToZero();
   
      /**
      * Compute and return average of all elements.
      */
      T average() const;
   
   };
   
}
}
#endif

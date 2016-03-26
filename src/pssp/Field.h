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
   class Field : public DArray<T>
   {
   public:

      /**
      * Constructor.
      */
      Field();
 
      using Array<T>::operator [];
      using Array<T>::capacity;
  
      /** 
      * Assignment operator.
      */
      Field<T>& operator = (Field<T>& other);
   
      /** 
      * Assignment - assign all elements to a common scalar.
      */
      Field<T>& operator = (T& scalar);
   
      /** 
      * Increment operator - add one field by another.
      */
      Field<T>& operator += (Field<T>& other);
   
      /** 
      * Decrement operator - subtract one field from another.
      */
      Field<T>& operator -= (Field<T>& other);
   
      /** 
      * Multiplication operator - multiply one field by a scalar.
      */
      Field<T>& operator *= (T scalar);
   
      /**
      * Pointwise multipication of field operator.
      */
      Field<T>& operator *= (Field<T>& other);
   
      /**
      * Set to zero.
      */
      void setToZero();
   
      /**
      * Compute and return average of all elements.
      */
      T average() const;
  
   protected:
 
      using Array<T>::capacity_;
      using Array<T>::data_;

   };
   
   /*
   * Assignment& Field<T>::operator.
   */
   template <typename T>
   Field<T>& Field<T>::operator = (Field<T>& other)
   {
      for (int i = 0; i < capacity_; ++i) {
         data_[i] = other[i];
      }
      return *this;
   }

   /*
   * Assignment - assign all elements to a common scalar.
   */
   template <typename T>
   Field<T>& Field<T>::operator = (T& scalar)
   {
      for (int i = 0; i < capacity_; ++i) {
         data_[i] = scalar;
      }
      return *this;
   }

   /*
   * Increment& Field<T>::operator - add one field by another.
   */
   template <typename T>
   Field<T>& Field<T>::operator += (Field<T>& other)
   {
      for (int i = 0; i < capacity_; ++i) {
         data_[i] += other[i];
      }
      return *this;
   }

   /*
   * Decrement& Field<T>::operator - subtract one field from another.
   */
   template <typename T>
   Field<T>& Field<T>::operator -= (Field<T>& other)
   {
      for (int i = 0; i < capacity_; ++i) {
         data_[i] -= other[i];
      }
      return *this;
   }

   /*
   * Multiplication& Field<T>::operator - multiply one field by a scalar.
   */
   template <typename T>
   Field<T>& Field<T>::operator *= (T scalar)
   {
      for (int i = 0; i < capacity_; ++i) {
         data_[i] *= scalar;
      }
      return *this;
   }

   /*
   * Pointwise multipication of field& Field<T>::operator.
   */
   template <typename T>
   Field<T>& Field<T>::operator *= (Field<T>& other)
   {
      for (int i = 0; i < capacity_; ++i) {
         data_[i] *= other[i];
      }
      return *this;
   }

   /*
   * Set to zero.
   */
   template <typename T>
   void Field<T>::setToZero()
   {
      for (int i = 0; i < capacity_; ++i) {
         data_[i] *= 0.0;
      }
   }

   /*
   * Compute and return average of all elements.
   */
   template <typename T>
   T Field<T>::average() const
   {
      double value = 0.0;
      for (int i = 0; i < capacity_; ++i) {
         value += data_[i];
      }
      return value/T(capacity_);
   }

}
}
#endif

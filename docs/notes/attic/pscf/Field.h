#ifndef PSCF_FIELD_H
#define PSCF_FIELD_H

/*
* PSCF - Polymer Self-Consistent Field Theory
*
* Copyright 2016 - 2022, The Regents of the University of Minnesota
* Distributed under the terms of the GNU General Public License.
*/

#include <util/containers/DArray.h>

namespace Pscf { 

   using namespace Util;

   /**
   * Base class template for a field defined on a spatial grid.
   *
   * Derived from DArray<T>, and provides useful arithmetic operations.
   *
   * \ingroup Pscf_Math_Module
   */
   template <typename T = double>
   class Field : public DArray<T>
   {
   public:

      /**
      * Constructor.
      */
      Field();
 
      /** 
      * Copy constructor.
      */
      Field(Field<T> const & other);
  
      /** 
      * Assignment operator.
      */
      Field<T>& operator = (Field<T> const & other);
   
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
      * Pointwise multipication of one field by another.
      */
      Field<T>& operator *= (Field<T>& other);
   
      /**
      * Set all elements to zero.
      */
      void setToZero();
   
      /**
      * Compute and return average of all elements.
      */
      T average() const;
  
      using Array<T>::operator [];
      using Array<T>::capacity;
      using DArray<T>::allocate;
      using DArray<T>::deallocate;
      using DArray<T>::isAllocated;

   protected:
 
      using Array<T>::capacity_;
      using Array<T>::data_;

   };
   
   /*
   * Copy constructor.
   *
   * Allocates new memory and copies all elements by value.
   *
   *\param other the Field to be copied.
   */
   template <class T>
   Field<T>::Field(Field<T> const & other)
    : DArray<T>()
   {
      if (!other.isAllocated()) {
         UTIL_THROW("Other Field not allocated.");
      }
      Memory::allocate(data_, other.capacity_);
      capacity_ = other.capacity_;
      for (int i = 0; i < capacity_; ++i) {
         data_[i] = other.data_[i];
      }
   }

   /*
   * Assignment from another Field, element-by-element.
   *
   * This operator will allocate memory if not allocated previously.
   *
   * \throw Exception if other Field is not allocated.
   * \throw Exception if both Fields are allocated with unequal capacities.
   *
   * \param other the rhs Field
   */
   template <class T>
   Field<T>& Field<T>::operator = (Field<T> const & other)
   {
      // Check for self assignment
      if (this == &other) return *this;

      // Precondition
      if (!other.isAllocated()) {
         UTIL_THROW("Other Field not allocated.");
      }

      if (!isAllocated()) {
         allocate(other.capacity());
      } else if (capacity_ != other.capacity_) {
         UTIL_THROW("Fields of unequal capacity");
      }

      // Copy elements
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
      if (!isAllocated()) {
         UTIL_THROW("Field not allocated.");
      }
      for (int i = 0; i < capacity_; ++i) {
         data_[i] = scalar;
      }
      return *this;
   }

   /*
   * Increment& Field<T>::operator, add one field by another.
   */
   template <typename T>
   Field<T>& Field<T>::operator += (Field<T>& other)
   {
      if (!other.isAllocated()) {
         UTIL_THROW("Other Field no allocated.");
      }
      if (!isAllocated()) {
         UTIL_THROW("This Field not allocated.");
      }
      if (capacity_ != other.capacity_) {
         UTIL_THROW("Fields of unequal capacity");
      }
      for (int i = 0; i < capacity_; ++i) {
         data_[i] += other.data_[i];
      }
      return *this;
   }

   /*
   * Decrement& Field<T>::operator, subtract one field from another.
   */
   template <typename T>
   Field<T>& Field<T>::operator -= (Field<T>& other)
   {

      // Preconditions
      if (!other.isAllocated()) {
         UTIL_THROW("Other Field not allocated.");
      }
      if (!isAllocated()) {
         UTIL_THROW("This Field not allocated.");
      }
      if (capacity_ != other.capacity_) {
         UTIL_THROW("Fields of unequal capacity");
      }

      for (int i = 0; i < capacity_; ++i) {
         data_[i] -= other.data_[i];
      }
      return *this;
   }

   /*
   * Multiplication& Field<T>::operator - multiply one field by a scalar.
   */
   template <typename T>
   Field<T>& Field<T>::operator *= (T scalar)
   {
      // Precondition
      if (!isAllocated()) {
         UTIL_THROW("Field not allocated.");
      }

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
      // Preconditions
      if (!other.isAllocated()) {
         UTIL_THROW("Other Field not allocated.");
      }
      if (!isAllocated()) {
         UTIL_THROW("This Field not allocated.");
      }
      if (capacity_ != other.capacity_) {
         UTIL_THROW("Unequal capacity");
      }

      for (int i = 0; i < capacity_; ++i) {
         data_[i] *= other.data_[i];
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
#endif

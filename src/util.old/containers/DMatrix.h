#ifndef UTIL_D_MATRIX_H
#define UTIL_D_MATRIX_H

/*
* Simpatico - Simulation Package for Polymeric and Molecular Liquids
*
* Copyright 2010 - 2014, The Regents of the University of Minnesota
* Distributed under the terms of the GNU General Public License.
*/

#include <util/containers/Matrix.h>
#include <util/misc/Memory.h>
#include <util/global.h>

namespace Util
{

   /**
   * Dynamically allocated Matrix.
   *
   * \ingroup Matrix_Module
   */
   template <typename Data>
   class DMatrix : public Matrix<Data>
   {

      using Matrix<Data>::data_;
      using Matrix<Data>::capacity1_;
      using Matrix<Data>::capacity2_;

   public:

      /**
      * Constructor.
      */
      DMatrix();

      /**
      * Copy constructor.
      */
      DMatrix(const DMatrix<Data>& other);

      /**
      * Assignment.
      *
      * \throw Exception if LHS and RHS dimensions do not match.
      */
      DMatrix<Data>& operator= (const DMatrix<Data>& other);

      /**
      * Destructor.
      *
      * Delete dynamically allocated C array.
      */
      ~DMatrix();

      /**
      * Allocate memory for a matrix.
      *
      * \param capacity1 number of rows (range of first index)
      * \param capacity2 number of columns (range of second index)
      */
      void allocate(int capacity1, int capacity2);

      /**
      * Serialize a DMatrix to/from an Archive.
      *
      * \param ar       archive 
      * \param version  archive version id
      */
      template <class Archive>
      void serialize(Archive& ar, const unsigned int version);

      /**
      * Return true if the DMatrix has been allocated, false otherwise.
      */
      bool isAllocated() const;

   };

   // Method definitions

   /*
   * Default constructor.
   */
   template <typename Data>
   DMatrix<Data>::DMatrix() :
      Matrix<Data>()
   {}
      
   /*
   * Copy constructor.
   */
   template <typename Data>
   DMatrix<Data>::DMatrix(const DMatrix<Data>& other) 
     : Matrix<Data>()
   {
      // Precondition
      if (other.data_ == 0) {
         UTIL_THROW("Other DMatrix must be allocated");
      }

      // Allocate and copy elements
      allocate(other.capacity1_, other.capacity2_);
      for (int i = 0; i < capacity1_*capacity2_; ++i) {
         data_[i] = other.data_[i];
      }
   }
   
   /*
   * Assignment.
   */
   template <typename Data>
   DMatrix<Data>& DMatrix<Data>::operator = (const DMatrix<Data>& other) 
   {
      // Check for self assignment.
      if (this == &other) return *this;

      // Precondition
      if (other.data_ == 0) {
         UTIL_THROW("RHS DMatrix must be allocated in assignment");
      }

      if (data_ == 0) {
         // If this DMatrix if not allocated, allocate now.
         allocate(other.capacity1_, other.capacity2_);
      } else {
         // If this is allocated, check that capacities are equal.
         if (capacity1_ != other.capacity1_ || 
             capacity2_ != other.capacity2_) {
            UTIL_THROW("Unequal capacities in assignment");
         }
      }

      // Copy elements
      for (int i = 0; i < capacity1_*capacity2_; ++i) {
         data_[i] = other.data_[i];
      }

      return *this;
   }

   /*
   * Destructor.
   *
   * Delete dynamically allocated C array.
   */
   template <typename Data>
   DMatrix<Data>::~DMatrix()
   {
      if (data_) {
         Memory::deallocate<Data>(data_, capacity1_*capacity2_);
      }
   }

   /*
   * Allocate memory for a matrix.
   *
   * \param capacity1 number of rows (range of first index)
   * \param capacity2 number of columns (range of second index)
   */
   template <typename Data>
   void DMatrix<Data>::allocate(int capacity1, int capacity2)
   {
      // Preconditions
      if (capacity1 <= 0) UTIL_THROW("Capacity1 must be positive");
      if (capacity2 <= 0) UTIL_THROW("Capacity2 must be positive");
      if (data_  != 0) UTIL_THROW("Attempt to re-allocate a Matrix");

      Memory::allocate<Data>(data_, capacity1*capacity2);
      capacity1_ = capacity1;
      capacity2_ = capacity2;
   }

   /*
   * Serialize a DMatrix to/from an Archive.
   */
   template <class Data>
   template <class Archive>
   void DMatrix<Data>::serialize(Archive& ar, const unsigned int version)
   {
      int capacity1, capacity2;
      if (Archive::is_saving()) {
         if (!isAllocated()) {
            UTIL_THROW("Cannot save unallocated DMatrix.");
         }
         capacity1 = capacity1_;
         capacity2 = capacity2_;
      }
      ar & capacity1;
      ar & capacity2;
      if (Archive::is_loading()) {
         if (!isAllocated()) {
            if (capacity1 > 0 && capacity2 > 0) {
               allocate(capacity1, capacity2);
            } else {
               UTIL_THROW("Unallocated DMatrix, invalid capacity values");
            }
         } else {
            if (capacity1 != capacity1_) {
               UTIL_THROW("Inconsistent DMatrix capacity1 values");
            }
            if (capacity2 != capacity2_) {
               UTIL_THROW("Inconsistent DMatrix capacity2 values");
            }
         }
      }
      for (int i = 0; i < capacity1_*capacity2_; ++i) {
         ar & data_[i];
      }
   }

   /*
   * Return true if the DMatrix has been allocated, false otherwise.
   */
   template <class Data>
   inline bool DMatrix<Data>::isAllocated() const 
   {  return !(data_ == 0); }

}
#endif

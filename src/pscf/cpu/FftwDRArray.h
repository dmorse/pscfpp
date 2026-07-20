#ifndef PSCF_FFTW_DR_ARRAY_H
#define PSCF_FFTW_DR_ARRAY_H

/*
* Util Package - C++ Utilities for Scientific Computation
*
* Copyright 2010 - 2017, The Regents of the University of Minnesota
* Distributed under the terms of the GNU General Public License.
*/

#include <util/containers/Array.h>        // base class
#include <util/misc/ReferenceCounter.h>   // member
#include <util/misc/CountedReference.h>   // member
#include <util/misc/Memory.h>             // member
#include <util/global.h>

#include <fftw3.h>

using namespace Util;

namespace Pscf {

   /**
   * Dynamic reference counted array for use with FFTW library.
   *
   * The allocate and deallocate functions of this class use functions
   * provided by the FFTW library to allocate and free aligned memory.
   * The class is otherwise identical tot he Util::DRArray class.
   *
   * A FftwDRArray may be in any of three states:
   *
   *   (1) Unallocated: In this state, there is no associated memory,
   *   so capacity() returns 0, while isAllocated(), isOwner() and
   *   isAssociated() all return false.
   *
   *   (2) A data owner: In this case, this object owns a C array that
   *   it is responsible for de-allocating. In this state, capacity()
   *   returns a positive integer, isAllocated() and isOwner() return
   *   true, and isAssociated() returns false.
   *
   *   (3) A data user: In this case, this object has a pointer to a C
   *   array that is owned by a different FftwDRArray object. We describe
   *   this by saying that this FftwDRArray (the data user) is "associated"
   *   with a C array that is owned by another object (the data owner), or
   *   that the data user "references" that array. In this state, capacity()
   *   returns a positive value, isAllocated() and isAssociated() return
   *   true, and isOwner() returns false.
   *
   * A FftwDRArray that owns a C array that is referenced by one or more
   * other associated FftwDRArray objects maintains a count of how many
   * other such objects reference its data. This counter is automatically
   * incremented when a reference is created and decremented when an
   * existing reference is destroyed.
   *
   * When a FftwDArray is allocated (i.e., either a data owner or user)
   * array elements may be accessed via a subscript operator (an
   * overloaded operator []) that is inherited from the Array<Data> base
   * class.  Member functions for memory management allow a FftwDRArray
   * to allocate or deallocate a C array that it owns, or to create or
   * release an association with an array slice that it does not own.
   *
   * It is a logical error to invoke the deallocate() member function of
   * a FftwDRArray that is unallocated or that references data that it does
   * not own. In either case, an Exception is thrown.  It is also an error
   * to attempt to deallocate an FftwDRArray that is referenced by one
   * other associated FftwDRArray data users.
   *
   * \ingroup Prdc_Cpu_Module
   */
   template <typename Data>
   class FftwDRArray : public Array<Data>
   {

   public:

      /**
      * Data type of each array element.
      */
      using ValueType = Data;

      /**
      * Default constructor.
      */
      FftwDRArray();

      /**
      * Allocating constructor.
      *
      * This function calls allocate(capacity) internally.
      *
      * \param capacity number of elements to allocate
      */
      FftwDRArray(int capacity);

      /**
      * Copy constructor.
      *
      * \param other  the FftwDRArray to be copied
      */
      FftwDRArray(FftwDRArray<Data> const & other);

      /**
      * Destructor.
      *
      * Deletes any C array that is owned by this object, and releases any
      * association with a C Array that is referred to but not owned by
      * this object. If this object owns an array that is referred to by
      * one or more other FftwDRArray objects, an error message is written
      * to std::cout.
      */
      ~FftwDRArray();

      /**
      * Assignment from another FftwDRArray<Data> container.
      *
      * Performs a deep copy, by copying values of all elements of another
      * FftwDRArray<Data> container. If this LHS array is already allocated
      * on entry, it must have the same capacity as the other RHS array.
      * If this LHS array is not allocated on entry, required memory is
      * allocated before copying values. After exit, isAllocated() and
      * isOwner() return true, while isAssociated() returns false.
      *
      * \throw Exception if other array is not allocated
      * \throw Exception if arrays are allocated with unequal capacities
      *
      * \param other  array container on RHS of assigment (input)
      */
      FftwDRArray<Data>& operator = (FftwDRArray<Data> const & other);

      /**
      * Assignment from an Array<Data> container.
      *
      * Performs a deep copy, by copying values of all elements of an
      * Array<Data> container. If this (LHS) array is already allocated
      * on entry, it must have the same capacity as the other (RHS) array.
      * If this LHS array is not allocated on entry, required memory is
      * allocated before copying values. After exit, isAllocated() and
      * isOwner() return true, while isAssociated() returns false.
      *
      * \throw Exception if other array is not allocated
      * \throw Exception if arrays are allocated with unequal capacities
      *
      * \param other  array container on RHS of assigment (input)
      */
      FftwDRArray<Data>& operator = (Array<Data> const & other);

      /**
      * Allocate an underlying C array, which this container then owns.
      *
      * On entry, this object must be unallocated, i.e., it must not have
      * data that it either owns or references.  After exit, isAllocated()
      * and isOwner() will return true, while isAssociated() will return
      * false.
      *
      * \throw Exception if this array is allocated on entry.
      *
      * \param capacity number of elements to allocate
      */
      void allocate(int capacity);

      /**
      * Dellocate an underlying C array that is owned by this container.
      *
      * After exit, isAllocated(), isOwner(), and isAssociated() will
      * all return false.
      *
      * \throw Exception if this object does not own data.
      */
      void deallocate();

      /**
      * Associate this object with a slice of a different FftwDRArray.
      *
      * On entry, this object must be not be allocated, i.e., it must
      * not have data that it either owns or references.  After exit,
      * isAllocated() and isAssociated() will return true, while isOwner()
      * will return false.
      *
      * \throw Exception if this array is allocated on entry.
      *
      * \param arr  parent array that owns the data
      * \param beginId  index in the parent array at which this array starts
      * \param capacity  number of elements associated with this container
      */
      void associate(FftwDRArray<Data>& arr, int beginId, int capacity);

      /**
      * Dissociate this object from an externally owned array slice.
      *
      * After exit, isAllocated(), isOwner(), and isAssociated() will
      * all return false.
      *
      * \throw Exception if this is not associated with external data
      */
      void dissociate();

      /**
      * Serialize an FftwDRArray to/from an Archive.
      *
      * \param ar       archive
      * \param version  archive version id
      */
      template <class Archive>
      void serialize(Archive& ar, const unsigned int version);

      /**
      * Return true if this container has data, false otherwise.
      *
      * A FftwDRArray is considered allocated if it has non-null pointer
      * to a C array, which may either be an array that it owns or a
      * slice of an array that is owned by another FftwDRArray object.
      */
      bool isAllocated() const;

      /**
      * Does this container own a dynamically allocated C array?
      *
      * If isAllocated() is false, isOwner() is also false.
      * If isAllocated() is true, either isOwner() or isAsssociated()
      * must be true, but not both.
      */
      bool isOwner() const;

      /**
      * Is this container associated with a C array it does not own?
      *
      * If isAllocated() is false, isAssociated() is also false.
      */
      bool isAssociated() const;

   protected:

      using Array<Data>::data_;
      using Array<Data>::capacity_;

   private:

      /// Counter for any containers that reference data owned by this.
      ReferenceCounter refCounter_;

      /// Reference to a container that owns memory referenced by this.
      CountedReference ref_;

   };

   // Inline member function definitions

   /*
   * Does this FftwDRArray have data (either owned or associated) ?
   */
   template <typename Data> inline
   bool FftwDRArray<Data>::isAllocated() const
   {  return (bool)data_; }

   /*
   * Does this object own data?
   */
   template <typename Data> inline
   bool FftwDRArray<Data>::isOwner() const
   {  return ((bool) data_ && !ref_.isAssociated()); }

   /*
   * Does this object reference data that it does not own?
   */
   template <typename Data> inline
   bool FftwDRArray<Data>::isAssociated() const
   {  return ((bool) data_ && ref_.isAssociated()); }

   /*
   * Serialize a FftwDArray to/from an Archive.
   */
   template <typename Data>
   template <class Archive>
   void FftwDRArray<Data>::serialize(Archive& ar, 
                                     const unsigned int version)
   {
      int capacity;
      if (Archive::is_saving()) {
         capacity = capacity_;
         if (capacity > 0) {
            UTIL_CHECK(isOwner());
         }
      }
      ar & capacity;
      if (Archive::is_loading()) {
         if (!isAllocated()) {
            if (capacity > 0) {
               allocate(capacity);
            }
         } else {
            UTIL_CHECK(capacity == capacity_);
            UTIL_CHECK(isOwner());
         }
      }
      UTIL_CHECK(capacity == capacity_);
      if (capacity > 0) {
         UTIL_CHECK(isOwner());
         for (int i = 0; i < capacity_; ++i) {
            ar & data_[i];
         }
      }
   }

   #if 0
   // Non-inline member function definitions

   /*
   * Default constructor.
   */
   template <typename Data>
   FftwDRArray<Data>::FftwDRArray()
    : Array<Data>()
   {}

   /*
   * Allocating constructor.
   */
   template <typename Data>
   FftwDRArray<Data>::FftwDRArray(int capacity)
    : Array<Data>()
   {  allocate(capacity); }

   /*
   * Copy constructor.
   */
   template <typename Data>
   FftwDRArray<Data>::FftwDRArray(FftwDRArray<Data> const & other)
    : Array<Data>()
   {
      if (!other.isAllocated()) {
         UTIL_THROW("Other FftwDRArray must be allocated.");
      }
      allocate(other.capacity_);
      for (int i = 0; i < capacity_; ++i) {
         data_[i] = other.data_[i];
      }
   }

   /*
   * Destructor.
   */
   template <typename Data>
   FftwDRArray<Data>::~FftwDRArray()
   {
      if (data_) {
         if (ref_.isAssociated()) {
            ref_.dissociate();
         } else {
            if (refCounter_.hasRefs()) {
               int nRef = refCounter_.nRef();
               std::cout
                   << std::endl
                   << "Error: Destroying a FftwDRArray that is referenced by "
                   << nRef << " other(s)" << std::endl;
            }
            fftw_free(data_);
         }
      }
      data_ = nullptr;
      capacity_ = 0;
   }

   /*
   * Assignment from another FftwDRArray<Data> (deep copy).
   */
   template <typename Data>
   FftwDRArray<Data>&
   FftwDRArray<Data>::operator = (FftwDRArray<Data> const & other)
   {
      // Check for self assignment
      if (this == &other) return *this;

      // Precondition - other array must be allocated
      UTIL_CHECK (other.isAllocated());

      // If this array is not allocated, then allocate
      if (!isAllocated()) {
         allocate(other.capacity());
      }

      // Copy elements
      UTIL_CHECK(capacity_ == other.capacity_);
      for (int i = 0; i < capacity_; ++i) {
         data_[i] = other[i];
      }

      return *this;
   }

   /*
   * Assignment from an Array<Data> (deep copy).
   */
   template <typename Data>
   FftwDRArray<Data>&
   FftwDRArray<Data>::operator = (Array<Data> const & other)
   {
      // Check for self assignment
      if (dynamic_cast< Array<Data>* >(this) == &other) return *this;

      // Precondition - other array must be allocated
      UTIL_CHECK(other.capacity() > 0);

      // If this array is not allocated, then allocate
      if (!isAllocated()) {
         allocate(other.capacity());
      }

      // Copy all elements
      UTIL_CHECK(capacity_ == other.capacity());
      for (int i = 0; i < capacity_; ++i) {
         data_[i] = other[i];
      }

      return *this;
   }

   /*
   * Allocate an underlying C array, which this container then owns.
   */
   template <typename Data>
   void FftwDRArray<Data>::allocate(int capacity)
   {
      if (capacity <= 0) {
         UTIL_THROW("Attempt to allocate with capacity <= 0");
      }
      if (isAllocated()) {
         UTIL_THROW("Attempt to allocate a FftwDRArray that has data.");
      }
      data_ = (Data*) fftw_malloc(sizeof(Data) * capacity);
      capacity_ = capacity;
   }

   /*
   * Deallocate a C array that is owned by this container.
   */
   template <typename Data>
   void FftwDRArray<Data>::deallocate()
   {
      UTIL_CHECK(data_);
      UTIL_CHECK(!ref_.isAssociated());
      fftw_free(data_);
      data_ = nullptr;
      capacity_ = 0;
      UTIL_CHECK(!refCounter_.hasRefs());
   }

   /*
   * Associate this object with a slice of a different FftwDRArray.
   */
   template <typename Data>
   void FftwDRArray<Data>::associate(FftwDRArray<Data>& owner,
                                 int beginId, int capacity)
   {
      UTIL_CHECK(owner.isAllocated());
      UTIL_CHECK(owner.isOwner());
      UTIL_CHECK(beginId >= 0);
      UTIL_CHECK(capacity > 0);
      UTIL_CHECK(beginId + capacity <= owner.capacity());
      UTIL_CHECK(!data_);
      UTIL_CHECK(!ref_.isAssociated());

      // Copy data pointer and capacity
      data_ = owner.cArray() + beginId;
      capacity_ = capacity;

      // Associate private ReferencecCounter of the data owner with the
      // CountedReference ref_ member variable of this data user.
      ref_.associate(owner.refCounter_);

      // On exit from CountedReference::associate, the ReferenceCounter
      // of the data owner is incremented and the CountedReference of
      // this data user has a pointer to that ReferenceCounter.
   }

   /*
   * Dissociate this object from array slice owned by another object.
   */
   template <typename Data>
   void FftwDRArray<Data>::dissociate()
   {
      UTIL_CHECK(data_);
      UTIL_CHECK(ref_.isAssociated());

      data_ = nullptr;
      capacity_ = 0;
      ref_.dissociate(); // decrements counter mainained by owner
   }
   #endif

} // namespace Pscf
#include "FftwDRArray.tpp"
#endif

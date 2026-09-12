#ifndef PSCF_FFTW_DR_ARRAY_H
#define PSCF_FFTW_DR_ARRAY_H

/*
* Util Package - C++ Utilities for Scientific Computation
*
* Copyright 2010 - 2017, The Regents of the University of Minnesota
* Distributed under the terms of the GNU General Public License.
*/

#include <util/containers/ArraySource.h>  // base class
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
   * This class is otherwise identical to the Util::DRArray class.
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
   * \ingroup Pscf_Backend_Cpp_Module
   */
   template <typename Data>
   class FftwDRArray : public ArraySource<Data>
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
      * \throw Exception if other array is not a data owner on entry.
      *
      * \param other  parent array that owns the data
      * \param beginId  index in of other array at which slice begins
      * \param capacity  number of elements in the slice
      */
      void associate(FftwDRArray<Data>& other, int beginId, int capacity);

      /**
      * Associate this object with all of a different FftwDRArray.
      *
      * This function associates this FftwDRArray with all of another
      * FftwDRArray that is a data owner. The function associate(other)
      * is equivalent to associate(other, 0, other.capacity()).
      *
      * On entry, this object must be not be allocated, i.e., it must
      * not have data that it either owns or references, while the other
      * array must own data. After exit, isAllocated() and isAssociated()
      * will return true, while isOwner() will return false.
      *
      * \throw Exception if this array is allocated
      * \throw Exception if other array is not a data owner
      *
      * \param other  array that owns the data
      */
      void associate(FftwDRArray<Data>& other);

      /**
      * Dissociate this object from an externally owned array slice.
      *
      * After exit, isAllocated(), isOwner(), and isAssociated() will
      * all return false.
      *
      * \throw Exception if this is not associated with another array
      */
      void dissociate();

      /**
      * Serialize an FftwDRArray to/from an Archive.
      *
      * Precondition: To serialize to a write archive, this array must
      * be unallocated or a data owner.
      *
      * \throw Exception if a write is attempted for a data user
      *
      * \param ar       archive
      * \param version  archive version id
      */
      template <class Archive>
      void serialize(Archive& ar, const unsigned int version);

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

      /*
      * A FftwDRArray is considered allocated if it has non-null pointer
      * to a C array, which may either be an array that it owns or a
      * slice of an array that is owned by another FftwDRArray object.
      */
      using Array<Data>::isAllocated;

   protected:

      using Array<Data>::data_;
      using Array<Data>::capacity_;

   private:

      /// Reference to a container that owns memory referenced by this.
      CountedReference ref_;

      // Prohibit public access to the reference counter.
      using ReferenceCounter::hasRefs;
      using ReferenceCounter::nRef;

      // Note: ReferenceCounter::nRef_ is declared mutable, so changes
      // to this variable should not be visible in the public interface.

   };

   // Inline member function definitions

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

} // namespace Pscf
#include "FftwDRArray.tpp"
#endif

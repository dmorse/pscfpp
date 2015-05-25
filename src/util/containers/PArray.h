#ifndef UTIL_P_ARRAY_H
#define UTIL_P_ARRAY_H

/*
* Simpatico - Simulation Package for Polymeric and Molecular Liquids
*
* Copyright 2010 - 2014, The Regents of the University of Minnesota
* Distributed under the terms of the GNU General Public License.
*/

#include <util/containers/PArrayIterator.h>
#include <util/containers/ConstPArrayIterator.h>
#include <util/global.h>

namespace Util
{

   /**
   * An array that only holds pointers to its elements.
   *
   * A PArray<Data> is an array that is implemented by storing pointers 
   * to Data objects, rather than actual Data objects. The array suscript 
   * [] operator returns a reference to an associated Data object, as for 
   * Array<Data>. A PArray<Data> is not responsible for destroying the
   * associated Data objects.
   *
   * A PArray cannot be instantiated, because its constructor is protected.
   * PArray is a base class for DPArray and for ArraySet.
   *
   * \ingroup Pointer_Array_Module
   */
   template <typename Data>
   class PArray
   {

   public:

      // Protected default constructor

      /**
      * Destructor.
      */
      virtual ~PArray();

      /**
      * Return allocated size.
      *
      * \return Number of elements allocated in array.
      */
      int capacity() const;

      /**
      * Return logical size.
      *
      * \return logical size of this array.
      */
      int size() const;

      /**
      * Set a PArrayIterator to the beginning of this PArray.
      *
      * \param iterator PArrayIterator, initialized on output.
      */
      void begin(PArrayIterator<Data> &iterator) const;

      /**
      * Set a ConstPArrayIterator to the beginning of this PArray.
      *
      * \param iterator PArrayIterator, initialized on output.
      */
      void begin(ConstPArrayIterator<Data> &iterator) const;

      /**
      * Mimic C array subscripting.
      *
      * \param  i array index
      * \return reference to element i
      */
      Data& operator[] (int i) const;

   protected:

      /// Constructor (protected to prevent instantiation).
      PArray();

      /// PArray of of pointers to Data objects.
      Data** ptrs_;

      /// Allocated size of ptrs_ array.
      int capacity_;

      /// Logical size (number of elements with initialized data).
      int size_;

   private:

      /**
      * Copy constructor, private to prohibit copy construction.
      */
      PArray(const PArray<Data>& other);

      /**
      * Assignment, private to prohibit assignment.
      */
      PArray<Data>& operator = (const PArray<Data>& other);

   };

   /*
   * Constructor.
   */
   template <typename Data>
   inline PArray<Data>::PArray()
    : ptrs_(0),
      capacity_(0),
      size_(0)
   {}

   /*
   * Destructor.
   */
   template <typename Data>
   PArray<Data>::~PArray()
   {}

   /*
   * Return allocated size.
   */
   template <typename Data>
   inline int PArray<Data>::capacity() const
   {  return capacity_; }

   /*
   * Return logical size.
   */
   template <typename Data>
   inline int PArray<Data>::size() const
   {  return size_; }

   /**
   * Set an PArrayIterator to the beginning of this PArray.
   *
   * \param iterator PArrayIterator, initialized on output. 
   */
   template <typename Data>
   inline void PArray<Data>::begin(PArrayIterator<Data> &iterator) const
   {
      if (ptrs_ && size_ > 0) {
         iterator.setCurrent(ptrs_);
         iterator.setEnd(ptrs_ + size_);
      } else {
         iterator.setNull();
      }
   }

   /**
   * Set an ConstPArrayIterator to the beginning of this PArray.
   *
   * \param iterator ConstPArrayIterator, initialized on output. 
   */
   template <typename Data>
   inline void PArray<Data>::begin(ConstPArrayIterator<Data> &iterator) const
   {
      if (ptrs_ && size_ > 0) {
         iterator.setCurrent(ptrs_);
         iterator.setEnd(ptrs_ + size_);
      } else {
         iterator.setNull();
      }
   }

   /**
   * Subscript - return a reference.
   *
   * \param  i array index
   * \return reference to element i
   */
   template <typename Data>
   inline Data& PArray<Data>::operator[] (int i) const
   {
      assert(ptrs_);
      assert(i >= 0);
      assert(i < size_);
      return *(ptrs_[i]);
   }

} 
#endif

#ifndef PSSP_RFIELD_H
#define PSSP_RFIELD_H

/*
* PSCF++ Package 
*
* Copyright 2010 - 2017, The Regents of the University of Minnesota
* Distributed under the terms of the GNU General Public License.
*/

#include <util/global.h>

namespace Pssp
{

   using namespace Util;

   /**
   * Field of real double precision values for FFT transform.
   */
   class RField
   {

   public:

      /**
      * Default constructor.
      */
      RField();

      /**
      * Copy constructor.
      *
      * Allocates new memory and copies all elements by value.
      *
      *\param other the RField to be copied.
      */
      RField(const RField& other);

      /**
      * Destructor.
      *
      * Deletes underlying C array, if allocated previously.
      */
      virtual ~RField();

      /**
      * Assignment operator.
      *
      * If this RField is not allocated, allocates and copies all elements.
      *
      * If this and the other RField are both allocated, the capacities must
      * be exactly equal. If so, this method copies all elements.
      *
      * \param other the RHS RField
      */
      RField& operator = (const RField& other);

      /**
      * Allocate the underlying C array.
      *
      * \throw Exception if the RField is already allocated.
      *
      * \param capacity number of elements to allocate.
      */
      void allocate(int capacity);

      /**
      * Dellocate the underlying C array.
      *
      * \throw Exception if the RField is not allocated.
      */
      void deallocate();

      /**
      * Return true if the RField has been allocated, false otherwise.
      */
      bool isAllocated() const;

      /**
      * Return allocated size.
      *
      * \return Number of elements allocated in array.
      */
      int capacity() const;

      /**
      * Get an element by non-const reference.
      *
      * Mimic C-array subscripting.
      *
      * \param  i array index
      * \return non-const reference to element i
      */
      double& operator[] (int i);

      /**
      * Get an element by const reference.
      *
      * Mimics C-array subscripting.
      *
      * \param i array index
      * \return const reference to element i
      */
      const double& operator[] (int i) const;

      /**
      * Return pointer to underlying C array.
      */
      double* cRField();

      /**
      * Return pointer to const to underlying C array.
      */
      const double* cRField() const;

      /**
      * Serialize a RField to/from an Archive.
      *
      * \param ar       archive
      * \param version  archive version id
      */
      template <class Archive>
      void serialize(Archive& ar, const unsigned int version);

   protected:

      /// Pointer to an array of double elements.
      double* data_;

      /// Allocated size of the data_ array.
      int capacity_;

   };

   /*
   * Return allocated size.
   */
   inline int RField::capacity() const
   {  return capacity_; }

   /*
   * Get an element by reference (C-array subscripting)
   */
   inline double& RField::operator[] (int i)
   {
      assert(data_ != 0);
      assert(i >= 0);
      assert(i < capacity_);
      return *(data_ + i);
   }

   /*
   * Get an element by const reference (C-array subscripting)
   */
   inline const double& RField::operator[] (int i) const
   {
      assert(data_ != 0);
      assert(i >= 0 );
      assert(i < capacity_);
      return *(data_ + i);
   }

   /*
   * Get a pointer to the underlying C array.
   */
   inline double* RField::cRField()
   {  return data_; }

   /*
   * Get a pointer to const to the underlying C array.
   */
   inline const double* RField::cRField() const
   {  return data_; }

   /*
   * Return true if the RField has been allocated, false otherwise.
   */
   inline bool RField::isAllocated() const
   {  return (bool)data_; }

   /*
   * Serialize a RField to/from an Archive.
   */
   template <class Archive>
   void RField::serialize(Archive& ar, const unsigned int version)
   {
      int capacity;
      if (Archive::is_saving()) {
         capacity = capacity_;
      }
      ar & capacity;
      if (Archive::is_loading()) {
         if (!isAllocated()) {
            if (capacity > 0) {
               allocate(capacity);
            }
         } else {
            if (capacity != capacity_) {
               UTIL_THROW("Inconsistent RField capacities");
            }
         }
      }
      if (isAllocated()) {
         for (int i = 0; i < capacity_; ++i) {
            ar & data_[i];
         }
      }
   }

}
#endif

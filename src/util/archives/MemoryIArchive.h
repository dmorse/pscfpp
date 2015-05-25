#ifndef UTIL_MEMORY_I_ARCHIVE_H
#define UTIL_MEMORY_I_ARCHIVE_H

/*
* Simpatico - Simulation Package for Polymeric and Molecular Liquids
*
* Copyright 2010 - 2014, The Regents of the University of Minnesota
* Distributed under the terms of the GNU General Public License.
*/

#include "Byte.h"
#include "serialize.h"

#include <util/space/Vector.h>
#include <util/space/IntVector.h>

#include <complex>
#include <string>
#include <vector>

namespace Util
{

   class MemoryOArchive;

   /**
   * Input archive for packed heterogeneous binary data.
   *
   * \ingroup Serialize_Module
   */
   class MemoryIArchive 
   {

   public:

      /// Returns true;
      static bool is_saving();

      /// Returns false;
      static bool is_loading();

      /**
      * Constructor.
      */
      MemoryIArchive();

      /**
      * Destructor.
      */
      ~MemoryIArchive();

      /**
      * Allocate memory block.
      *
      * \param capacity sizeof of block, in Bytes.
      */
      void allocate(size_t capacity);

      /**
      * Assignment from MemoryOArchive.
      */
      MemoryIArchive& operator = (MemoryOArchive& other);

      /**
      * Reset the cursor to the beginning (for rereading).
      */
      void reset();

      /**
      * Reset to empty state.
      *
      * Resets cursor and end pointers to beginning of memory block. 
      */
      void clear();

      /**
      * Release memory obtained by assignment.
      */
      void release();

      /**
      * Load one object.
      *
      * \param data object to be loaded from this archive.
      */
      template <typename T>
      MemoryIArchive& operator & (T& data);

      /**
      * Load one object.
      *
      * \param data object to be loaded from this archive.
      */
      template <typename T>
      MemoryIArchive& operator >> (T& data);

      /**
      * Unpack one object of type T.
      *
      * \param data object to be loaded from this archive.
      */
      template <typename T>
      void unpack(T& data);
   
      /**
      * Read a C-array of objects of type T.
      *
      * \param array array into which data should be loaded.
      * \param n     expected number of elements in the array.
      */
      template <typename T>
      void unpack(T* array, int n);

      /**
      * Unpack a 2D C array.
      * 
      * Unpack m rows of n elements into array of type T array[mp][np],
      * with m <= mp and n <= np. 
      *
      * \param array pointer to [0][0] element of 2D array
      * \param m  logical number of rows
      * \param n  logical number of columns
      * \param np physical number of columns
      */
      template <typename T> 
      void unpack(T* array, int m, int n, int np);

      #ifdef UTIL_MPI
      /**
      * Receive packed data via MPI.
      *
      * \param comm   MPI communicator
      * \param source rank of processor from which data is sent.
      */
      void recv(MPI::Intracomm& comm, int source);
      #endif

      /**
      * Return pointer to beginning of block.
      */
      Byte* begin() const;

      /**
      * Return pointer to current position (cursor).
      */
      Byte* cursor() const;

      /**
      * Return pointer to end of packed block (one Byte past the last).
      */
      Byte* end() const;

      /**
      * Return capacity in Bytes.
      */
      size_t capacity() const;

      /**
      * Has memory been allocated?
      */
      bool isAllocated() const;

   private:

      /// Pointer to first byte of allocated memory in block. 
      Byte* buffer_;

      /// Pointer to first element in block. 
      Byte* begin_;

      /// Current element (read/write cursor).
      Byte* cursor_;

      /// End of block containing data (one byte past last).
      Byte* end_;

      /// Pointer to one Byte past last in allocated block.
      Byte* endAllocated_;

      /// Pointer to associated MemoryOArchive, if any.
      MemoryOArchive* oArchivePtr_;

      /// Allocated size of send and recv buffers, in Bytes.
      size_t capacity_;

      /// Archive version number.
      unsigned int version_;

      /// Did this archive allocate the memory block?
      bool ownsData_;

   };

   // Inline static methods

   inline bool MemoryIArchive::is_saving()
   { return false; }

   inline bool MemoryIArchive::is_loading()
   { return true; }

   // Inline non-static methods

   /*
   * Return pointer to beginning of block.
   */
   inline Byte* MemoryIArchive::begin() const
   {  return begin_; }

   /*
   * Return pointer to cursor position.
   */
   inline Byte* MemoryIArchive::cursor() const
   {  return cursor_; }

   /*
   * Return end of packed block (one Byte past the last).
   */
   inline Byte* MemoryIArchive::end() const
   {  return end_; }

   /*
   * Return capacity in Bytes.
   */
   inline size_t MemoryIArchive::capacity() const
   {  return capacity_; }

   /*
   * Has a memory block been allocated?
   */
   inline bool MemoryIArchive::isAllocated() const
   {  return (bool) begin_; }

   // Template methods

   /*
   * Load one T object from a MemoryIArchive.
   */
   template <typename T>
   inline MemoryIArchive& MemoryIArchive::operator & (T& data)
   {  
      serialize(*this, data, version_); 
      return *this;
   }

   /*
   * Load one T object from this MemoryIArchive.
   */
   template <typename T>
   inline MemoryIArchive& MemoryIArchive::operator >> (T& data)
   {   
      serialize(*this, data, version_); 
      return *this;
   }

   /*
   * Load one T object from this MemoryIArchive.
   */
   template <typename T>
   void MemoryIArchive::unpack(T& data)
   {
      if (cursor_ + sizeof(data) > end_) {
         UTIL_THROW("Attempted read past end of packed block");
      }
      T* ptr = (T *)cursor_;
      data = *ptr;
      ++ptr;
      cursor_ = (Byte *)ptr;
   }

   /*
   * Load a C-array of objects of type T.
   */
   template <typename T>
   void MemoryIArchive::unpack(T* array, int n)
   {
      if (cursor_ + n*sizeof(T) > end_) {
         UTIL_THROW("Attempted read past end of data");
      }
      T* ptr = (T *)cursor_;
      for (int i=0; i < n; ++i) {
         array[i] = *ptr;
         ++ptr;
      }
      cursor_ = (Byte *)ptr;
   }

   /*
   * Bitwise pack a 2D C-array of objects of type T.
   */
   template <typename T>
   void MemoryIArchive::unpack(T* array, int m, int n, int np)
   {
      if (cursor_ + m*n*sizeof(T) > end_) {
         UTIL_THROW("Attempted read past end of data");
      }
      int i, j;
      T* ptr = (T *)cursor_;
      for (i = 0; i < m; ++i) {
         for (j = 0; j < n; ++j) {
            array[i*np + j] = *ptr;
            ++ptr;
         }
      }
      cursor_ = (Byte *)ptr;
   }

   // Explicit specializations of serialize function

   /*
   * Load a char from a MemoryIArchive.
   */
   template <>
   inline void serialize(MemoryIArchive& ar, char& data, 
                         const unsigned int version)
   {  ar.unpack(data); }

   /*
   * Load a bool from a MemoryIArchive.
   */
   template <>
   inline void serialize(MemoryIArchive& ar, bool& data, 
                         const unsigned int version)
   {  ar.unpack(data); }

   /*
   * Load an unsigned int from a MemoryIArchive.
   */
   template <>
   inline void serialize(MemoryIArchive& ar, unsigned int& data, 
                         const unsigned int version)
   {  ar.unpack(data); }

   /*
   * Load an int from a MemoryIArchive.
   */
   template <>
   inline void serialize(MemoryIArchive& ar, int& data, 
                         const unsigned int version)
   {  ar.unpack(data); }

   /*
   * Load an unsigned long int from a MemoryIArchive.
   */
   template <>
   inline void serialize(MemoryIArchive& ar, unsigned long& data, 
                         const unsigned int version)
   {  ar.unpack(data); }

   /*
   * Load a long int from a MemoryIArchive.
   */
   template <>
   inline void serialize(MemoryIArchive& ar, long& data, 
                         const unsigned int version)
   {  ar.unpack(data); }

   /*
   * Load a float from a MemoryIArchive.
   */
   template <>
   inline void serialize(MemoryIArchive& ar, float& data, 
                         const unsigned int version)
   {  ar.unpack(data); }

   /*
   * Load a double from a MemoryIArchive.
   */
   template <>
   inline void serialize(MemoryIArchive& ar, double& data, 
                         const unsigned int version)
   {  ar.unpack(data); }

   /*
   * Load a std::vector from a MemoryIArchive.
   */
   template <typename T>
   void serialize(MemoryIArchive& ar, std::vector<T>& data, 
                  const unsigned int version)
   {
      T element;
      std::size_t size;
      ar.unpack(size);
      data.reserve(size);
      data.clear();
      for (size_t i = 0; i < size; ++i) {
         ar & element;
         data.push_back(element);
      }
   }

   // Explicit serialize methods for std library types.

   /*
   * Load a std::complex<float> from a MemoryIArchive.
   */
   template <>
   inline void serialize(MemoryIArchive& ar, std::complex<float>& data,        
                         const unsigned int version)
   {  ar.unpack(data); }

   /*
   * Load a std::complex<double> from a MemoryIArchive.
   */
   template <>
   inline void serialize(MemoryIArchive& ar, std::complex<double>& data, 
                         const unsigned int version)
   {  ar.unpack(data); }

   /*
   * Load a std::string from MemoryIArchive.
   */
   template <>
   void serialize(MemoryIArchive& ar, std::string& data, 
                         const unsigned int version);

   // Explicit serialize methods for Util namespace types.

   /*
   * Load a Util::Vector from a MemoryIArchive.
   */
   template <>
   inline void serialize(MemoryIArchive& ar, Vector& data, 
                         const unsigned int version)
   {  ar.unpack(data); }

   /*
   * Load a Util::IntVector from a MemoryIArchive.
   */
   template <>
   inline void serialize(MemoryIArchive& ar, IntVector& data, 
                         const unsigned int version)
   {  ar.unpack(data); }

}
#endif

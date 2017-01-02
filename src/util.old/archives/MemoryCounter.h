#ifndef UTIL_MEMORY_COUNTER_H
#define UTIL_MEMORY_COUNTER_H

/*
* Simpatico - Simulation Package for Polymeric and Molecular Liquids
*
* Copyright 2010 - 2014, The Regents of the University of Minnesota
* Distributed under the terms of the GNU General Public License.
*/

#include "serialize.h"

#include <util/space/Vector.h>
#include <util/space/IntVector.h>

#include <complex>

namespace Util
{

   /**
   * Archive to computed packed size of a sequence of objects, in bytes.
   *
   * This class computes the number of bytes required to pack a sequence 
   * of objects within a MemoryOArchive. The interface is that of a
   * loading Archive, but the << and & operators are overloaded to
   * compute the size required for an object and to increment a size
   * counter, rather than to actually save data. 
   *
   * The size() method returns the number of bytes required to pack all
   * of the objects serialized thus far. The size counter is set to zero
   * upon construction. The clear() method resets the size counter to
   * zero. 
   */
   class MemoryCounter
   {

   public:

      /**
      * Returns true.
      */
      static bool is_saving();

      /**
      * Returns false.
      */
      static bool is_loading();

      /**
      * Constructor.
      */
      MemoryCounter();

      /**
      * Destructor.
      */
      ~MemoryCounter();

      /**
      * Resets the size counter to zero.
      */
      void clear();

      /**
      * Add packed size of one object.
      */
      template <typename T>
      MemoryCounter& operator & (T& data);

      /**
      * Add packed size of one object.
      */
      template <typename T>
      MemoryCounter& operator << (T& data);

      /**
      * Add size of one object in memory.
      *
      * This method just increments the size by sizeof(T). It is
      * appropriate only for primitive C++ variables and POD types
      * for which a bitwise copy is appropriate.
      */
      template <typename T> 
      void count(const T& data);

      /**
      * Compute the size in memory of a C array.
      *
      * This method increments the size by n*sizeof(T). It is
      * appropriate for C arrays of primitive variables and of
      * POD types for which a bitwise copy is appropriate.
      */
      template <typename T> 
      void count(T* array, int n);

      /**
      * Return size required for archive, in Bytes.
      */
      size_t size() const;

   private:

      /// Packed size of sequence of objects, in Bytes.
      size_t size_;

      /// Version id.
      int version_;

      /**
      * Copy constructor (pivate and not implemented).
      */
      MemoryCounter (const MemoryCounter& other);

      /**
      * Assignment (private and not implemented).
      */
      MemoryCounter& operator = (const MemoryCounter& other);

   // friends:

      friend void serialize<>(MemoryCounter&, std::string&, const unsigned int);

   };

   /**
   * Function template to compute memory size of one object.
   */
   template <typename T>
   int memorySize(T& data)
   { 
      MemoryCounter counter;
      counter & data;
      return counter.size();
   }

   // Inline methods

   inline bool MemoryCounter::is_saving()
   { return true; }

   inline bool MemoryCounter::is_loading()
   { return false; }

   /*
   * Return capacity in Bytes.
   */
   inline size_t MemoryCounter::size() const
   {  return size_; }

   /*
   * Compute size of one object, default implementation.
   */
   template <typename T>
   inline MemoryCounter& MemoryCounter::operator & (T& data)
   {  
      serialize(*this, data, version_); 
      return *this;
   }

   /*
   * Compute size of one object, default implementation.
   */
   template <typename T>
   inline MemoryCounter& MemoryCounter::operator << (T& data)
   {   
      serialize(*this, data, version_); 
      return *this;
   }

   // Inline method templates

   /*
   * Compute size of one object.
   */
   template <typename T>
   inline void MemoryCounter::count(const T& data)
   {  size_ += sizeof(T); }

   /*
   * Bitwise pack a C-array of objects of type T.
   */
   template <typename T>
   inline void MemoryCounter::count(T* array, int n)
   {  size_ += n*sizeof(T); }


   // Explicit specializations of serialize function

   // Serialize functions for primitive C++ types

   /*
   * Compute size of an bool.
   */
   template <>
   inline void serialize(MemoryCounter& ar, bool& data, 
                         const unsigned int version)
   {  ar.count(data); }

   /*
   * Compute size of a char.
   */
   template <>
   inline void serialize(MemoryCounter& ar, char& data, 
                         const unsigned int version)
   {  ar.count(data); }

   /*
   * Compute size of an unsigned int.
   */
   template <>
   inline void serialize(MemoryCounter& ar, unsigned int& data, 
                         const unsigned int version)
   {  ar.count(data); }

   /*
   * Compute size of an int.
   */
   template <>
   inline void serialize(MemoryCounter& ar, int& data, 
                         const unsigned int version)
   {  ar.count(data); }

   /*
   * Compute size of an unsigned long int.
   */
   template <>
   inline void serialize(MemoryCounter& ar, unsigned long& data,  
                         const unsigned int version)
   {  ar.count(data); }

   /*
   * Compute size of a long int.
   */
   template <>
   inline void serialize(MemoryCounter& ar, long& data,  
                         const unsigned int version)
   {  ar.count(data); }

   /*
   * Compute size of a float.
   */
   template <>
   inline void serialize(MemoryCounter& ar, float& data, 
                         const unsigned int version)
   {  ar.count(data); }

   /*
   * Compute size of a double.
   */
   template <>
   inline void serialize(MemoryCounter& ar, double& data, 
                         const unsigned int version)
   {  ar.count(data); }

   // Serialize functions for std library types

   /*
   * Compute size of a std::complex<float>.
   */
   template <>
   inline void serialize(MemoryCounter& ar, std::complex<float>& data, 
                         const unsigned int version)
   {  ar.count(data); }

   /*
   * Compute size of a std::complex<double>.
   */
   template <>
   inline void serialize(MemoryCounter& ar, std::complex<double>& data, 
                         const unsigned int version)
   {  ar.count(data); }

   /*
   * Compute size of a std::string.
   */
   template <>
   inline void serialize(MemoryCounter& ar, std::string& data, 
                         const unsigned int version)
   {  ar.size_ += sizeof(size_t) + (data.size() + 1)*sizeof(char); }

   // Serialize functions for Util namespace types

   /*
   * Compute size of a Util::Vector.
   */
   template <>
   inline void serialize(MemoryCounter& ar, Vector& data, 
                         const unsigned int version)
   {  ar.count(data); }

   /*
   * Compute size of a Util::IntVector.
   */
   template <>
   inline void serialize(MemoryCounter& ar, IntVector& data, 
                         const unsigned int version)
   {  ar.count(data); }

}
#endif

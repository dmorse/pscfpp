#ifndef UTIL_BINARY_FILE_I_ARCHIVE_H
#define UTIL_BINARY_FILE_I_ARCHIVE_H

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
#include <iostream>

namespace Util
{

   /**
   * Saving archive for binary istream.
   *
   * \ingroup Serialize_Module
   */
   class BinaryFileIArchive
   {

   public:

      /// Returns true;
      static bool is_saving();

      /// Returns false;
      static bool is_loading();

      /**
      * Constructor.
      */
      BinaryFileIArchive();

      /**
      * Constructor.
      *
      * \param filename name of file to open for reading.
      */
      BinaryFileIArchive(std::string filename);

      /**
      * Constructor.
      *
      * \param file output file
      */
      BinaryFileIArchive(std::ifstream& file);

      /**
      * Destructor.
      */
      virtual ~BinaryFileIArchive();

      /**
      * Get the underlying ifstream by reference.
      */
      std::ifstream& file();

      /**
      * Read one object.
      */
      template <typename T>
      BinaryFileIArchive& operator & (T& data);

      /**
      * Read one object.
      */
      template <typename T>
      BinaryFileIArchive& operator >> (T& data);

      /**
      * Unpack a single T object.
      */
      template <typename T> 
      void unpack(T& data);

      /**
      * Unpack a C array.
      *
      * \param array pointer to array (or first element)
      * \param n number of elements
      */
      template <typename T> 
      void unpack(T* array, int n);

      /**
      * Unpack a 2D C array.
      *
      * This unpacks the elements of an m x n logical array into
      * a physical 2D C array of type array[][np], where np is 
      * the physical length of a row, i.e., the amount of memory 
      * allocated per row.
      *
      * \param array pointer to first row
      * \param m number of rows
      * \param n logical number of columns
      * \param np physical number of columns
      */
      template <typename T> 
      void unpack(T* array, int m, int n, int np);

   private:

      /// Pointer to output stream file.
      std::ifstream* filePtr_;

      /// Archive version id.
      unsigned int  version_;

      /// Was the associated file created by this object?
      bool createdFile_;

   };

   // Inline methods

   inline bool BinaryFileIArchive::is_saving()
   {  return false; }

   inline bool BinaryFileIArchive::is_loading()
   {  return true; }

   // Inline non-static method templates

   /*
   * Read one object.
   */
   template <typename T>
   inline BinaryFileIArchive& BinaryFileIArchive::operator & (T& data)
   {  
      serialize(*this, data, version_); 
      return *this;
   }

   /*
   * Read one object.
   */
   template <typename T>
   inline BinaryFileIArchive& BinaryFileIArchive::operator >> (T& data)
   {   
      serialize(*this, data, version_); 
      return *this;
   }

   // Method templates

   /*
   * Load a single object of type T.
   */
   template <typename T>
   inline void BinaryFileIArchive::unpack(T& data)
   {  filePtr_->read( (char*)(&data), sizeof(T) ); }

   /*
   * Load a C-array of objects of type T.
   */
   template <typename T>
   inline void BinaryFileIArchive::unpack(T* array, int n)
   {
      for (int i=0; i < n; ++i) {
         filePtr_->read( (char*)(&array[i]), sizeof(T));
      }
   }

   /*
   * Unpack a 2D C-array of objects of type T.
   */
   template <typename T>
   inline void BinaryFileIArchive::unpack(T* array, int m, int n, int np)
   {
      int i, j;
      for (i = 0; i < m; ++i) {
         for (j = 0; j < n; ++j) {
            filePtr_->read( (char*)(&array[i*np + j]), sizeof(T));
         }
      }
   }

   // Explicit serialize functions for primitive types

   /*
   * Load a bool from a BinaryFileIArchive.
   */
   template <>
   inline void serialize(BinaryFileIArchive& ar, bool& data, 
                         const unsigned int version)
   {  ar.unpack(data); }

   /*
   * Load a char from a BinaryFileIArchive.
   */
   template <>
   inline void serialize(BinaryFileIArchive& ar, char& data, 
                         const unsigned int version)
   {  ar.unpack(data); }

   /*
   * Load an unsigned int from a BinaryFileIArchive.
   */
   template <>
   inline void serialize(BinaryFileIArchive& ar, unsigned int& data, 
                         const unsigned int version)
   {  ar.unpack(data); }

   /*
   * Load an int from a BinaryFileIArchive.
   */
   template <>
   inline void serialize(BinaryFileIArchive& ar, int& data, 
                         const unsigned int version)
   {  ar.unpack(data); }

   /*
   * Load an unsigned long int from a BinaryFileIArchive.
   */
   template <>
   inline void serialize(BinaryFileIArchive& ar, unsigned long& data,  
                         const unsigned int version)
   {  ar.unpack(data); }

   /*
   * Load a long int from a BinaryFileIArchive.
   */
   template <>
   inline void serialize(BinaryFileIArchive& ar, long& data,  
                         const unsigned int version)
   {  ar.unpack(data); }

   /*
   * Load a float from a BinaryFileIArchive.
   */
   template <>
   inline void serialize(BinaryFileIArchive& ar, float& data, 
                         const unsigned int version)
   {  ar.unpack(data); }

   /*
   * Load a double from a BinaryFileIArchive.
   */
   template <>
   inline void serialize(BinaryFileIArchive& ar, double& data, 
                         const unsigned int version)
   {  ar.unpack(data); }

   /*
   * Load a std::vector from a BinaryFileIArchive.
   */
   template <typename T>
   void serialize(BinaryFileIArchive& ar, std::vector<T>& data, 
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

   // Explicit serialize functions for std library types

   /*
   * Load a std::complex<float> from a BinaryFileIArchive.
   */
   template <>
   inline 
   void serialize(BinaryFileIArchive& ar, std::complex<float>& data, 
                  const unsigned int version)
   {  ar.unpack(data); }

   /*
   * Load a std::complex<double> from a BinaryFileIArchive.
   */
   template <>
   inline 
   void serialize(BinaryFileIArchive& ar, std::complex<double>& data, 
                  const unsigned int version)
   {  ar.unpack(data); }

   /*
   * Load a std::string from a BinaryFileIArchive.
   */
   template <>
   void serialize(BinaryFileIArchive& ar, std::string& data, 
                  const unsigned int version);

   // Explicit serialize functions for namespace Util types

   /*
   * Load a Util::Vector from a BinaryFileIArchive.
   */
   template <>
   inline void serialize(BinaryFileIArchive& ar, Vector& data, 
                         const unsigned int version)
   {  ar.unpack(data); } 

   /*
   * Load a Util::IntVector from a BinaryFileIArchive.
   */
   template <>
   inline void serialize(BinaryFileIArchive& ar, IntVector& data, 
                         const unsigned int version)
   {  ar.unpack(data); }

}
#endif

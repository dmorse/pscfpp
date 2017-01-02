#ifndef UTIL_BINARY_FILE_O_ARCHIVE_H
#define UTIL_BINARY_FILE_O_ARCHIVE_H

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
#include <vector>
#include <string>
#include <iostream>

namespace Util
{

   /**
   * Saving / output archive for binary ostream.
   *
   * \ingroup Serialize_Module
   */
   class BinaryFileOArchive
   {

   public:

      /// Returns true;
      static bool is_saving();

      /// Returns false;
      static bool is_loading();

      /**
      * Constructor.
      */
      BinaryFileOArchive();

      /**
      * Constructor.
      *
      * \param filename name of file to open for reading.
      */
      BinaryFileOArchive(std::string filename);

      /**
      * Constructor.
      *
      * \param file output file
      */
      BinaryFileOArchive(std::ofstream& file);

      /**
      * Destructor.
      */
      virtual ~BinaryFileOArchive();

      /**
      * Get the underlying ifstream by reference.
      */
      std::ofstream& file();

      /**
      * Save one object.
      */
      template <typename T>
      BinaryFileOArchive& operator & (T& data);

      /**
      * Save one object.
      */
      template <typename T>
      BinaryFileOArchive& operator << (T& data);

      /**
      * Pack one object of type T.
      */
      template <typename T> 
      void pack(const T& data);

      /**
      * Pack a C array.
      * 
      * \param array address of first element
      * \param n     number of elements
      */
      template <typename T> 
      void pack(const T* array, int n);

      /**
      * Pack a 2D C array.
      *
      * This packs m rows of length n within a 2D C array allocated
      * as array[][np], where np is the physical length of one row.
      *
      * \param array pointer to [0][0] element in 2D array
      * \param m     number of rows
      * \param n     logical number of columns
      * \param np    physical number of columns
      */
      template <typename T> 
      void pack(const T* array, int m, int n, int np);

   private:

      /// Pointer to output file.
      std::ofstream* filePtr_;

      /// Archive version id.
      unsigned int  version_;

      /// Did this object instantiated the associated file?
      bool createdFile_;

   };

   // Inline methods

   inline bool BinaryFileOArchive::is_saving()
   {  return true; }

   inline bool BinaryFileOArchive::is_loading()
   {  return false; }

   // Inline non-static methods

   /*
   * Write one object.
   */
   template <typename T>
   inline BinaryFileOArchive& BinaryFileOArchive::operator & (T& data)
   {   
      serialize(*this, data, version_); 
      return *this;
   }

   /*
   * Write one object.
   */
   template <typename T>
   inline BinaryFileOArchive& BinaryFileOArchive::operator << (T& data)
   {  
      serialize(*this, data, version_);
      return *this;
   }

   // Method templates

   /*
   * Bitwise pack a single object of type T.
   */
   template <typename T>
   inline void BinaryFileOArchive::pack(const T& data)
   {  filePtr_->write( (char*)(&data), sizeof(T)); }

   /*
   * Bitwise pack a C-array of objects of type T.
   */
   template <typename T>
   inline void BinaryFileOArchive::pack(const T* array, int n)
   {
      for (int i=0; i < n; ++i) {
         filePtr_->write( (char*)(&array[i]), sizeof(T));
      }
   }

   /*
   * Bitwise pack a 2D C-array of objects of type T.
   */
   template <typename T>
   inline void BinaryFileOArchive::pack(const T* array, int m, int n, int np)
   {
      int i, j;
      for (i=0; i < m; ++i) {
         for (j=0; j < n; ++j) {
            filePtr_->write( (char*)(&array[i*np + j]), sizeof(T));
         }
      }
   }

   // Explicit serialize functions for primitive types

   /*
   * Save a bool to a BinaryFileOArchive.
   */
   template <>
   inline void serialize(BinaryFileOArchive& ar, bool& data, 
                         const unsigned int version)
   {  ar.pack(data); }

   /*
   * Save a char to a BinaryFileOArchive.
   */
   template <>
   inline void serialize(BinaryFileOArchive& ar, char& data, 
                         const unsigned int version)
   {  ar.pack(data); }

   /*
   * Save an unsigned int to a BinaryFileOArchive.
   */
   template <>
   inline void serialize(BinaryFileOArchive& ar, unsigned int& data, 
                         const unsigned int version)
   {  ar.pack(data); }

   /*
   * Save an int to a BinaryFileOArchive.
   */
   template <>
   inline void serialize(BinaryFileOArchive& ar, int& data, 
                         const unsigned int version)
   {  ar.pack(data); }

   /*
   * Save an unsigned long int to a BinaryFileOArchive.
   */
   template <>
   inline void serialize(BinaryFileOArchive& ar, unsigned long& data,  
                         const unsigned int version)
   {  ar.pack(data); }

   /*
   * Save a long int to a BinaryFileOArchive.
   */
   template <>
   inline void serialize(BinaryFileOArchive& ar, long& data,  
                         const unsigned int version)
   {  ar.pack(data); }

   /*
   * Save a float to a BinaryFileOArchive.
   */
   template <>
   inline void serialize(BinaryFileOArchive& ar, float& data, 
                         const unsigned int version)
   {  ar.pack(data); }

   /*
   * Save an double to a BinaryFileOArchive.
   */
   template <>
   inline void serialize(BinaryFileOArchive& ar, double& data, 
                         const unsigned int version)
   {  ar.pack(data); }

   /*
   * Save a std::complex<double> to a BinaryFileOArchive.
   */
   template <typename T>
   void serialize(BinaryFileOArchive& ar, std::vector<T>& data, 
                  const unsigned int version)
   {
      size_t size = data.size();
      ar.pack(size);
      for (size_t i = 0; i < size; ++i) {
         ar & data[i];
      }
   }

   // Explicit serialize functions for std library types

   /*
   * Save a std::complex<float> to a BinaryFileOArchive.
   */
   template <>
   inline 
   void serialize(BinaryFileOArchive& ar, std::complex<float>& data, 
                  const unsigned int version)
   {  ar.pack(data); }

   /*
   * Save a std::complex<double> to a BinaryFileOArchive.
   */
   template <>
   inline 
   void serialize(BinaryFileOArchive& ar, std::complex<double>& data, 
                  const unsigned int version)
   {  ar.pack(data); }

   /*
   * Save a std::string to a BinaryFileOArchive.
   */
   template <>
   inline void serialize(BinaryFileOArchive& ar, std::string& data, 
                         const unsigned int version)
   {
      size_t size = data.size() + 1; // the +1 is for the NULL
      ar.pack(size);
      const char* temp = data.c_str();
      ar.pack(temp, size);
   }

   // Explicit serialize functions for namespace Util

   /*
   * Save a Util::Vector to a BinaryFileOArchive.
   */
   template <>
   inline void serialize(BinaryFileOArchive& ar, Vector& data, 
                         const unsigned int version)
   {  ar.pack(data); } 

   /*
   * Save a Util::IntVector to a BinaryFileOArchive.
   */
   template <>
   inline void serialize(BinaryFileOArchive& ar, IntVector& data, 
                         const unsigned int version)
   {  ar.pack(data); }

}
#endif

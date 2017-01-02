#ifndef UTIL_TEXT_FILE_O_ARCHIVE_H
#define UTIL_TEXT_FILE_O_ARCHIVE_H

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
#include <fstream>

namespace Util
{

   /**
   * Saving archive for character based ostream.
   *
   * \ingroup Serialize_Module
   */
   class TextFileOArchive
   {

   public:

      /// Returns true;
      static bool is_saving();

      /// Returns false;
      static bool is_loading();

      /**
      * Constructor.
      */
      TextFileOArchive();

      /**
      * Constructor.
      *
      * \param filename name of file to open for reading.
      */
      TextFileOArchive(std::string filename);

      /**
      * Constructor.
      *
      * \param file output file
      */
      TextFileOArchive(std::ofstream& file);

      /**
      * Destructor.
      */
      virtual ~TextFileOArchive();

      /**
      * Get the underlying ifstream by reference.
      */
      std::ofstream& file();

      /**
      * Save one T object to this archive.
      */
      template <typename T>
      TextFileOArchive& operator & (T& data);

      /**
      * Save one T object to this archive.
      */
      template <typename T>
      TextFileOArchive& operator << (T& data);

      /**
      * Save one T object to this archive.
      *
      * \param data object to be written to file
      */
      template <typename T> 
      void pack(const T& data);

      /**
      * Save a C-array of T objects to this archive.
      *
      * \param array  C array of T objects (pointer to first element)
      * \param n  number of elements
      */
      template <typename T> 
      void pack(const T* array, int n);

      /**
      * Save a 2D C array to this archive.
      *
      * \param array address of first element array[0][0] of 2D array
      * \param m  logical number of rows
      * \param n  logical number of columns
      * \param np physical number of columns
      */
      template <typename T> 
      void pack(const T* array, int m, int n, int np);

   private:

      /// Pointer to output stream file.
      std::ofstream* filePtr_;

      /// Archive version id.
      unsigned int  version_;

      /// Did this object instantiated the associated file?
      bool createdFile_;

   };

   // Inline methods

   inline bool TextFileOArchive::is_saving()
   {  return true; }

   inline bool TextFileOArchive::is_loading()
   {  return false; }

   // Inline non-static method templates.

   /*
   * Save one object.
   */
   template <typename T>
   inline TextFileOArchive& TextFileOArchive::operator & (T& data)
   {   
      serialize(*this, data, version_); 
      return *this;
   }

   /*
   * Save one object to this archive.
   */
   template <typename T>
   inline TextFileOArchive& TextFileOArchive::operator << (T& data)
   {   
      serialize(*this, data, version_); 
      return *this;
   }

   // Method templates

   /*
   * Save a single object of type T.
   */
   template <typename T>
   inline void TextFileOArchive::pack(const T& data)
   {  *filePtr_ << data << std::endl; }

   /*
   * Save a single object of type double.
   */
   template <>
   inline void TextFileOArchive::pack(const double& data)
   {  
      filePtr_->setf(std::ios::scientific);
      filePtr_->width(25);
      filePtr_->precision(17);
      *filePtr_ << data << std::endl; 
   }

   /*
   * Save a C-array of objects of type T.
   */
   template <typename T>
   inline void TextFileOArchive::pack(const T* array, int n)
   {
      for (int i=0; i < n; ++i) {
        *filePtr_ << array[i] << "  ";
      }
      *filePtr_ << std::endl;
   }

   /*
   * Pack a 2D C-array of objects of type T, from array type T array[][np].
   */
   template <typename T>
   inline void TextFileOArchive::pack(const T* array, int m, int n, int np)
   {
      int i, j;
      for (i=0; i < m; ++i) {
         for (j=0; j < n; ++j) {
            *filePtr_ << array[i*np + j] << "  ";
         }
         *filePtr_ << std::endl;
      }
   }

   /*
   * Pack a C-array of double values.
   */
   template <>
   inline void TextFileOArchive::pack(const double* array, int n)
   {
      filePtr_->setf(std::ios::scientific);
      filePtr_->precision(16);
      for (int i=0; i < n; ++i) {
        filePtr_->width(25);
        *filePtr_ << array[i] << "  ";
      }
      *filePtr_ << std::endl;
   }

   // Explicit serialize functions for primitive types

   /*
   * Save a bool to a TextFileOArchive.
   */
   template <>
   inline void serialize(TextFileOArchive& ar, bool& data, 
                         const unsigned int version)
   {  ar.pack(data); }

   /*
   * Save a char to a TextFileOArchive.
   */
   template <>
   inline void serialize(TextFileOArchive& ar, char& data, 
                         const unsigned int version)
   {  ar.pack(data); }

   /*
   * Save an unsigned int to a TextFileOArchive.
   */
   template <>
   inline void serialize(TextFileOArchive& ar, unsigned int& data, 
                         const unsigned int version)
   {  ar.pack(data); }

   /*
   * Save an int to a TextFileOArchive.
   */
   template <>
   inline void serialize(TextFileOArchive& ar, int& data, 
                         const unsigned int version)
   {  ar.pack(data); }

   /*
   * Save an unsigned long int to a TextFileOArchive.
   */
   template <>
   inline void serialize(TextFileOArchive& ar, unsigned long& data,  
                         const unsigned int version)
   {  ar.pack(data); }

   /*
   * Save a long int to a TextFileOArchive.
   */
   template <>
   inline void serialize(TextFileOArchive& ar, long& data,  
                         const unsigned int version)
   {  ar.pack(data); }

   /*
   * Save a float to a TextFileOArchive.
   */
   template <>
   inline void serialize(TextFileOArchive& ar, float& data, 
                         const unsigned int version)
   {  ar.pack(data); }

   /*
   * Save an double to a TextFileOArchive.
   */
   template <>
   inline void serialize(TextFileOArchive& ar, double& data, 
                         const unsigned int version)
   {  ar.pack(data); }

   /*
   * Save a std::vector to a TextFileOArchive.
   */
   template <typename T>
   void serialize(TextFileOArchive& ar, std::vector<T>& data, 
                  const unsigned int version)
   {
      size_t size = data.size();
      ar.pack(size);
      for (size_t i = 0; i < size; ++i) {
         ar & data[i];
      }
   }

   // Explicit serialize functions for standard library types

   /*
   * Save a std::complex<float> to a TextFileOArchive.
   */
   template <>
   inline 
   void serialize(TextFileOArchive& ar, std::complex<float>& data, 
                  const unsigned int version)
   {  ar.pack(data); }

   /*
   * Save a std::complex<double> to a TextFileOArchive.
   */
   template <>
   inline 
   void serialize(TextFileOArchive& ar, std::complex<double>& data, 
                  const unsigned int version)
   {  ar.pack(data); }

   /*
   * Save a std::string to a TextFileOArchive.
   */
   template <>
   inline void serialize(TextFileOArchive& ar, std::string& data, 
                         const unsigned int version)
   {
      int size = data.size();
      ar.pack(size);
      if (size > 0) {
         ar.pack(data);
      }
   }

   // Explicit serialize functions for namespace Util

   /*
   * Save a Util::Vector to a TextFileOArchive.
   */
   template <>
   inline void serialize(TextFileOArchive& ar, Vector& data, 
                         const unsigned int version)
   {  ar.pack(data); } 

   /*
   * Save a Util::IntVector to a TextFileOArchive.
   */
   template <>
   inline void serialize(TextFileOArchive& ar, IntVector& data, 
                         const unsigned int version)
   {  ar.pack(data); }

}
#endif

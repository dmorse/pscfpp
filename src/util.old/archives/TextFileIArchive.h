#ifndef UTIL_TEXT_FILE_I_ARCHIVE_H
#define UTIL_TEXT_FILE_I_ARCHIVE_H

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
   * Loading archive for text istream.
   *
   * \ingroup Serialize_Module
   */
   class TextFileIArchive
   {

   public:

      /// Returns true;
      static bool is_saving();

      /// Returns false;
      static bool is_loading();

      /**
      * Constructor.
      */
      TextFileIArchive();

      /**
      * Constructor.
      *
      * \param filename name of file to open for reading.
      */
      TextFileIArchive(std::string filename);

      /**
      * Constructor.
      *
      * \param file output file
      */
      TextFileIArchive(std::ifstream& file);

      /**
      * Destructor.
      */
      virtual ~TextFileIArchive();

      /**
      * Get the underlying ifstream by reference.
      */
      std::ifstream& file();

      /**
      * Load one object.
      */
      template <typename T>
      TextFileIArchive& operator & (T& data);

      /**
      * Load one object.
      */
      template <typename T>
      TextFileIArchive& operator >> (T& data);

      /**
      * Load a single T object.
      *
      * \param data object to be loaded from this archive.
      */
      template <typename T> 
      void unpack(T& data);

      /**
      * Load a C-array of T objects.
      *
      * \param array pointer to array of T objecs.
      * \param n     number of elements in array
      */
      template <typename T> 
      void unpack(T* array, int n);

      /**
      * Load a 2D C array.
      *
      * \param array pointer to first row or element
      * \param m  logical number of rows
      * \param n  logical number of columns
      * \param np physical number of columns (elements allocated per row)
      */
      template <typename T> 
      void unpack(T* array, int m, int n, int np);

   private:

      /// Pointer to input text file.
      std::ifstream* filePtr_;

      /// Archive version id.
      unsigned int  version_;

      /// Was the associated file created by this object?
      bool createdFile_;

   };

   // Inline methods

   inline bool TextFileIArchive::is_saving()
   {  return false; }

   inline bool TextFileIArchive::is_loading()
   {  return true; }

   // Inline non-static methods

   /*
   * Load one object.
   */
   template <typename T>
   inline TextFileIArchive& TextFileIArchive::operator & (T& data)
   {   
      serialize(*this, data, version_); 
      return *this;
   }

   /*
   * Load one object.
   */
   template <typename T>
   inline TextFileIArchive& TextFileIArchive::operator >> (T& data)
   {   
      serialize(*this, data, version_); 
      return *this;
   }

   // Method templates

   /*
   * Load a single object of type T.
   */
   template <typename T>
   inline void TextFileIArchive::unpack(T& data)
   {  *filePtr_ >> data; }

   /*
   * Load a C array.
   */
   template <typename T>
   inline void TextFileIArchive::unpack(T* array, int n)
   {
      for (int i=0; i < n; ++i) {
         *filePtr_ >> array[i];
      }
   }

   /*
   * Bitwise pack a 2D C-array of objects of type T.
   */
   template <typename T>
   void TextFileIArchive::unpack(T* array, int m, int n, int np)
   {
      int i, j;
      for (i = 0; i < m; ++i) {
         for (j = 0; j < n; ++j) {
            *filePtr_ >> array[i*np + j];
         }
      }
   }

   /*
   * Load a single char.
   */
   template <>
   inline void TextFileIArchive::unpack(char& data)
   {   filePtr_->get(data); }

   /*
   * Load a C-array of characters.
   */
   template <>
   inline void TextFileIArchive::unpack(char* array, int n)
   {
      filePtr_->get(array, n+1,'\0');
   }

   // Explicit serialize functions for primitive types

   /*
   * Load a bool from a TextFileIArchive.
   */
   template <>
   inline void serialize(TextFileIArchive& ar, bool& data, 
                         const unsigned int version)
   {  ar.unpack(data); }

   /*
   * Load a char from a TextFileIArchive.
   */
   template <>
   inline void serialize(TextFileIArchive& ar, char& data, 
                         const unsigned int version)
   {  ar.unpack(data); }

   /*
   * Load an unsigned int from a TextFileIArchive.
   */
   template <>
   inline void serialize(TextFileIArchive& ar, unsigned int& data, 
                         const unsigned int version)
   {  ar.unpack(data); }

   /*
   * Load an int from a TextFileIArchive.
   */
   template <>
   inline void serialize(TextFileIArchive& ar, int& data, 
                         const unsigned int version)
   {  ar.unpack(data); }

   /*
   * Load an unsigned long int from a TextFileIArchive.
   */
   template <>
   inline void serialize(TextFileIArchive& ar, unsigned long& data,  
                         const unsigned int version)
   {  ar.unpack(data); }

   /*
   * Load a long int from a TextFileIArchive.
   */
   template <>
   inline void serialize(TextFileIArchive& ar, long& data,  
                         const unsigned int version)
   {  ar.unpack(data); }

   /*
   * Load a float from a TextFileIArchive.
   */
   template <>
   inline void serialize(TextFileIArchive& ar, float& data, 
                         const unsigned int version)
   {  ar.unpack(data); }

   /*
   * Load a double from a TextFileIArchive.
   */
   template <>
   inline void serialize(TextFileIArchive& ar, double& data, 
                         const unsigned int version)
   {  ar.unpack(data); }

   /*
   * Load a std::vector from a TextFileIArchive.
   */
   template <typename T>
   void serialize(TextFileIArchive& ar, std::vector<T>& data, 
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
   * Load a std::complex<float> from a TextFileIArchive.
   */
   template <>
   inline 
   void serialize(TextFileIArchive& ar, std::complex<float>& data, 
                  const unsigned int version)
   {  ar.unpack(data); }

   /*
   * Load a std::complex<double> from a TextFileIArchive.
   */
   template <>
   inline 
   void serialize(TextFileIArchive& ar, std::complex<double>& data, 
                  const unsigned int version)
   {  ar.unpack(data); }

   /*
   * Load a std::string from a TextFileIArchive.
   */
   template <>
   void serialize(TextFileIArchive& ar, std::string& data, 
                  const unsigned int version);

   // Explicit serialize functions for namespace Util types

   /*
   * Load a Util::Vector from a TextFileIArchive.
   */
   template <>
   inline void serialize(TextFileIArchive& ar, Vector& data, 
                         const unsigned int version)
   {  ar.unpack(data); } 

   /*
   * Load a Util::IntVector from a TextFileIArchive.
   */
   template <>
   inline void serialize(TextFileIArchive& ar, IntVector& data, 
                         const unsigned int version)
   {  ar.unpack(data); }

}
#endif

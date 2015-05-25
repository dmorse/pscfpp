#ifndef UTIL_XDR_FILE_I_ARCHIVE_H
#define UTIL_XDR_FILE_I_ARCHIVE_H

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
#include <vector>
#include <string>

#include <rpc/rpc.h>
#include <rpc/xdr.h>
#include <stdio.h>

namespace Util
{

   /**
   * Loading / input archive for binary XDR file. 
   *
   * XDR is a standard protocol for writing and reading
   * binary in a portable format.  This archive saves
   * data to an associated file in XDR format. It depends
   * on the unix xdr library <rpc/xdr.h>. Because this
   * library is written in C (not C++), this archive 
   * uses a standard C library file handle, not a C++
   * iostream.
   *
   * \ingroup Serialize_Module
   */
   class XdrFileIArchive
   {

   public:

      /// Returns false.
      static bool is_saving();

      /// Returns true.
      static bool is_loading();

      /**
      * Constructor.
      */
      XdrFileIArchive();

      /**
      * Constructor.
      *
      * \param filename name of file to open for reading.
      */
      XdrFileIArchive(std::string filename);

      /**
      * Constructor.
      *
      * \param file output file
      */
      XdrFileIArchive(std::ofstream& file);

      /**
      * Destructor.
      */
      virtual ~XdrFileIArchive();

      /**
      * Initialize by associating with an open file.
      *
      * \param file C library file handle, must be open for reading.
      */
      void init(FILE* file);

      /**
      * Load one object.
      */
      template <typename T>
      XdrFileIArchive& operator & (T& data);

      /**
      * Load one object.
      */
      template <typename T>
      XdrFileIArchive& operator >> (T& data);

      /**
      * Get the underlying file handle.
      */
      FILE* file();

      /**
      * Get a pointer to the enclosed XDR object.
      */
      XDR* xdrPtr();

   private:

      /// XDR (external data representation) handle.
      XDR  xdr_;      

      /// Pointer to standard C library file handle.
      FILE* filePtr_;       
 
      /// Archive version id.
      unsigned int  version_;

      /// Did this object create the associated file?
      bool createdFile_;

   };

   // Inline methods

   inline bool XdrFileIArchive::is_saving()
   {  return false; }

   inline bool XdrFileIArchive::is_loading()
   {  return true; }

   // Inline non-static methods

   /*
   * Write one object.
   */
   template <typename T>
   inline XdrFileIArchive& XdrFileIArchive::operator & (T& data)
   {   
      serialize(*this, data, version_); 
      return *this;
   }

   /*
   * Write one object.
   */
   template <typename T>
   inline XdrFileIArchive& XdrFileIArchive::operator >> (T& data)
   {  
      serialize(*this, data, version_);
      return *this;
   }

   /*
   * Get the underlying file handle.
   */
   inline FILE* XdrFileIArchive::file()
   {  return filePtr_; }

   /*
   * Get a pointer to XDR object.
   */
   inline XDR* XdrFileIArchive::xdrPtr()
   {  return &xdr_; }

   // Explicit serialize functions for primitive types

   /*
   * Load a bool from a XdrFileIArchive.
   */
   template <>
   inline void serialize(XdrFileIArchive& ar, bool& data, 
                         const unsigned int version)
   {
      bool_t temp = data;
      xdr_bool(ar.xdrPtr(), &temp);  
   }

   /*
   * Load a char from a XdrFileIArchive.
   */
   template <>
   inline void serialize(XdrFileIArchive& ar, char& data, 
                         const unsigned int version)
   {  xdr_char(ar.xdrPtr(), &data);  }

   /*
   * Load an unsigned int from a XdrFileIArchive.
   */
   template <>
   inline void serialize(XdrFileIArchive& ar, unsigned int& data, 
                         const unsigned int version)
   {  xdr_u_int(ar.xdrPtr(), &data);  }

   /*
   * Load an int from a XdrFileIArchive.
   */
   template <>
   inline void serialize(XdrFileIArchive& ar, int& data, 
                         const unsigned int version)
   {  xdr_int(ar.xdrPtr(), &data);  }

   #if 0
   /*
   * Load an unsigned long int from a XdrFileIArchive.
   */
   template <>
   inline void serialize(XdrFileIArchive& ar, unsigned long& data,  
                         const unsigned int version)
   {  xdr_u_long(ar.xdrPtr(), &data);  }

   /*
   * Load a long int from a XdrFileIArchive.
   */
   template <>
   inline void serialize(XdrFileIArchive& ar, long& data,  
                         const unsigned int version)
   {  xdr_long(ar.xdrPtr(), &data);  }
   #endif

   /*
   * Load a float from a XdrFileIArchive.
   */
   template <>
   inline void serialize(XdrFileIArchive& ar, float& data, 
                         const unsigned int version)
   {  xdr_float(ar.xdrPtr(), &data);  }

   /*
   * Load a double from a XdrFileIArchive.
   */
   template <>
   inline void serialize(XdrFileIArchive& ar, double& data, 
                         const unsigned int version)
   {  xdr_double(ar.xdrPtr(), &data);  }

   /*
   * Load a std::vector from a XdrFileIArchive.
   */
   template <typename T>
   void serialize(XdrFileIArchive& ar, std::vector<T>& data, 
                  const unsigned int version)
   {
      unsigned int size = data.size();
      xdr_u_int(ar.xdrPtr(), &size);
      for (size_t i = 0; i < size; ++i) {
         ar & data[i];
      }
   }

   // Explicit serialize functions for std library types

   /*
   * Load a std::complex<float> from a XdrFileIArchive.
   */
   template <>
   inline 
   void serialize(XdrFileIArchive& ar, std::complex<float>& data, 
                  const unsigned int version)
   {
      float a;
      float b;
      xdr_float(ar.xdrPtr(), &a);  
      xdr_float(ar.xdrPtr(), &b);  
      data = std::complex<float>(a, b);
   }

   /*
   * Load a std::complex<double> from a XdrFileIArchive.
   */
   template <>
   inline 
   void serialize(XdrFileIArchive& ar, std::complex<double>& data, 
                  const unsigned int version)
   {
      double a;
      double b;
      xdr_double(ar.xdrPtr(), &a);  
      xdr_double(ar.xdrPtr(), &b);  
      data = std::complex<double>(a, b);
   }

   /*
   * Load a std::string from a XdrFileIArchive.
   */
   template <>
   void serialize(XdrFileIArchive& ar, std::string& data, 
                         const unsigned int version);

   // Explicit serialize functions for namespace Util

   /*
   * Load a Util::Vector from a XdrFileIArchive.
   */
   template <>
   inline void serialize(XdrFileIArchive& ar, Vector& data, 
                         const unsigned int version)
   {
      ar & data[0]; 
      ar & data[1]; 
      ar & data[2]; 
   } 

   /*
   * Load a Util::IntVector from a XdrFileIArchive.
   */
   template <>
   inline void serialize(XdrFileIArchive& ar, IntVector& data, 
                         const unsigned int version)
   {  
      ar & data[0];
      ar & data[1];
      ar & data[2];
   }

}
#endif

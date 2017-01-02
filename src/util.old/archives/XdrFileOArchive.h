#ifndef UTIL_XDR_FILE_O_ARCHIVE_H
#define UTIL_XDR_FILE_O_ARCHIVE_H

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
   * Saving / output archive for binary XDR file.
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
   class XdrFileOArchive
   {

   public:

      /// Returns true;
      static bool is_saving();

      /// Returns false;
      static bool is_loading();

      /**
      * Constructor.
      */
      XdrFileOArchive();

      /**
      * Constructor.
      *
      * \param filename name of file to open for reading.
      */
      XdrFileOArchive(std::string filename);

      /**
      * Destructor.
      */
      virtual ~XdrFileOArchive();

      /**
      * Associate with an open file and initialize.
      *
      * \param file C file handle, must be open for writing.
      */
      void init(FILE* file);

      /**
      * Get the underlying ifstream by reference.
      */
      FILE* file();

      /**
      * Save one object.
      */
      template <typename T>
      XdrFileOArchive& operator & (T& data);

      /**
      * Save one object.
      */
      template <typename T>
      XdrFileOArchive& operator << (T& data);

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

   inline bool XdrFileOArchive::is_saving()
   {  return true; }

   inline bool XdrFileOArchive::is_loading()
   {  return false; }

   // Inline non-static methods

   /*
   * Write one object.
   */
   template <typename T>
   inline XdrFileOArchive& XdrFileOArchive::operator & (T& data)
   {   
      serialize(*this, data, version_); 
      return *this;
   }

   /*
   * Write one object.
   */
   template <typename T>
   inline XdrFileOArchive& XdrFileOArchive::operator << (T& data)
   {  
      serialize(*this, data, version_);
      return *this;
   }

   /*
   * Get a C file handle.
   */
   inline FILE* XdrFileOArchive::file()
   {  return filePtr_; }

   /*
   * Get a pointer to XDR object.
   */
   inline XDR* XdrFileOArchive::xdrPtr()
   {  return &xdr_; }

   // Explicit serialize functions for primitive types

   /*
   * Save a bool to a XdrFileOArchive.
   */
   template <>
   inline void serialize(XdrFileOArchive& ar, bool& data, 
                         const unsigned int version)
   {
      bool_t temp = data;
      xdr_bool(ar.xdrPtr(), &temp);  
   }

   /*
   * Save a char to a XdrFileOArchive.
   */
   template <>
   inline void serialize(XdrFileOArchive& ar, char& data, 
                         const unsigned int version)
   {  xdr_char(ar.xdrPtr(), &data);  }

   /*
   * Save an unsigned int to a XdrFileOArchive.
   */
   template <>
   inline void serialize(XdrFileOArchive& ar, unsigned int& data, 
                         const unsigned int version)
   {  xdr_u_int(ar.xdrPtr(), &data);  }

   /*
   * Save an int to a XdrFileOArchive.
   */
   template <>
   inline void serialize(XdrFileOArchive& ar, int& data, 
                         const unsigned int version)
   {  xdr_int(ar.xdrPtr(), &data);  }

   #if 0
   /*
   * Save an unsigned long int to a XdrFileOArchive.
   */
   template <>
   inline void serialize(XdrFileOArchive& ar, unsigned long& data,  
                         const unsigned int version)
   {  xdr_u_long(ar.xdrPtr(), &data);  }

   /*
   * Save a long int to a XdrFileOArchive.
   */
   template <>
   inline void serialize(XdrFileOArchive& ar, long& data,  
                         const unsigned int version)
   {  xdr_long(ar.xdrPtr(), &data);  }
   #endif

   /*
   * Save a float to a XdrFileOArchive.
   */
   template <>
   inline void serialize(XdrFileOArchive& ar, float& data, 
                         const unsigned int version)
   {  xdr_float(ar.xdrPtr(), &data);  }

   /*
   * Save a double to a XdrFileOArchive.
   */
   template <>
   inline void serialize(XdrFileOArchive& ar, double& data, 
                         const unsigned int version)
   {  xdr_double(ar.xdrPtr(), &data);  }

   /*
   * Save a std::vector to a XdrFileOArchive.
   */
   template <typename T>
   void serialize(XdrFileOArchive& ar, std::vector<T>& data, 
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
   * Save a std::complex<float> to a XdrFileOArchive.
   */
   template <>
   inline 
   void serialize(XdrFileOArchive& ar, std::complex<float>& data, 
                  const unsigned int version)
   {
      float a = data.real();
      float b = data.imag();
      xdr_float(ar.xdrPtr(), &a);  
      xdr_float(ar.xdrPtr(), &b);  
   }

   /*
   * Save a std::complex<double> to a XdrFileOArchive.
   */
   template <>
   inline 
   void serialize(XdrFileOArchive& ar, std::complex<double>& data, 
                  const unsigned int version)
   {
      double a = data.real();
      double b = data.imag();
      xdr_double(ar.xdrPtr(), &a);  
      xdr_double(ar.xdrPtr(), &b);  
   }

   /*
   * Save a std::string to a XdrFileOArchive.
   */
   template <>
   void serialize(XdrFileOArchive& ar, std::string& data, 
                         const unsigned int version);

   // Explicit serialize functions for namespace Util

   /*
   * Save a Util::Vector to a XdrFileOArchive.
   */
   template <>
   inline void serialize(XdrFileOArchive& ar, Vector& data, 
                         const unsigned int version)
   {
      ar & data[0]; 
      ar & data[1]; 
      ar & data[2]; 
   } 

   /*
   * Save a Util::IntVector to a XdrFileOArchive.
   */
   template <>
   inline void serialize(XdrFileOArchive& ar, IntVector& data, 
                         const unsigned int version)
   {  
      ar & data[0]; 
      ar & data[1]; 
      ar & data[2]; 
   }

}
#endif

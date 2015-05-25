/*
* Simpatico - Simulation Package for Polymeric and Molecular Liquids
*
* Copyright 2010 - 2014, The Regents of the University of Minnesota
* Distributed under the terms of the GNU General Public License.
*/

#include "XdrFileIArchive.h"
#include <string.h>

namespace Util
{

   /*
   * Constructor.
   */
   XdrFileIArchive::XdrFileIArchive()
    : xdr_(),
      filePtr_(0),
      version_(0),
      createdFile_(false)
   {}

   /*
   * Constructor.
   */
   XdrFileIArchive::XdrFileIArchive(std::string filename)
    : xdr_(),
      filePtr_(0),
      version_(0),
      createdFile_(true)
   {
      filePtr_ = fopen(filename.c_str(), "rb"); 
      if (filePtr_ == NULL) {
         std::string msg = "Failure to open C file: ";
         msg += filename.c_str();
         UTIL_THROW(msg.c_str());
      }
      xdrstdio_create(&xdr_, filePtr_, XDR_DECODE);
   }

   /*
   * Destructor.
   */
   XdrFileIArchive::~XdrFileIArchive()
   {}

   /*
   * Initialize if default constructed.
   */
   void XdrFileIArchive::init(FILE* file)
   {
      filePtr_ = file;  
      xdrstdio_create(&xdr_, filePtr_, XDR_DECODE); 
   }

   /*
   * Load a std::string from a XdrFileIArchive.
   */
   template <>
   void serialize(XdrFileIArchive& ar, std::string& data, 
                  const unsigned int version)
   {
      static char* temp = 0;
      static unsigned int tempsize = 0;

      unsigned int size;
      xdr_u_int(ar.xdrPtr(), &size);

      if (temp) {
         if (size > tempsize) {
            tempsize = size;
            delete [] temp;
            temp = new char[tempsize];
         }
      } else {
         tempsize = 256;
         if (size > tempsize) {
            tempsize = size;
         }
         temp = new char[tempsize];
      }
      xdr_string(ar.xdrPtr(), &temp, size);
      data = temp;
   }

}

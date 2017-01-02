/*
* Simpatico - Simulation Package for Polymeric and Molecular Liquids
*
* Copyright 2010 - 2014, The Regents of the University of Minnesota
* Distributed under the terms of the GNU General Public License.
*/

#include "XdrFileOArchive.h"
#include <string.h>

namespace Util
{

   /*
   * Constructor.
   */
   XdrFileOArchive::XdrFileOArchive()
    : xdr_(),
      filePtr_(0),
      version_(0),
      createdFile_(false)
   {}

   /*
   * Constructor.
   */
   XdrFileOArchive::XdrFileOArchive(std::string filename)
    : xdr_(),
      filePtr_(0),
      version_(0),
      createdFile_(true)
   {
      filePtr_ = fopen(filename.c_str(), "wb+"); 
      if (filePtr_ == NULL) {
         std::string msg = "Failure to open C file: ";
         msg += filename.c_str();
         UTIL_THROW(msg.c_str());
      }
      xdrstdio_create(&xdr_, filePtr_, XDR_ENCODE);
   }

   /*
   * Destructor.
   */
   XdrFileOArchive::~XdrFileOArchive()
   {}

   /*
   * Initialize if default constructed.
   */
   void XdrFileOArchive::init(FILE* file)
   {
      filePtr_ = file;  
      xdrstdio_create(&xdr_, filePtr_, XDR_ENCODE); 
   }

   /*
   * Save a std::string to a XdrFileOArchive.
   */
   template <>
   void serialize(XdrFileOArchive& ar, std::string& data, 
                  const unsigned int version)
   {
      static char* temp = 0;
      static unsigned int tempsize = 0;

      unsigned int size = data.size() + 1; // +1 for the '\0'
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
      strcpy(temp, data.c_str());
      xdr_string(ar.xdrPtr(), &temp, size);
   }

}

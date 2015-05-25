/*
* Simpatico - Simulation Package for Polymeric and Molecular Liquids
*
* Copyright 2010 - 2014, The Regents of the University of Minnesota
* Distributed under the terms of the GNU General Public License.
*/

#include "TextFileIArchive.h"

#include <vector>

namespace Util
{

   /*
   * Constructor.
   */
   TextFileIArchive::TextFileIArchive()
    : filePtr_(0),
      version_(0),
      createdFile_(true)
   {  filePtr_ = new std::ifstream(); }

   /*
   * Constructor.
   */
   TextFileIArchive::TextFileIArchive(std::string filename)
    : filePtr_(0),
      version_(0),
      createdFile_(true)
   {  filePtr_ = new std::ifstream(filename.c_str()); }


   /*
   * Constructor.
   */
   TextFileIArchive::TextFileIArchive(std::ifstream& file)
    : filePtr_(&file),
      version_(0),
      createdFile_(false)
   {  
      if (!file.is_open()) {
         UTIL_THROW("File not open");
      }  
   }


   /*
   * Destructor.
   */
   TextFileIArchive::~TextFileIArchive()
   {  
      if (filePtr_ && createdFile_) {  
         delete filePtr_; 
      }
   }

   /*
   * Return underlying file by reference.
   */
   std::ifstream& TextFileIArchive::file()
   {  return *filePtr_; }

   /*
   * Load a std::string from TextFileIArchive.
   */
   template <>
   void serialize(TextFileIArchive& ar, std::string& data, 
                  const unsigned int version)
   {
      size_t size;
      ar.unpack(size);
      if (size > 0) {
         static std::vector<char> charvec;
         if (size > charvec.capacity()) {
            charvec.reserve(size + 8);
         }
         // Read endline character after size
         char   endline;
         ar.unpack(endline);
         // Read actual string.
         ar.unpack(&charvec[0], size);
         data = &charvec[0];
      } else {
         data = "";
      }
   }

}

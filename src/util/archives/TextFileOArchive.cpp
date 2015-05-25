/*
* Simpatico - Simulation Package for Polymeric and Molecular Liquids
*
* Copyright 2010 - 2014, The Regents of the University of Minnesota
* Distributed under the terms of the GNU General Public License.
*/

#include "TextFileOArchive.h"

namespace Util
{

   /*
   * Constructor.
   */
   TextFileOArchive::TextFileOArchive()
    : filePtr_(0),
      version_(0),
      createdFile_(true)
   {  filePtr_ = new std::ofstream(); }

   /*
   * Constructor.
   */
   TextFileOArchive::TextFileOArchive(std::string filename)
    : filePtr_(0),
      version_(0),
      createdFile_(true)
   {  filePtr_ = new std::ofstream(filename.c_str()); }


   /*
   * Constructor.
   */
   TextFileOArchive::TextFileOArchive(std::ofstream& file)
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
   TextFileOArchive::~TextFileOArchive()
   {
      if (filePtr_ && createdFile_) {  
         delete filePtr_; 
      }
   }

   /*
   * Return underlying file by reference.
   */
   std::ofstream& TextFileOArchive::file()
   {  return *filePtr_; }

}

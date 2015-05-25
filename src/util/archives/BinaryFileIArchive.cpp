/*
* Simpatico - Simulation Package for Polymeric and Molecular Liquids
*
* Copyright 2010 - 2014, The Regents of the University of Minnesota
* Distributed under the terms of the GNU General Public License.
*/

#include "BinaryFileIArchive.h"

#include <vector>

namespace Util
{

   /*
   * Constructor.
   */
   BinaryFileIArchive::BinaryFileIArchive()
    : filePtr_(0),
      version_(0),
      createdFile_(true)
   {  filePtr_ = new std::ifstream(); }

   /*
   * Constructor.
   */
   BinaryFileIArchive::BinaryFileIArchive(std::string filename)
    : filePtr_(0),
      version_(0),
      createdFile_(true)
   {  filePtr_ = new std::ifstream(filename.c_str()); }

   /*
   * Constructor.
   */
   BinaryFileIArchive::BinaryFileIArchive(std::ifstream& file)
    : filePtr_(&file),
      version_(0),
      createdFile_(false)
   {}

   /*
   * Destructor.
   */
   BinaryFileIArchive::~BinaryFileIArchive()
   {
      if (filePtr_ && createdFile_) {  
         delete filePtr_; 
      }
   }

   /*
   * Return underlying file by reference.
   */
   std::ifstream& BinaryFileIArchive::file()
   {  return *filePtr_; }

   /*
   * Load a std::string from BinaryFileIArchive.
   */
   template <>
   void serialize(BinaryFileIArchive& ar, std::string& data, 
                  const unsigned int version)
   {
      static std::vector<char> charvec;
      size_t size;
      ar.unpack(size);
      if (size > charvec.capacity()) {
         charvec.reserve(size + 8);
      }
      ar.unpack(&charvec[0], size);
      data = &charvec[0];
   }


}

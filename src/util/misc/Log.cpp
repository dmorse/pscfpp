/*
* Simpatico - Simulation Package for Polymeric and Molecular Liquids
*
* Copyright 2010 - 2014, The Regents of the University of Minnesota
* Distributed under the terms of the GNU General Public License.
*/

#include "Log.h"

namespace Util
{

   // Define and nullify streamPtr_ and filePtr_
   std::ostream*  Log::streamPtr_ = 0;
   std::ofstream* Log::filePtr_ = 0;

   /*
   * This static method exists to guarantee initialization of static 
   * constants that are defined in the same file.  Call it somewhere 
   * in the program to guarantee that the contents of this file will 
   * be linked, rather than optimized away. It may only be called once.
   */
   void Log::initStatic() 
   {
      static int nCall = 0;
      if (nCall == 0) {
         Log::streamPtr_ = 0;
         Log::filePtr_ = 0;
      }
      ++nCall;
   }

   /*
   *  Set the log file to a new file.
   */
   void Log::setFile(std::ofstream& file)
   {  
      streamPtr_ = &file; 
      filePtr_   = &file; 
   }

   /*
   * Close the log file, if any.
   */
   void Log::close() 
   {  
      if (filePtr_) {
         filePtr_->close();  
      }
      filePtr_ = 0;
      streamPtr_ = &std::cout;
   }

   /*
   * Get the log ostream by reference.
   */
   std::ostream& Log::file() 
   {
      if (!streamPtr_) {
         streamPtr_ = &std::cout;
      }  
      return *streamPtr_; 
   }

} 

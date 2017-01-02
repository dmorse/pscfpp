#ifndef UTIL_LOG_H
#define UTIL_LOG_H

/*
* Simpatico - Simulation Package for Polymeric and Molecular Liquids
*
* Copyright 2010 - 2014, The Regents of the University of Minnesota
* Distributed under the terms of the GNU General Public License.
*/

#include <iostream>
#include <fstream>

namespace Util
{

   /**
   * A static class that holds a log output stream.
   *
   * The Log class has one a static pointer member that points to an
   * ostream that should be used by all other classes to output log 
   * and error messages. This stream is accessed by the file() method. 
   *
   * The log file initialized to point to std::cout. It may be reset 
   * to point to a ofstream file object using the static setFile() 
   * method. 
   *
   * \ingroup Misc_Module
   */
   class Log
   {
   
   public:

      /**
      * Initialize static members.
      */
      static void initStatic();

      /**
      * Set the log ostream to a file.
      *
      * \param file ofstream open for writing.
      */
      static void setFile(std::ofstream& file);
 
      /**
      * Close log file, if any.
      */
      static void close();

      /**
      * Get log ostream by reference.
      */
      static std::ostream& file();

   private:

      /// Pointer to log output stream
      static std::ostream* streamPtr_;

      /// Pointer to log file (if it is a file)
      static std::ofstream* filePtr_;

      /// Default constructor private to prevent instantiation.
      Log();

   };

}
#endif

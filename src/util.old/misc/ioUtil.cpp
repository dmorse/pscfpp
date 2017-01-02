/*
* Simpatico - Simulation Package for Polymeric and Molecular Liquids
*
* Copyright 2010 - 2014, The Regents of the University of Minnesota
* Distributed under the terms of the GNU General Public License.
*/

#include <util/global.h>
#include "ioUtil.h"
#include "Log.h"

namespace Util
{

   /*
   * Strip trailing whitespace from a string.
   */
   int rStrip(std::string& str)
   {
      size_t found;
      std::string whitespaces(" \t\n\r");
      found = str.find_last_not_of(whitespaces);
      if (found != std::string::npos) {
        str.erase(found + 1);
        return int(found + 1);
      } else {
        str.clear();
        return 0;
      }
   }

   /*
   * Read string, and compare to expected value. 
   *
   * Throw Exception if input value differs from expected value.
   */
   void checkString(std::istream& in, const std::string& expected)
   {
      std::string actual;
      in >> actual;
      if (actual != expected) {
         Log::file() << "Error in checkString"     << std::endl;
         Log::file() << "Expected: " <<  expected  << std::endl;
         Log::file() << "Scanned:  " <<  actual    << std::endl;
         UTIL_THROW("Incorrect string");
      };
   }

   /*
   * Return std::string representation of an integer.
   */
   std::string toString(int n)
   {
      std::stringstream ss;
      ss << n;
      return ss.str();
   }

   /*
   * Transfer the next line of input to a stringstream
   */
   bool getLine(std::istream& in, std::stringstream& line)
   {
      std::string string;
      line.str(string);
      line.clear();
      if (!in.eof()) {
         getline(in, string);
         line.str(string);
         return true;
      } else {
         return false;
      }
   }

   /*
   * Get the next non-blank line of input, strip trailing whitespace.
   */
   bool getNextLine(std::istream& in, std::string& line)
   {
      while (true) {
         if (!in.eof()) {
            getline(in, line);
            rStrip(line);
            if (!line.empty()) {
               // Return true indicates eof not reached
               return true;
            }
         } else {
            line.clear();
            // Return false indicates eof was reached
            return false;
         }
      }
   }

   /*
   * Assign the next non-blank line of input to a stringstream
   */
   bool getNextLine(std::istream& in, std::stringstream& line)
   {
      std::string string;
      line.str(string);
      line.clear();
      while (true) {
         if (!in.eof()) {
            getline(in, string);
            rStrip(string);
            if (!string.empty()) {
               line.str(string);
               // Return true indicates eof not reached
               return true;
            }
         } else {
            // Return false indicates eof reached
            return false;
         }
      }
   }

}

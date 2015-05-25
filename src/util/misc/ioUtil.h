#ifndef UTIL_IO_UTIL_H
#define UTIL_IO_UTIL_H

/*
* Simpatico - Simulation Package for Polymeric and Molecular Liquids
*
* Copyright 2010 - 2014, The Regents of the University of Minnesota
* Distributed under the terms of the GNU General Public License.
*/

#include <string>
#include <iostream>
#include <sstream>

namespace Util
{

   /**
   * Return string representation of an integer.
   *
   * \ingroup Misc_Module
   *
   * \param n integer to be converted.
   */
   std::string toString(int n);

   /**
   * Strip trailing whitespace from a string.
   *
   * \ingroup Misc_Module
   *
   * \param  string string (stripped upon return).
   * \return length of stripped string.
   */
   int rStrip(std::string& string);

   /**
   * Extract string from stream, and compare to expected value. 
   *
   * \throw Exception if input value differs from expected value.
   *
   * \ingroup Misc_Module
   *
   * \param in input stream
   * \param expected expected value of string read from stream
   */
   void checkString(std::istream& in, const std::string& expected);

   /**
   * Read the next line into a stringstream.
   *
   * Variant of std::getline(). Does not strip trailing whitespace.
   *
   * \ingroup Misc_Module
   *
   * \param in input stream from which to read.
   * \param line stringstream containing line, on output.
   */
   bool getLine(std::istream& in, std::stringstream& line);

   /**
   * Read the next non-empty line into a string, strip trailing whitespace.
   *
   * Variant of std::getline() that skips empty lines.
   *
   * \ingroup Misc_Module
   *
   * \param in input stream from which to read.
   * \param line string with next non-empty line, on output.
   * \return true if not end-of-file, false if end-of-file.
   */
   bool getNextLine(std::istream& in, std::string& line);

   /**
   * Read the next non-empty line into a stringstream, strip trailing whitespace.
   *
   * Variant of std::getline() that skips empty lines and uses stringstream.
   *
   * \ingroup Misc_Module
   *
   * \param in   input stream from which to read.
   * \param line stringstream containing next non-empty line, on output.
   * \return true if not end-of-file, false if end-of-file.
   */
   bool getNextLine(std::istream& in, std::stringstream& line);

}
#endif

#ifndef COMMAND_LINE_H
#define COMMAND_LINE_H

/*
* Simpatico - Simulation Package for Polymeric and Molecular Liquids
*
* Copyright 2010 - 2014, The Regents of the University of Minnesota
* Distributed under the terms of the GNU General Public License.
*/

#include "TestException.h"

#include <string>
#include <vector>

/**
* Abstraction of a C array of command line arguments.
*
* \ingroup Test_Module
*/
class CommandLine
{

public:

   /**
   * Constructor.
   */
   CommandLine()
   {  clear(); }

   /**
   * Add a command line argument string.
   */
   void append(const char* arg)
   {  strings_.push_back(std::string(arg)); }

   /**
   * Clear all arguments.
   */
   void clear()
   {   
      argv_.clear(); 
      strings_.clear(); 
      append("");
      // Empty string is a substitute for the executable name.
      // The C standard says that the first argument argv[0] should
      // be the executable name, or an empty string if unknown.
   }

   /**
   * Return number of command line arguments.
   */
   int argc()
   {  return strings_.size(); }

   /**
   * Return array of C-string command line argument strings.
   *
   * \return pointer to C array of null terminated C strings.
   */
   char** argv()
   {
      argv_.clear();
      for (unsigned int i = 0; i < strings_.size(); ++i) {
         argv_.push_back(const_cast<char*>(strings_[i].c_str()));
      }
      return &(argv_[0]);
   }

private:

   std::vector<std::string> strings_;
   std::vector<char*> argv_;

};
#endif

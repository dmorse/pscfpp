/*
* Simpatico - Simulation Package for Polymeric and Molecular Liquids
*
* Copyright 2010 - 2014, The Regents of the University of Minnesota
* Distributed under the terms of the GNU General Public License.
*/

#include "XmlEndTag.h"
#include <util/misc/ioUtil.h>

namespace Util
{

   XmlEndTag::XmlEndTag()
   {}

   XmlEndTag::~XmlEndTag()
   {}

   bool XmlEndTag::match(const std::string& line, int begin)
   {
      label_ = "";
      setString(line, begin);

      // Skip leading white space
      if (isEnd()) return false;
      skip();
      if (isEnd()) return false;

      // Read opening bracket
      if (c() != '<') {
         //std::cout << "Missing opening bracket" << std::endl;
         return false;
      }
      next();
      if (isEnd()) return false;
      if (c() != '/') {
         //std::cout << "Missing slash" << std::endl;
         return false;
      }

      // Read label
      int beginLabel, endLabel;
      next();
      skip();
      if (isEnd()) return false;
      beginLabel = cursor();
      while (c() != '>' && c() != ' ') {
         next();
         if (isEnd()) return false;
      }
      endLabel = cursor();
      label_ = string().substr(beginLabel, endLabel - beginLabel);
    
      skip();
      if (isEnd()) return false;
      if (c() == '>') {
         return true;
      } else {
         return false;
      }

   }
   
   void XmlEndTag::match(const std::string expected,
                         const std::string& line, int begin)
   {
      if (!match(line, begin)) {
         Log::file() << "line = " << line << std::endl;
         UTIL_THROW("No end tag");
      } 
      if (label() != expected) {
         Log::file() << "line     = " << line << std::endl;
         Log::file() << "expected = " << expected << std::endl; 
         Log::file() << "label    = " << label()  << std::endl; 
         UTIL_THROW("Incorrect end tag label");
      }
   }

}

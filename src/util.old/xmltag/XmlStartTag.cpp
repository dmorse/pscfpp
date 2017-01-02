/*
* Simpatico - Simulation Package for Polymeric and Molecular Liquids
*
* Copyright 2010 - 2014, The Regents of the University of Minnesota
* Distributed under the terms of the GNU General Public License.
*/

#include "XmlStartTag.h"
#include <util/misc/ioUtil.h>

namespace Util
{

   XmlStartTag::XmlStartTag()
    : endBracket_(false)
   {}

   XmlStartTag::~XmlStartTag()
   {}

   bool XmlStartTag::matchLabel(const std::string& line, int begin)
   {
      label_ = "";
      endBracket_ = false;
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

      // Read label
      int beginLabel, endLabel;
      next();
      if (isEnd()) return false;
      if (c() == '/') return false; // Indicates an end tag
      skip();
      if (isEnd()) return false;
      beginLabel = cursor();
      while (c() != ' ' && c() != '>') {
         next();
         if (isEnd()) return false;
      }
      endLabel = cursor();
      label_ = string().substr(beginLabel, endLabel - beginLabel);
      if (c() == '>') {
         endBracket_ = true;
      }
      return true;
   }
         
   void XmlStartTag::matchLabel(const std::string expected,
                                const std::string& line, int begin)
   {
      if (!matchLabel(line, begin)) {
         UTIL_THROW("Unable to match start tag label");
      }
      if (label() != expected) {
          Log::file() << "found    = " << label() << std::endl;
          Log::file() << "expected = " << expected << std::endl;
          UTIL_THROW("Incorrect start tag label");
      }
   }

   bool XmlStartTag::matchAttribute(XmlAttribute& attribute)
   {
      skip();
      if (isEnd()) return false;
      if (c() == '>') {
         endBracket_ = true;
         return false;
      } else if (c() == '/') {
         next();
         if (isEnd()) return false;
         if (c() == '>') {
            endBracket_ = true;
         }
         return false;
      } else {
         return attribute.match(*this);
      }
   }

   /*
   * Throw exception if now end bracket was found.
   */
   void XmlStartTag::finish()
   {
      if (!endBracket()) {
         Log::file() << "line = " << string() << std::endl;
         UTIL_THROW("Missing end bracket");
      }
   }

}

/*
* Simpatico - Simulation Package for Polymeric and Molecular Liquids
*
* Copyright 2010 - 2014, The Regents of the University of Minnesota
* Distributed under the terms of the GNU General Public License.
*/

#include "XmlXmlTag.h"
#include <util/misc/ioUtil.h>

namespace Util
{

   XmlXmlTag::XmlXmlTag()
    : endBracket_(false)
   {}

   XmlXmlTag::~XmlXmlTag()
   {}

   bool XmlXmlTag::match(const std::string& line, int begin)
   {
      // Check label
      if (!matchLabel(line, begin)) return false;

      // Find attributes
      XmlAttribute attribute;
      bool version = false;
      bool encoding = false;
      while (matchAttribute(attribute)) {
         if (attribute.label() != "version") {
            if (version) return false;
            version = true;
         } else
         if (attribute.label() != "encoding") {
            if (encoding) return false;
            encoding = true;
         } else {
            return false;
         }
      }
      if (!version) return false;
      if (!encoding) return false;

      // Check for end bracket
      if (!endBracket()) {
         return false;
      }
      return true;
   }

   bool XmlXmlTag::matchLabel(const std::string& line, int begin)
   {
      label_ = "";
      endBracket_ = false;
      setString(line, begin);

      // Skip leading white space
      if (isEnd()) return false;
      skip();
      if (isEnd()) return false;

      // Read opening <? bracket
      if (c() != '<') {
         return false;
      }
      next();
      if (c() != '?') {
         return false;
      }

      // Read "xml" label
      int beginLabel, endLabel;
      next();
      skip();
      if (isEnd()) return false;
      beginLabel = cursor();
      while (c() != ' ') {
         next();
         if (isEnd()) return false;
      }
      endLabel = cursor();
      label_ = string().substr(beginLabel, endLabel - beginLabel);
      if (label_ != "xml") {
         return false;
      } else {
         return true;
      }

   }
         
   bool XmlXmlTag::matchAttribute(XmlAttribute& attribute)
   {
      skip();
      if (isEnd()) return false;
      if (c() == '?') {
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

}

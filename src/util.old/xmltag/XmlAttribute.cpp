/*
* Simpatico - Simulation Package for Polymeric and Molecular Liquids
*
* Copyright 2010 - 2014, The Regents of the University of Minnesota
* Distributed under the terms of the GNU General Public License.
*/

#include "XmlAttribute.h"
#include <util/misc/ioUtil.h>

namespace Util
{

   XmlAttribute::XmlAttribute()
    : XmlBase(),
      label_(""),
      value_()
   {
      value_.str("");
      value_.clear(); 
   }

   XmlAttribute::~XmlAttribute()
   {} 

   bool XmlAttribute::match(const std::string& line, int begin)
   {
      // Clear members 
      label_ = "";
      value_.str("");
      value_.clear(); 

      setString(line, begin);

      // Skip white space
      if (isEnd()) return false;
      skip();
      if (isEnd()) return false;

      // Read label
      int beginLabel, endLabel;
      if (c() == '=') return false;
      if (c() == '>') return false;
      beginLabel = cursor();
      while (c() != '=') {
         next();
         if (isEnd()) return false;
         if (c() == '>') return false;
      }
      endLabel = cursor();

      // Advance past '=' and white space
      next();
      skip();
      if (isEnd()) return false;

      // Read value
      int beginValue, endValue;
      char quote;
      if (c() == '\'' || c() == '\"') {
         quote = c();
      } else {
         return false;
      }
      next();
      if (isEnd()) return false;
      beginValue = cursor();
      while (c() != quote) {
         next();
         if (isEnd()) return false;
      }
      endValue = cursor();
      next();
  
      label_ = string().substr(beginLabel, endLabel - beginLabel); 
      rStrip(label_);
      value_.str(string().substr(beginValue, endValue - beginValue)); 
 
      return true;
   }

   bool XmlAttribute::match(XmlBase& parser)
   {
      bool status = match(parser.string(), parser.cursor());
      if (status) {
         parser.setCursor(cursor());
      }
      return status;
   }

}

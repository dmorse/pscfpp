#ifndef UTIL_XML_ATTRIBUTE_H
#define UTIL_XML_ATTRIBUTE_H

/*
* Simpatico - Simulation Package for Polymeric and Molecular Liquids
*
* Copyright 2010 - 2014, The Regents of the University of Minnesota
* Distributed under the terms of the GNU General Public License.
*/

#include "XmlBase.h"
#include <sstream>
#include <vector>

namespace Util
{
   
   /**
   * Parser for an XML attribute.
   * 
   * \ingroup XmlTag_Module
   */
   class XmlAttribute : public XmlBase
   {

   public:

      /**
      * Constructor
      */ 
      XmlAttribute();

      /**
      * Destructor
      */ 
      virtual ~XmlAttribute();

      /**
      * Return true if an attribute is found, false otherwise.
      */
      bool match(const std::string& string, int begin);

      /**
      * If successful return true and advance cursor or parent parser.
      *
      * \param parser parent parser object
      */
      bool match(XmlBase& parser);

      /**
      * Return label string.
      */
      const std::string& label()
      {  return label_; }

      /**
      * Return value string, without quotes.
      */
      std::stringstream& value()
      {  return value_; }

   private:

      std::string label_;
      std::stringstream value_;

   };

}
#endif

#ifndef UTIL_XML_XML_TAG_H
#define UTIL_XML_XML_TAG_H

/*
* Simpatico - Simulation Package for Polymeric and Molecular Liquids
*
* Copyright 2010 - 2014, The Regents of the University of Minnesota
* Distributed under the terms of the GNU General Public License.
*/

#include "XmlBase.h"
#include "XmlAttribute.h"
#include <sstream>
#include <vector>

namespace Util
{
   
   /**
   * Parser for an XML file declaration tag (first line in file).
   *
   * The match function attempts to match an xml file declaration tag, 
   * such as: <?xml version="1.0" encoding="UTF-8"?>. 
   * 
   * \ingroup XmlTag_Module
   */
   class XmlXmlTag : public XmlBase
   {

   public:

      /**
      * Constructor
      */ 
      XmlXmlTag();

      /**
      * Destructor
      */ 
      ~XmlXmlTag();

      /**
      * Attempt to match entire xml tag.
      *
      * \param string containing text of XML tag
      * \param begin  index of first character
      * \return true on match, false otherwise
      */
      bool match(const std::string& string, int begin);

   private:

      /**
      * Label string (name of XML element).
      */
      std::string label_;

      /**
      * Set true when end bracket found, false until then.
      */
      bool endBracket_;

      /**
      * Attempt to match opening bracket and xml label.
      *
      * \param string containing text of XML tag
      * \param begin  index of first character
      */
      bool matchLabel(const std::string& string, int begin);

      /**
      * Attempt to match an attribute.
      *
      * \param  attribute on return, matched attribute, if any
      * \return true if an attribute is found, false otherwise
      */
      bool matchAttribute(XmlAttribute& attribute);

      /**
      * True if a closing bracket was found.
      */
      bool endBracket()
      {  return endBracket_; }

   };

}
#endif

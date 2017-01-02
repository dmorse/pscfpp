#ifndef UTIL_XML_END_TAG_H
#define UTIL_XML_END_TAG_H

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
   * Parser for an XML end tag.
   *
   * \ingroup XmlTag_Module
   */
   class XmlEndTag : public XmlBase
   {

   public:

      /**
      * Constructor
      */ 
      XmlEndTag();

      /**
      * Destructor
      */ 
      ~XmlEndTag();

      /**
      * Attempt to match any end tag.
      * 
      * Return true if end tag found, false otherwise.
      *
      * \param string containing text of XML tag
      * \param begin  index of first character
      */
      bool match(const std::string& string, int begin);

      /**
      * Match a required end tag.
      * 
      * Throw exception is specified end tag does not match.
      *
      * \param expected expected label string
      * \param string containing text of XML tag
      * \param begin  index of first character
      */
      void match(const std::string expected,
                 const std::string& string, int begin);

      /**
      * Label string.
      */
      const std::string label()
      {  return label_; }

   private:

      /**
      * Label string (name of XML element).
      */
      std::string label_;

   };

}
#endif

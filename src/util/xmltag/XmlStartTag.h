#ifndef UTIL_XML_START_TAG_H
#define UTIL_XML_START_TAG_H

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
   * Parser for an XML start tag.
   * 
   * Usage:
   * \code
   *    XmlStartTag tag;
   *    XmlAttribute attribute;
   *    std::string line;
   *    tag.matchLabel(line, 0);
   *    while (matchAttribute(attribute)) {
   *       // process attribute;
   *    } 
   *    tag.finish();
   * \endcode
   *
   * \ingroup XmlTag_Module
   */
   class XmlStartTag : public XmlBase
   {

   public:

      /**
      * Constructor
      */ 
      XmlStartTag();

      /**
      * Destructor
      */ 
      ~XmlStartTag();

      /**
      * Match opening bracket and any label.
      *
      * \param string containing text of XML tag
      * \param begin  index of first character
      * \return true if match, false otherwise
      */
      bool matchLabel(const std::string& string, int begin);

      /**
      * Match opening bracket and a specific required label.
      *
      * Throws exception if no match.
      *
      * \param expected  expected label string
      * \param string  containing text of XML tag
      * \param begin index of first character
      */
      void matchLabel(const std::string expected, 
                      const std::string& string, int begin);

      /**
      * Attempt to match an attribute.
      *
      * \param  attribute on return, matched attribute, if any
      * \return true if an attribute is found, false otherwise
      */
      bool matchAttribute(XmlAttribute& attribute);

      /**
      * Check if end bracket was found.
      * 
      * Throws exception if no end bracket was found.
      */
      void finish();

      /**
      * Label string.
      */
      const std::string label()
      {  return label_; }

      /**
      * True if a closing bracket was found.
      */
      bool endBracket()
      {  return endBracket_; }

   private:

      /**
      * Label string (name of XML element).
      */
      std::string label_;

      /**
      * Set true when end bracket found, false until then.
      */
      bool endBracket_;

   };

}
#endif

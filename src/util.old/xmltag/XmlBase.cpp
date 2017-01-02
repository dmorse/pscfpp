/*
* Simpatico - Simulation Package for Polymeric and Molecular Liquids
*
* Copyright 2010 - 2014, The Regents of the University of Minnesota
* Distributed under the terms of the GNU General Public License.
*/

#include "XmlBase.h"

namespace Util
{

   /*
   * Constructor.
   */
   XmlBase::XmlBase() 
    : stringPtr_(0),
      end_(0),
      cursor_(0),
      c_('\0')
   {}

   /*
   * Destructor.
   */
   XmlBase::~XmlBase() 
   {}

   /*
   * Initialize string and cursor.
   */
   void XmlBase::setString(const std::string& string, int cursor)
   {
      stringPtr_ = &string;
      end_ = string.length();
      setCursor(cursor);
   }

   /*
   * Set cursor. String must already be set.
   */
   void XmlBase::setCursor(int cursor)
   {
      if (!stringPtr_) {
         UTIL_THROW("No string");
      }
      if (cursor > end_) {
         UTIL_THROW("Error: cursor > end_");
      }
      if (cursor < 0) {
         UTIL_THROW("Error: cursor < 0");
      }
      cursor_ = cursor;
      if (cursor < end_) {
         c_ = (*stringPtr_)[cursor];
      } else {
         c_ = '\0';
      }
   }

}

/*
* Simpatico - Simulation Package for Polymeric and Molecular Liquids
*
* Copyright 2010 - 2014, The Regents of the University of Minnesota
* Distributed under the terms of the GNU General Public License.
*/

#include "Bool.h"
#include "Format.h"

namespace Util
{

   // Default constructor.
   Bool::Bool()
    : value_(0),
      width_(Format::defaultWidth())
   {}

   // Constructor, value only.
   Bool::Bool(bool value)
    : value_(value),
      width_(Format::defaultWidth())
   {}

   /// Constructor, value and width.
   Bool::Bool(bool value, int width)
    : value_(value),
      width_(width)
   {}

   void Bool::setValue(bool value)
   {  value_ = value; }

   void Bool::setWidth(int width)
   { width_ = width; }

   bool Bool::value()
   { return value_; }

   int Bool::width()
   { return width_; }

   /*
   * Input stream extractor for an Bool object.
   */
   std::istream& operator>>(std::istream& in, Bool &object)
   {
      in >> object.value_;
      return in;
   }

   /*
   * Output stream inserter for an Bool object.
   */
   std::ostream& operator<<(std::ostream& out, const Bool &object)
   {
      out.width(object.width_);
      out << object.value_;
      return out;
   }

} 

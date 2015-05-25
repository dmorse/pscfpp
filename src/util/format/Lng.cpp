/*
* Simpatico - Simulation Package for Polymeric and Molecular Liquids
*
* Copyright 2010 - 2014, The Regents of the University of Minnesota
* Distributed under the terms of the GNU General Public License.
*/

#include "Lng.h"
#include "Format.h"

namespace Util
{

   // Default constructor.
   Lng::Lng()
    : value_(0),
      width_(20)
      //width_(Format::defaultWidth())
   {}

   // Constructor, value only.
   Lng::Lng(long int value)
    : value_(value),
      width_(20)
      //width_(Format::defaultWidth())
   {}

   /// Constructor, value and width.
   Lng::Lng(long int value, int width)
    : value_(value),
      width_(width)
   {}

   void Lng::setValue(long int value)
   {  value_ = value; }

   void Lng::setWidth(int width)
   {  width_ = width; }

   long int Lng::value()
   {  return value_; }

   int Lng::width()
   {  return width_; }

   /**
   * Input stream extractor for an Lng object.
   *
   * \param in      input stream
   * \param object  Lng object to be read from stream
   * \return modified input stream
   */
   std::istream& operator>>(std::istream& in, Lng &object)
   {
      in >> object.value_;
      return in;
   }

   /**
   * Output stream inserter for an Lng object.
   *
   * \param  out  output stream
   * \param  object   Lng to be written to stream
   * \return modified output stream
   */
   std::ostream& operator<<(std::ostream& out, const Lng &object)
   {
      out.width(object.width_);
      out << object.value_;
      return out;
   }

} 

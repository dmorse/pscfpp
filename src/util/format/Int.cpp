/*
* Simpatico - Simulation Package for Polymeric and Molecular Liquids
*
* Copyright 2010 - 2014, The Regents of the University of Minnesota
* Distributed under the terms of the GNU General Public License.
*/

#include "Int.h"
#include "Format.h"

#include <iostream>

namespace Util
{

   /* 
   * Default constructor.
   */
   Int::Int()
    : value_(0),
      width_(Format::defaultWidth())
   {}

   /* 
   * Constructor, value only.
   */
   Int::Int(int value)
    : value_(value),
      width_(Format::defaultWidth())
   {}

   /* 
   * Constructor, value and width.
   */
   Int::Int(int value, int width)
    : value_(value),
      width_(width)
   {}

   /*
   * Set value.
   */
   void Int::setValue(int value)
   {  value_ = value; }

   /*
   * Set field Width.
   */
   void Int::setWidth(int width)
   { width_ = width; }

   /*
   * Get value.
   */
   int Int::value()
   { return value_; }

   /*
   * Get field width.
   */
   int Int::width()
   { return width_; }

   /**
   * Input stream extractor for an Int object.
   *
   * \param in  input stream
   * \param object  Int object to be read from stream
   * \return modified input stream
   */
   std::istream& operator>>(std::istream& in, Int &object)
   {
      in >> object.value_;
      return in;
   }

   /**
   * Output stream inserter for an Int object.
   *
   * \param  out  output stream
   * \param  object   Int to be written to stream
   * \return modified output stream
   */
   std::ostream& operator<<(std::ostream& out, const Int &object)
   {
      out.width(object.width_);
      out << object.value_;
      return out;
   }
  
} 

#ifndef UTIL_INT_H
#define UTIL_INT_H

/*
* Simpatico - Simulation Package for Polymeric and Molecular Liquids
*
* Copyright 2010 - 2014, The Regents of the University of Minnesota
* Distributed under the terms of the GNU General Public License.
*/

#include <iostream>

namespace Util
{

   /**
   * Wrapper for an int, for formatted ostream output.
   *
   * An Int object has a int numerical value, and a minimum output 
   * field width. The << operator for an Int uses the specified width.
   * The numerical value and width may both be specified as parameters 
   * to a constructor. If the width is not specified as a constructor
   * parameter, it is set within the constructor to a default value
   * equal to Format::defaultWidth().
   *
   * An Int object may be passed to an ostream as a temporary object.
   * For example, the expression:
   * \code
   *    std::cout << Int(13) << Int(25, 10) << std::endl;
   * \endcode
   * outputs the number 13 using the default width, followed by the
   * number 25 in a field of minimum width 10.
   *
   * \ingroup Format_Module
   */
   class Int
   {

   public:

      /// \name Constructors
      //@{ 
      
      /// Default constructor.
      Int();

      /// Constructor, value only.
      explicit Int(int value);

      /// Constructor, value and width.
      Int(int value, int width);

      //@}
      ///\name Setters
      //@{

      /**
      * Set the integer value.
      *
      * \param value value of the associated int variable
      */      
      void setValue(int value);

      /**
      * Set the output field width.
      *
      * \param width output field width
      */      
      void setWidth(int width);

      //@}
      ///\name Accessors
      //@{
      
      /**
      * Get the integer value.
      */      
      int value();

      /**
      * Get the minimum field width.
      */      
      int width();

      //@}
      
   private:

      /// Value of associated integer.
      int value_;

      /// Minimum field width.
      int width_;

   //friends:

      friend std::istream& operator>>(std::istream& in, Int &object);
      friend std::ostream& operator<<(std::ostream& out, const Int &object);

   };

   /**
   * Input stream extractor for an Int object.
   *
   * \param in  input stream
   * \param object  Int object to be read from stream
   * \return modified input stream
   */
   std::istream& operator>>(std::istream& in, Int &object);

   /**
   * Output stream inserter for an Int object.
   *
   * \param  out  output stream
   * \param  object   Int to be written to stream
   * \return modified output stream
   */
   std::ostream& operator<<(std::ostream& out, const Int &object);

} 
#endif

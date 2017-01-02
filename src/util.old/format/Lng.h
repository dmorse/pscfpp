#ifndef UTIL_LNG_H
#define UTIL_LNG_H

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
   * Wrapper for a long int, for formatted ostream output.
   *
   * An Lng object has a long int numerical value, and a minimum output 
   * field width. The << operator for an Lng uses the specified width.
   * The numerical value and width may both be optionally specified as 
   * parameters to a constructor. If the width is not specified, it is 
   * is is set to a default value equal to Format::defaultWidth().
   *
   * An Lng object may be passed to an ostream as a temporary object.
   * For example, the expression:
   * \code
   *    std::cout << Lng(13) << Lng(25, 10) << std::endl;
   * \endcode
   * outputs the number 13 using the default width, followed by the
   * number 25 in a field of minimum width 10.
   *
   * \ingroup Format_Module
   */
   class Lng
   {

   public:

      /// \name Constructors
      //@{ 
      
      /**
      * Default constructor.
      */
      Lng();

      /**
      * Constructor, value only.
      *
      * \param value associated long int
      */
      explicit Lng(long int value);

      /**
      * Constructor, value and width.
      *
      * \param value associated long int
      * \param width field width
      */
      Lng(long int value, int width);

      //@}
      ///\name Setters
      //@{
      
      /**
      * Set value of long int.
      *
      * \param value associated long int
      */
      void setValue(long int value);

      /**
      * Set field width.
      *
      * \param width field width
      */
      void setWidth(int width);

      //@}
      ///\name Accessors
      //@{
      
      /**
      * Get value of long int.
      */
      long int value();

      /**
      * Get field width.
      */
      int width();

      //@}
      
   private:

      /// Associated long int.
      long int value_;

      /// Output field width.
      int      width_;

   //friends:

      friend std::istream& operator>>(std::istream& in, Lng &object);
      friend std::ostream& operator<<(std::ostream& out, const Lng &object);

   };

   /**
   * Input stream extractor for an Lng object.
   *
   * \param in  input stream
   * \param object  Lng object to be read from stream
   * \return modified input stream
   */
   std::istream& operator>>(std::istream& in, Lng &object);

   /**
   * Output stream inserter for an Lng object.
   *
   * \param  out  output stream
   * \param  object   Lng to be written to stream
   * \return modified output stream
   */
   std::ostream& operator<<(std::ostream& out, const Lng &object);

} 
#endif

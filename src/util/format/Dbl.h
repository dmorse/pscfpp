#ifndef UTIL_DBL_H
#define UTIL_DBL_H

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
   * Wrapper for a double precision number, for formatted ostream output.
   *
   * An Dbl object has a double precision numerical value, as well as
   * members (width and precision) that control its output format. The
   * << operator for an Dbl object uses the specified width and precision.
   * The double precision number, the width and the precision may all be
   * specified as parameters to one of several constructors. Values of
   * width and precision that are not specified as parameters of a
   * constructor are set within the constructor to default values given
   * by Format::defaultWidth() and Format::defaultPrecision(), respectively. 
   *
   * A Dbl object may be passed to an ostream as a temporary object.
   * For example, the expression:
   * \code
   *    std::cout << Dbl(2.0) << Dbl(3.0, 15, 8) << std::endl;
   * \endcode
   * outputs the number 2.0 using the default width and precision,
   * followed by the number 3.0 in a field of minimum width 15 and 
   * precision 8. 
   *
   * \ingroup Format_Module
   */
   class Dbl
   {

   public:

      /// \name Constructors
      //@{ 

      /**      
      * Default constructor.
      */
      Dbl();

      /**      
      * Constructor, value only.
      */
      explicit Dbl(double value);

      /**      
      * Constructor, value and width.
      */
      Dbl(double value, int width);

      /** 
      * Constructor: value, width, precision, and format.
      */
      Dbl(double value, int width, int precision, bool isFixed = false);

      //@}
      ///\name Mutators
      //@{
      
      /**      
      * Set value of associated double.
      */
      void setValue(double value);

      /**      
      * Set output field width.
      */
      void setWidth(int width);

      /**      
      * Set output floating point precision.
      */
      void setPrecision(int precision);

      //@}
      ///\name Accessors
      //@{
      
      /**      
      * Get value of associated double.
      */
      double value();

      /**      
      * Get field width.
      */
      int width();

      /**      
      * Get floating point precision.
      */
      int precision();
     
      //@}
      
   private:

      /// Value of number to be output.
      double value_;

      /// Output field width.
      int    width_;

      /// Output floating point precision.
      int    precision_;

      /// Is format fixed precision?
      bool   isFixed_;

   //friends:

      friend std::istream& operator>>(std::istream& in, Dbl &object);
      friend std::ostream& operator<<(std::ostream& out, const Dbl &object);

   };

   /**
   * Input stream extractor for an Dbl object.
   *
   * \param in      input stream
   * \param object  Dbl object to be read from stream
   * \return modified input stream
   */
   std::istream& operator>>(std::istream& in, Dbl &object);

   /**
   * Output stream inserter for an Dbl object.
   *
   * \param  out  output stream
   * \param  object   Dbl to be written to stream
   * \return modified output stream
   */
   std::ostream& operator<<(std::ostream& out, const Dbl &object);

} 
#endif

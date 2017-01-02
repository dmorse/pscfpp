#ifndef UTIL_FORMAT_H
#define UTIL_FORMAT_H

/*
* Simpatico - Simulation Package for Polymeric and Molecular Liquids
*
* Copyright 2010 - 2014, The Regents of the University of Minnesota
* Distributed under the terms of the GNU General Public License.
*/

//#include <iostream>

namespace Util
{

   /**
   * Base class for output wrappers for formatted C++ ostream output.
   *
   * Public members are all getters and setters for static variables
   * defaultWidth and defaultPrecision.
   *
   * \ingroup Format_Module
   */
   class Format
   {

   public:

      /// Initialize or reset default width and precision values.
      static void initStatic();

      /// Set the default output field width.
      static void setDefaultWidth(int width);

      /// Set the default output precision.
      static void setDefaultPrecision(int precision);

      /// Return the default output field width.
      static int defaultWidth();

      /// Return the default output precision.
      static int defaultPrecision();

   private:

      /// Default width of field for formatted output
      static int defaultWidth_;

      /// Default precision for formatted output of floating point numbers.
      static int defaultPrecision_;

   };

} 
#endif

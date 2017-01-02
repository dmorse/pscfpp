#ifndef UTIL_WRITE_H
#define UTIL_WRITE_H

/*
* Simpatico - Simulation Package for Polymeric and Molecular Liquids
*
* Copyright 2010 - 2014, The Regents of the University of Minnesota
* Distributed under the terms of the GNU General Public License.
*/

#include <iostream>
#include <complex>

namespace Util
{

   /**
   * Function template for output in a standard format.
   *
   * The write function template is designed to simplify
   * formatted stream output of variables within class and
   * function template, when the typename of a variable is
   * a template parameter.
   *
   * The primary template implementation simply invokes the 
   * insertion << operator for the specified type. For types
   * controlled by the primary template (i.e., those for which
   * no explicit specialization is provided) the expression
   * write(out, data) is equivalent to out << data.
   *
   * Explicit specializations of this method are provided for
   * int, long, double, bool, and string. Each of these uses
   * an appropriate wrapper class (Int, Lng, Dbl, Bool, or 
   * Str) to format output. For example, if data is an int,
   * write(out, data) is equivalent to out << Int(data). The
   * width and (if appropriate) precision are controlled by
   * Format::defaultWidth() and Format::defaultWidth().
   *
   * \ingroup Format_Module
   */
   template <typename Type>
   void write(std::ostream& out, Type data);

   /**
   * Explicit specialization of write for double data.
   */
   template <> void write(std::ostream& out, double data);
   
   /**
   * Explicit specialization of write for double data.
   */
   template <> void write(std::ostream& out, std::complex<double> data);
   
   /**
   * Explicit specialization of write for int data.
   */
   template <> void write(std::ostream& out, int data);
   
   /**
   * Explicit specialization of write for long data.
   */
   template <> void write(std::ostream& out, long data);
   
   /**
   * Explicit specialization of write for bool data.
   */
   template <> void write(std::ostream& out, bool data);
   
   /**
   * Explicit specialization of write for std::string data.
   */
   template <> void write(std::ostream& out, std::string data);

   // Template definition

   /*
   * Primary template.
   */
   template <typename Type>
   inline void write(std::ostream& out, Type data)
   { out << data; }

}
#endif

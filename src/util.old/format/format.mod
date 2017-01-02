namespace Util
{

   /**
   * \defgroup Format_Module Output Format
   * \ingroup Util_NS_Module
   *
   * \brief Utilities to simplify formatted C++ stream output.
   *
   * This module provides wrapper classes that can simplify formatted
   * output of the primitive data types with controllable field width
   * and floating point precision. 
   *
   * \section Wrapper Classes
   *
   * The classes Int, Lng, Dbl, Bool, and Str are wrappers for outputting
   * the data types int, long double, bool, and std::string, respectively.
   * An inserter (<<) operator is defined for each such wrapper class that
   * produces formatted output of the enclosed data with a controllable 
   * field width and (for Dbl) precision. Each wrapper class has a member 
   * variable of the associated data type and an integer field width member. 
   * The Dbl class also has an integer precision member, to control floating 
   * point precision.
   *
   * Example: We wish to output the elements of two double precision
   * precision array named "A" and "B" in two column with a minimum field
   * width of 20 characters for elements of A, with 10 digit precision,
   * and 10 characters for elements of B, with 6 digit precision. The
   * following code accomplishes this:
   * \code
   * double A[10], B[10];
   *
   * // ... code that assigns values to elements of A and B ...
   *
   * for (int i=0; i< 10; ++i) {
   *    std::cout << Dbl(A[i], 20, 10) << Dbl(B[i], 10, 6) << std::endl;
   * }
   *
   * \endcode
   * The Dbl constructor used in this snippet has the interface
   * Dbl::Dbl(double value, int width, int precision). The use of
   * wrapper classes allows one to control output format using an
   * an interface that is more compact than the C++ iostream
   * interace, and only slightly more verbose than that of the C
   * fprint function.
   *
   * Two or more constructors are provide for each wrapper class.
   * Each class has a constructor that requires only the value of
   * of the variable, while others require the value and field width
   * or (as in the above example) the value, width and precision.
   * If a field width or precision is not specified as a parameter
   * to the constructor, it may be set after construction using
   * setter functions.
   *
   * When no value is specified for the field width or (for Dbl) the
   * precision, default values are used. The default width and
   * precision for all data types are given by Format::defaultWidth()
   * and Format::defaultPrecision(). These default values may be
   * modified using the static methods Format::setDefaultWidth()
   * and Format::setDefaultPrecision().
   *
   * Example: Suppose we wish to output the two column array
   * described in the previous example, but are willing to use
   * a 15 column field an 7 digits of precision for both columns.
   * This could also be accomplished as follows:
   * \code
   * double A[10], B[10];
   *
   * Format::setDefaultWidth(15);
   * Format::setDefaultPrecision(7);
   *
   * for (int i=0; i< 10; ++i) {
   *    std::cout << Dbl(A[i]) << Dbl(B[i]) << std::endl;
   * }
   *
   * \endcode
   * The setDefaultWidth() and setDefaultPrecision() functions are
   * not needed if one is happy with the initial default settings,
   * which are a width of 20 characters and a precision of 12.
   *
   * \section Write Function Template
   *
   * The write() function template provides a generic interface for
   * formatting ostream output, which can be used within a class or
   * function template to output data for which the type is a template
   * parameter. The wrapper classes cannot be used directly in this
   * situation, because they require that an object of the appropriate
   * wrapper class be specified explicitly.  To output a variable data
   * to an ostream out, one calls write(out, data). An explicit
   * specialization of write() is provided for each data type for which
   * there exists a wrapper class. Each explicit specialization uses
   * the corresponding wrapper class internally to format the output.
   * Thus, if variable data is an int, write(out, data) is equivalent
   * to out << Int(data). For other data types, for which there exists
   * no wrapper class, write(out, data) is equivalent out << data.
   */

}

#ifndef PSCF_NAN_EXCEPTION_H
#define PSCF_NAN_EXCEPTION_H

/*
* PSCF - Polymer Self-Consistent Field Theory
*
* Copyright 2015 - 2025, The Regents of the University of Minnesota
* Distributed under the terms of the GNU General Public License.
*/

#include <util/global.h>                  
#include <util/misc/Exception.h>

namespace Pscf {

   using namespace Util;

   /**
   * Exception thrown when not-a-number (NaN) is encountered.
   * 
   * An exception to be thrown when a required numerical parameter has a 
   * value of NaN. This may happen, for instance, as a result of dividing
   * by zero. A NanException is used rather than a generic Exception in 
   * these instances so that the NanException can be caught by a try-catch
   * block.
   * 
   * A NanException behaves identically to a generic Exception, but with a
   * pre-defined error message rather than a user-specified message. There 
   * is no preprocessor macro to throw a NanException, so they must be thrown
   * using the constructor. This will typically assume the following syntax,
   * where values of the file and line parameters are given by the built-in 
   * macros __FILE__ and __LINE__, respectively:
   * \code
   *    throw Exception("MyClass::myFunction", __FILE__, __LINE__ );
   * \endcode
   *
   * \ingroup Pscf_Iterator_Module
   */
   class NanException : public Exception
   {

   public:

      /**
      * Constructor. 
      *
      * Constructs error message that includes file and line number. Values 
      * of the file and line parameters should be given by the built-in macros
      * __FILE__ and __LINE__, respectively, in the calling function. 
      *
      * \param function name of the function from which the Exception was thrown
      * \param file     name of the file from which the Exception was thrown
      * \param line     line number in file
      * \param echo     if echo, then echo to Log::file() when constructed. 
      */
      NanException(const char *function, const char *file, int line, 
                   int echo = 1);

      /**
      * Constructor without function name parameter.
      *
      * \param file     name of the file from which the Exception was thrown
      * \param line     line number in file
      * \param echo     if echo, then echo to std out when constructed. 
      */
      NanException(const char *file, int line, int echo = 1);

      /**
      * Destructor
      */   
      ~NanException();

   };

}
#endif
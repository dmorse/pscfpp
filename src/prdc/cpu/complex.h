#ifndef PRDC_CPU_COMPLEX_H
#define PRDC_CPU_COMPLEX_H

/*
* PSCF Package 
*
* Copyright 2016 - 2022, The Regents of the University of Minnesota
* Distributed under the terms of the GNU General Public License.
*/

#include <fftw3.h>
#include <iostream>

namespace Pscf {
namespace Prdc {
namespace Cpu {

   /**
   * Stream extraction operator for fftw_complex
   *
   * \param is  input stream
   * \param z   complex number
   */
   std::istream& operator >> (std::istream& is, fftw_complex & z);

   /**
   * Stream insertion operator for fftw_complex
   *
   * \param os  output stream
   * \param z  complex number
   * \return modified output stream
   */
   std::ostream& operator << (std::ostream& os, fftw_complex const & z);

   #if 0
   /**
   * Serialization function template for fftw_complex number.
   *
   * Implementation serializes real part, then imaginary part.
   *
   * \param ar Archive object
   * \param z complex number
   * \param version version id
   */
   template <typename Archive>
   void serialize(Archive& ar, double z[2], const unsigned int version = 0)
   {
      serialize(ar, z[0], version);
      serialize(ar, z[1], version);
   }
   #endif

}
}
}
#endif

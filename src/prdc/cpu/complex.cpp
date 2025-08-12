/*
* PSCF - Polymer Self-Consistent Field 
*
* Copyright 2015 - 2025, The Regents of the University of Minnesota
* Distributed under the terms of the GNU General Public License.
*/

#include "complex.h"

namespace Pscf {
namespace Prdc {
namespace Cpu {

   /*
   * Stream extraction operator for fftw_complex.
   */
   std::istream& operator >> (std::istream& is, fftw_complex& z)
   {
      is >> z[0] >> z[1];
      return is;
   }

   /*
   * Stream insertion operator for fftw_complex.
   */
   std::ostream& operator << (std::ostream& os, fftw_complex const & z)
   {
      os <<  z[0] << "  " << z[1];
      return os;
   }

}
}
}

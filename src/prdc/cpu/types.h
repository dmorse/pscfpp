#ifndef PRDC_CPU_TYPES_H
#define PRDC_CPU_TYPES_H

/*
* PSCF Package - Polymer Self-Consistent Field 
*
* Copyright 2015 - 2025, The Regents of the University of Minnesota
* Distributed under the terms of the GNU General Public License.
*/

#include <fftw3.h>

namespace Pscf {
namespace Prdc {
namespace Cpu {

   /**
   * Complex number type used in CPU code that uses FFTW.
   */
   typedef fftw_complex Complex;

   /**
   * Real number type used in CPU code that uses FFTW.
   */
   typedef double Real;

}
}
}
#endif

#ifndef PSCF_CPU_TYPES_H
#define PSCF_CPU_TYPES_H

/*
* PSCF Package - Polymer Self-Consistent Field 
*
* Copyright 2016 - 2023, The Regents of the University of Minnesota
* Distributed under the terms of the GNU General Public License.
*/

#include <fftw3.h>
#include <complex>

namespace Pscf {
namespace Cpu {

   typedef double        Real;
   typedef fftw_complex  Complex;
   //typedef std::complex<Real> Complex;

}
}
#endif

/*
* PSCF - Polymer Self-Consistent Field
*
* Copyright 2015 - 2026, The Regents of the University of Minnesota
* Distributed under the terms of the GNU General Public License.
*/

#include "HostRecvArray.h"

namespace Pscf {

   // Explicit instantiation definitions
   template class HostRecvArray<double,CPT>;
   template class HostRecvArray<fftw_complex,CPT>;

}

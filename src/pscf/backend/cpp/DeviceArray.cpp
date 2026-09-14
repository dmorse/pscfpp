/*
* PSCF - Polymer Self-Consistent Field
*
* Copyright 2015 - 2026, The Regents of the University of Minnesota
* Distributed under the terms of the GNU General Public License.
*/

#include "DeviceArray.h"

namespace Pscf {

   // Explicit instantiation definitions
   template class DeviceArray<double,CPT>; 
   template class DeviceArray<fftw_complex,CPT>; 

}

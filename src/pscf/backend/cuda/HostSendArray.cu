/*
* PSCF - Polymer Self-Consistent Field
*
* Copyright 2015 - 2026, David Morse (morse012@um.edu) and co-authors
* Distributed under the terms of the GNU General Public License.
*/

#include "HostSendArray.h"

namespace Pscf {

   // Explicit instantiation definitions
   template class HostSendArray<cudaReal,CUT>;
   template class HostSendArray<cudaComplex,CUT>;

} // namespace Pscf

/*
* PSCF - Polymer Self-Consistent Field Theory
*
* Copyright 2016 - 2022, The Regents of the University of Minnesota
* Distributed under the terms of the GNU General Public License.
*/

#include <prdc/cpu/Field.tpp>
#include <util/global.h>
#include <fftw3.h>

namespace Pscf { 
namespace Prdc { 
namespace Cpu { 

   template class Field<double>;
   template class Field<fftw_complex>;

}
}
}

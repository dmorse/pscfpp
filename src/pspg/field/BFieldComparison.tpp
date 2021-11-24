#ifndef PSPC_B_FIELD_COMPARISON_TPP
#define PSPC_B_FIELD_COMPARISON_TPP

/*
* PSCF - Polymer Self-Consistent Field Theory
*
* Copyright 2016 - 2019, The Regents of the University of Minnesota
* Distributed under the terms of the GNU General Public License.
*/

#include "BFieldComparison.h"

namespace Pscf {
namespace Pspg {

   template <int D>
   BFieldComparison <D>::BFieldComparison(int begin)
    : FieldComparison< RDField<D> > (begin)
   {}; 

}
}

#endif

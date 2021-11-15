/*
* PSCF - Polymer Self-Consistent Field Theory
*
* Copyright 2016 - 2019, The Regents of the University of Minnesota
* Distributed under the terms of the GNU General Public License.
*/

#include "BFieldComparison.h"

namespace Pscf {
namespace Pspc {

   // Constructor
   BFieldComparison::BFieldComparison(int begin)
    : FieldComparison< DArray<double> > (begin)
   {};

}
}

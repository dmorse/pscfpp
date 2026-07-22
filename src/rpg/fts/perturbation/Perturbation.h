#ifndef RPG_PERTURBATION_H
#define RPG_PERTURBATION_H

#include <rp/fts/perturbation/Perturbation.h>  // base class template
#include <pscf/backends/CUT.h>                  // base class argument

// Explicit instantiation declarations
namespace Pscf {
   namespace Rp {
      extern template class Perturbation<1,CUT>;
      extern template class Perturbation<2,CUT>;
      extern template class Perturbation<3,CUT>;
   }
}
#endif

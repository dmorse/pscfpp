#include "SpaceSymmetry.h"

#include <util/format/Int.h>

namespace Pscf
{

   using namespace Util;

   // Define static members
   template<> SpaceSymmetry<2> SpaceSymmetry<2>::identity_;
   template<> bool             SpaceSymmetry<2>::hasIdentity_ = false;

   template<> SpaceSymmetry<3> SpaceSymmetry<3>::identity_;
   template<> bool             SpaceSymmetry<3>::hasIdentity_ = false;

   // Explicit instantiation.
   template class SpaceSymmetry<2>;
   template class SpaceSymmetry<3>;
}

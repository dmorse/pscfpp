#ifndef PSPG_FIELD_STATE_TPP
#define PSPG_FIELD_STATE_TPP

/*
* PSCF - Polymer Self-Consistent Field Theory
*
* Copyright 2016 - 2022, The Regents of the University of Minnesota
* Distributed under the terms of the GNU General Public License.
*/

#include "FieldState.h"

#include <pspg/System.h>                   
#include <prdc/cuda/FFT.h>                
#include <prdc/crystal/Basis.h>           
#include <pscf/mesh/Mesh.h>                
#include <util/misc/FileMaster.h>       

// #include <util/format/Str.h>
// #include <util/format/Int.h>
// #include <util/format/Dbl.h>

namespace Pscf {
namespace Rpg
{

   using namespace Util;
   using namespace Pscf::Prdc;

   /*
   * Constructor.
   */
   template <int D, class FT>
   FieldState<D, FT>::FieldState()
    : fields_(),
      unitCell_(),
      systemPtr_(0)
   {}
  
   /*
   * Constructor.
   */
   template <int D, class FT>
   FieldState<D, FT>::FieldState(System<D>& system)
    : fields_(),
      unitCell_(),
      systemPtr_(0)
   { setSystem(system); }

   /*
   * Destructor.
   */
   template <int D, class FT>
   FieldState<D, FT>::~FieldState()
   {}

   /*
   * Set association with system, after default construction.
   */
   template <int D, class FT>
   void FieldState<D, FT>::setSystem(System<D>& system)
   {
      if (hasSystem()) {
         UTIL_CHECK(systemPtr_ = &system);
      } else {
         systemPtr_ = &system;
      }
   }

} // namespace Rpg
} // namespace Pscf
#endif

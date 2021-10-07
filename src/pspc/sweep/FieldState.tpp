#ifndef PSPC_FIELD_STATE_TPP
#define PSPC_FIELD_STATE_TPP

/*
* PSCF - Polymer Self-Consistent Field Theory
*
* Copyright 2016 - 2019, The Regents of the University of Minnesota
* Distributed under the terms of the GNU General Public License.
*/

#include "FieldState.h"

#include <pspc/System.h>                   
#include <pspc/field/FFT.h>                
#include <pscf/mesh/Mesh.h>                
#include <pscf/crystal/Basis.h>           
#include <util/misc/FileMaster.h>       

// #include <util/format/Str.h>
// #include <util/format/Int.h>
// #include <util/format/Dbl.h>

namespace Pscf {
namespace Pspc
{

   using namespace Util;

   /*
   * Constructor.
   */
   template <int D, class FT>
   FieldState<D, FT>::FieldState()
    : fields_(),
      unitCell_(),
      fieldIo_(),
      systemPtr_(0)
   {}
  
   /*
   * Constructor.
   */
   template <int D, class FT>
   FieldState<D, FT>::FieldState(System<D>& system)
    : fields_(),
      unitCell_(),
      fieldIo_(),
      systemPtr_(&system)
   {
      std::string groupName = systemPtr_->groupName();
      fieldIo_.associate(unitCell_, systemPtr_->mesh(), 
                         systemPtr_->fft(), groupName,
                         systemPtr_->basis(), systemPtr_->fileMaster());
   }

   /*
   * Destructor.
   */
   template <int D, class FT>
   FieldState<D, FT>::~FieldState()
   {}

} // namespace Pspc
} // namespace Pscf
#endif

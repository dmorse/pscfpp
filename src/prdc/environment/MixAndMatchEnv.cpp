/*
* PSCF - Polymer Self-Consistent Field Theory
*
* Copyright 2015 - 2025, The Regents of the University of Minnesota
* Distributed under the terms of the GNU General Public License.
*/

#include "MixAndMatchEnv.h"

namespace Pscf {
namespace Prdc {

   /*
   * Constructor.
   */
   MixAndMatchEnv::MixAndMatchEnv()
   {  setClassName("MixAndMatchEnv"); }

   /*
   * Destructor.
   */
   MixAndMatchEnv::~MixAndMatchEnv()
   {}

   /*
   * Modify stress to account for Environment, for one lattice parameter.
   */
   double MixAndMatchEnv::modifyStress(int paramId, double stress) const
   {
      UTIL_CHECK(!EnvironmentBase::needsUpdate());

      // Add stress contributions from each FieldGenerator
      if (fieldGenPtr1_) {
         stress += fieldGenPtr1_->stress(paramId);
      }
      if (fieldGenPtr2_) {
         stress += fieldGenPtr2_->stress(paramId);
      }

      // Allow each FieldGenerator to modify the resulting stress
      if (fieldGenPtr1_) {
         stress = fieldGenPtr1_->modifyStress(paramId, stress);
      }
      if (fieldGenPtr2_) {
         stress = fieldGenPtr2_->modifyStress(paramId, stress);
      }

      return stress;
   }

} // namespace Prdc
} // namespace Pscf
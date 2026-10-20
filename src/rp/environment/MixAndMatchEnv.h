#ifndef RP_MIX_AND_MATCH_ENV_H
#define RP_MIX_AND_MATCH_ENV_H

/*
* PSCF - Polymer Self-Consistent Field
*
* Copyright 2015 - 2026, The Regents of the University of Minnesota
* Distributed under the terms of the GNU General Public License.
*/

#include "Environment.h"
#include "FieldGenerator.h"
#include <pscf/environment/MixAndMatchEnvTmpl.h>

namespace Pscf {
namespace Rp {

   using namespace Util;

   /**
   * Environment to mix and match field generators in variable-cell SCFT.
   * 
   * The classes Environment, FieldGenerator, and MixAndMatchEnv in the
   * Rp namespace are extensions of EnvironmentBase, FieldGeneratorBase,
   * and MixAndMatchEnvTmpl in the Pscf namespace, respectively. These
   * Rp classes are designed to simply add a feature that is not defined
   * in the corresponding Pscf classes: the ability to calculate, store, 
   * and write an environment-modified stress that will be used to optimize 
   * the lattice parameters of the periodic unit cell. 
   * 
   * This class is therefore a subclass of MixAndMatchEnvTmpl that defines 
   * Environment as its Environment base class, and uses FieldGenerator as 
   * the base class for its field generators.  It defines only one method, 
   * modifyStress, and otherwise adopts all the behavior of its parent 
   * MixAndMatchEnvTmpl class.
   * 
   * \ingroup Rp_Environment_Module
   */
   class MixAndMatchEnv 
      : public MixAndMatchEnvTmpl<Environment, FieldGenerator>
   {

   public:

      /**
      * Constructor.
      */
      MixAndMatchEnv();

      /**
      * Destructor.
      */
      ~MixAndMatchEnv();

      /**
      * Modify stress to account for Environment, for one lattice parameter.
      * 
      * This method performs only the Environment-specific aspects of the
      * stress calculation, while all other aspects are handled by the
      * computeStress method of the parent Rp::Environment class. 
      * 
      * \param paramId  index of the lattice parameter with this stress
      * \param stress  baseline SCFT stress computed by a Mixture object
      */
      double modifyStress(int paramId, double stress) const;

   };

} // namespace Rp
} // namespace Pscf

#endif

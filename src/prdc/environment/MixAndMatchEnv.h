#ifndef PRDC_MIX_AND_MATCH_ENV_H
#define PRDC_MIX_AND_MATCH_ENV_H

/*
* PSCF - Polymer Self-Consistent Field
*
* Copyright 2015 - 2025, The Regents of the University of Minnesota
* Distributed under the terms of the GNU General Public License.
*/

#include "Environment.h"
#include "FieldGenerator.h"
#include <pscf/environment/MixAndMatchEnvTmpl.h>

namespace Pscf {
namespace Prdc {

   using namespace Util;

   /**
   * Environment that can mix and match field generators in variable-cell SCFT.
   * 
   * The classes Environment, FieldGenerator, and MixAndMatchEnv in the
   * Prdc namespace are extensions of EnvironmentBase, FieldGeneratorBase,
   * and MixAndMatchEnvTmpl in the Pscf namespace, respectively. These
   * Prdc classes are designed to simply add a feature that is not 
   * defined in the corresponding Pscf classes: the ability to calculate, 
   * store, and write an environment-modified stress that will be used to 
   * optimize the lattice parameters of the periodic unit cell. 
   * 
   * This class is therefore a subclass of MixAndMatchEnvTmpl that 
   * defines Prdc::Environment as its Environment base class, and uses 
   * Prdc::FieldGenerator as the base class for its field generators.
   * It defines only one method, modifyStress, and otherwise adopts all
   * the behavior of its parent MixAndMatchEnvTmpl class.
   * 
   * \ingroup Prdc_Environment_Module
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
      * computeStress method of the parent Prdc::Environment class. 
      * 
      * \param paramId  index of the lattice parameter with this stress
      * \param stress  baseline SCFT stress computed by a Mixture object
      */
      double modifyStress(int paramId, double stress) const;

   };

} // namespace Prdc
} // namespace Pscf

#endif

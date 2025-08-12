#ifndef PRDC_FIELD_GENERATOR_H
#define PRDC_FIELD_GENERATOR_H

/*
* PSCF - Polymer Self-Consistent Field
*
* Copyright 2015 - 2025, The Regents of the University of Minnesota
* Distributed under the terms of the GNU General Public License.
*/

#include <pscf/environment/FieldGeneratorBase.h> // base class

namespace Pscf {
namespace Prdc {

   using namespace Util;

   /**
   * Base class field generator for MixAndMatchEnv in variable-cell SCFT.
   * 
   * The classes Environment, FieldGenerator, and MixAndMatchEnv in the
   * Prdc namespace are extensions of EnvironmentBase, FieldGeneratorBase,
   * and MixAndMatchEnvTmpl in the Pscf namespace, respectively. These
   * Prdc classes are designed to simply add a feature that is not 
   * defined in the corresponding Pscf classes: the ability to calculate, 
   * store, and write an environment-modified stress that will be used to 
   * optimize the lattice parameters of the periodic unit cell. 
   * 
   * Therefore, this class is simply a subclass of FieldGeneratorBase that 
   * defines two additional methods, stress() and modifyStress(), which 
   * are used by Prdc::MixAndMatchEnv.
   * 
   * \ingroup Prdc_Environment_Module
   */
   class FieldGenerator : public FieldGeneratorBase
   {

   public:

      /**
      * Constructor
      */
      FieldGenerator();

      /**
      * Destructor
      */
      ~FieldGenerator();

      /**
      * Get contribution to the stress from this imposed field.
      * 
      * If the imposed field(s) change in a non-affine manner under changes
      * in the lattice parameters, then the "stress" used to optimize the
      * lattice parameters must contain additional terms arising from the
      * imposed field(s). Thus, this method may return a nonzero value 
      * regardless of the Type of this object. A return value of zero 
      * indicates that the imposed field(s) stretch affinely under a change 
      * in the given lattice parameter. 
      * 
      * The default implementation of the method raises an error. If a
      * FieldGenerator subclass is designed to be used only in a rigid 
      * box, then the subclass does not need to redefine this method,
      * as it will not ever be called. But subclasses must define the 
      * method if they are to be used with non-rigid boxes, or else the
      * error message in this default implementation will be triggered.
      * 
      * \param paramId  index of the lattice parameter being varied
      */
      virtual double stress(int paramId) const;

      /**
      * Modify stress to minimize a property other than fHelmholtz. 
      * 
      * It may be preferable with certain imposed fields to minimize a  
      * property other than fHelmholtz with respect to the lattice  
      * parameters. For instance, in a thin film it is useful to minimize 
      * the excess free energy per unit area, (fHelmholtz - fRef) * Delta, 
      * where fRef is a reference free energy and Delta is the film 
      * thickness. The information needed to perform such a modification
      * is often contained within the FieldGenerator objects. Therefore,
      * this method allows subclasses of FieldGenerator to modify the 
      * stress.
      * 
      * The method is called by the MixAndMatchEnv object that owns this
      * object, and the return value should be used to optimize the 
      * lattice parameters. 
      * 
      * By default, this method will simply return the value of stress
      * that is provided as an input, without performing a modification,
      * which corresponds to a stress that will be used to minimize 
      * fHelmholtz.
      * 
      * \param paramId  index of the lattice parameter with this stress
      * \param stress  stress value calculated by Mixture object
      */
      virtual double modifyStress(int paramId, double stress) const;

   };

} // namespace Prdc
} // namespace Pscf
#endif
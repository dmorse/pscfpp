#ifndef PSCF_PARAMETER_TYPE_H
#define PSCF_PARAMETER_TYPE_H

#include <string>

/*
* PSCF - Polymer Self-Consistent Field Theory
*
* Copyright 2016 - 2022, The Regents of the University of Minnesota
* Distributed under the terms of the GNU General Public License.
*/

namespace Pscf {

   class ParameterModifier; // Forward declaration, avoids circular reference

   /**
   * Declaration of a specialized sweep parameter type.
   * 
   * The SweepTmpl class has a member list named parameterTypes_
   * which contains a set of ParameterType objects that are added
   * at run time. Each of these ParameterType objects represents
   * a sweepable parameter that is associated with an object that
   * is a subclass of ParameterModifier. This latter object is 
   * responsible for defining the setter/getter functions that
   * actually affect the sweepable parameter.
   * 
   * \ingroup Pscf_Sweep_Module
   */
   struct ParameterType 
   {

      /**
      * Constructor.
      */
      ParameterType();

      /**
      * Alternate constructor that sets all members
      * 
      * \param name  String representing the name of this parameter
      * \param nId  The number of indices needed to specify this parameter
      * \param modifier  The ParameterModifier that owns this parameter
      */
      ParameterType(std::string name, int nId, ParameterModifier& modifier);

      /**
      * Destructor.
      */
      ~ParameterType();

      // String identifier
      std::string name;

      // Number of associated integer indices
      int nId;

      // Pointer to object that can get and set the parameter
      ParameterModifier* modifierPtr;

   };

} // namespace Pscf
#endif
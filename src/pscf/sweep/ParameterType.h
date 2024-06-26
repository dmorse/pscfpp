#ifndef PSCF_PARAMETER_TYPE_H
#define PSCF_PARAMETER_TYPE_H

#include "pscf/sweep/ParameterModifier.h" // base class

/*
* PSCF - Polymer Self-Consistent Field Theory
*
* Copyright 2016 - 2022, The Regents of the University of Minnesota
* Distributed under the terms of the GNU General Public License.
*/

namespace Pscf {

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
   */
   struct ParameterType 
   {

      /**
      * Constructor.
      */
      ParameterType();

      /**
      * Destructor.
      */
      ~ParameterType();

      // String identifier
      std::string name;

      // Number of associated integer indices
      int nId;

      // Pointer to object that can get and set the parameter
      ParameterModifier* modifierPtr_;

   };

} // namespace Pscf
#endif
#ifndef PSCF_PARAMETER_MODIFIER_H
#define PSCF_PARAMETER_MODIFIER_H

#include <string>
#include "util/containers/DArray.h"
#include "pscf/sweep/ParameterType.h" // base class

/*
* PSCF - Polymer Self-Consistent Field Theory
*
* Copyright 2016 - 2022, The Regents of the University of Minnesota
* Distributed under the terms of the GNU General Public License.
*/

namespace Pscf {

   using namespace Util;

   /**
   * Base class allowing subclasses to define sweepable parameters.
   * 
   * This class, along with SweepTmpl, define an interface which allows 
   * "specialized" sweepable parameters to be added to the Sweep class 
   * at run time. Any object with a specialized sweep parameter must 
   * be a subclass of ParameterModifier, and must redefine the three 
   * virtual methods getParameterTypes, setParameter, and getParameter.
   * The Sweep object is then responsible for calling the method 
   * getParameterTypes for any ParameterModifier in the system and 
   * adding the parameters to its list of specialized sweep parameters.
   * This should happen in the Sweep class constructor.
   * 
   * This interface was designed to allow objects whose class is 
   * determined at run time (i.e., an object created by a Factory) to
   * have sweepable parameters, which are added to the Sweep object
   * at run time.
   *
   * \ingroup Pscf_Sweep_Module
   */
   class ParameterModifier
   {

   public:

      /**
      * Constructor.
      */
      ParameterModifier();

      /**
      * Destructor.
      */
      ~ParameterModifier();

      /**
      * Return specialized sweep parameter types to add to the Sweep object.
      * 
      * This method should be called by the Sweep object in its constructor.
      * If the ParameterModifier object does not have any specialized sweep 
      * parameters, this method can be left as implemented here, returning 
      * an empty array.
      */
      virtual DArray<ParameterType> getParameterTypes()
      {  return DArray<ParameterType>(); } // empty array

      /**
      * Set the value of a specialized sweep parameter.
      * 
      * \param name  name of the specialized parameter
      * \param ids  array of integer indices specifying the value to set
      * \param value  the value to which the parameter is set
      */
      virtual 
      void setParameter(std::string name, DArray<int> ids, double value)
      {  UTIL_THROW("Error: Unimplemented setParameter function"); }

      /**
      * Get the value of a specialized sweep parameter.
      * 
      * \param name  name of the specialized parameter
      * \param ids  array of integer indices specifying the value to set
      */
      virtual 
      double getParameter(std::string name, DArray<int> ids) const
      {  
         UTIL_THROW("Error: Unimplemented getParameter function"); 
         return 0.0;
      }

   };

} // namespace Pscf
#endif
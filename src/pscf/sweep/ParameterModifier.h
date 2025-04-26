#ifndef PSCF_PARAMETER_MODIFIER_H
#define PSCF_PARAMETER_MODIFIER_H

#include "pscf/sweep/ParameterType.h"
#include "util/containers/GArray.h"
#include "util/containers/DArray.h"
#include <string>

/*
* PSCF - Polymer Self-Consistent Field Theory
*
* Copyright 2015 - 2025, The Regents of the University of Minnesota
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
      virtual GArray<ParameterType> getParameterTypes();

      /**
      * Set the value of a specialized sweep parameter.
      * 
      * Subclass implementations of setParameter should attempt to set 
      * the desired parameter, and should not throw an error regardless 
      * of whether this attempt succeeded. A boolean input parameter 
      * named success should be set to true or false depending on 
      * whether the parameter was successfully set. This method should
      * be used when the caller does not know if the specialized sweep
      * parameter belongs to this ParameterModifier or not. 
      * 
      * \param name  name of the specialized parameter
      * \param ids  array of integer indices specifying the value to set
      * \param value  the value to which the parameter is set
      * \param success  boolean flag used to indicate if parameter was set
      */
      virtual 
      void setParameter(std::string name, DArray<int> ids, 
                                          double value, bool& success)
      {  success = false; }

      /**
      * Get the value of a specialized sweep parameter.
      * 
      * Subclass implementations of getParameter should attempt to get 
      * the desired parameter, and should not throw an error regardless 
      * of whether this attempt succeeded. A boolean input parameter 
      * named success should be set to true or false depending on 
      * whether the parameter was successfully gotten. This method 
      * should be used when the caller does not know if the specialized 
      * sweep parameter belongs to this ParameterModifier or not. 
      * 
      * \param name  name of the specialized parameter
      * \param ids  array of integer indices specifying the value to set
      * \param success  boolean flag used to indicate if parameter was gotten
      */
      virtual 
      double getParameter(std::string name, DArray<int> ids, bool& success) 
      const
      {  
         success = false; 
         return 0.0;
      }

      /**
      * Set the value of a specialized sweep parameter.
      * 
      * This is an overloaded version of the setParameter method above,
      * which should be used only when the caller is certain that the
      * specialized sweep parameter belongs to this ParameterModifier.
      * An error will be thrown if the specialized parameter is not 
      * settable using this class.
      * 
      * \param name  name of the specialized parameter
      * \param ids  array of integer indices specifying the value to set
      * \param value  the value to which the parameter is set
      */
      void setParameter(std::string name, DArray<int> ids, double value);

      /**
      * Get the value of a specialized sweep parameter.
      * 
      * This is an overloaded version of the getParameter method above,
      * which should be used only when the caller is certain that the
      * specialized sweep parameter belongs to this ParameterModifier.
      * An error will be thrown if the specialized parameter is not 
      * gettable using this class.
      * 
      * \param name  name of the specialized parameter
      * \param ids  array of integer indices specifying the value to set
      */
      double getParameter(std::string name, DArray<int> ids) const;

   };

} // namespace Pscf
#endif
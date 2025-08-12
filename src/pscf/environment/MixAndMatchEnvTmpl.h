#ifndef PSCF_MIX_AND_MATCH_ENV_TMPL_H
#define PSCF_MIX_AND_MATCH_ENV_TMPL_H

/*
* PSCF - Polymer Self-Consistent Field
*
* Copyright 2015 - 2025, The Regents of the University of Minnesota
* Distributed under the terms of the GNU General Public License.
*/

#include <util/param/ParamComposite.h>    // base class of parent Env
#include <pscf/sweep/ParameterModifier.h> // base class of parent Env
#include <string>
#include <iostream>

namespace Pscf {

   using namespace Util;

   /**
   * Template class for Environments that mix and match field generators.
   * 
   * An Environment in PSCF is an object that generates a mask and/or a 
   * set of external fields to impose upon the calculation. In some cases,
   * one may wish to design multiple types of Environment subclasses that
   * generate the same mask but different external fields, or vice versa.
   * This class is designed to allow for such cases without the need for 
   * extensive duplicate code. In a MixAndMatchEnvTmpl object, the mask
   * generation and external field generation are handled by separate 
   * FieldGenerator objects, which are owned by the MixAndMatchEnvTmpl. 
   * This way, multiple subclasses of MixAndMatchEnvTmpl can use the same 
   * FieldGenerator for their masks or external fields, a "mix-and-match" 
   * approach.
   * 
   * This class has two template parameters, Env and FG, both of which
   * are class types. Env is the parent Environment class type, which 
   * may be EnvironmentBase or a subclass of EnvironmentBase. Similarly, 
   * FG is the class type for the FieldGenerator objects, which may be 
   * FieldGeneratorBase or a subclass of FieldGeneratorBase. This template
   * structure allows for this class to be paired with a namespace-specific 
   * Environment base class, rather than the more generic EnvironmentBase,
   * in cases where the Environment must assume namespace-specific 
   * functionality that is not defined in EnvironmentBase. In such cases, 
   * it may also be necessary to use a namespace-specific FieldGenerator 
   * base class that defines additional methods to accommodate the 
   * Environment class's namespace-specific functionality. 
   * 
   * Fully defined subclasses of MixAndMatchEnvTmpl should have names that 
   * correspond to the type of Environment they generate (e.g., 
   * "FilmEnvironment" for thin-film systems). The only method that needs 
   * to be implemented by subclasses is createGenerators, which should 
   * dynamically allocate the types of FieldGenerator objects needed by 
   * that Environment and store pointers to them. Otherwise, all behavior 
   * is defined in this parent class.
   * 
   * \ingroup Pscf_Environment_Module
   */
   template <class Env, class FG>
   class MixAndMatchEnvTmpl : public Env
   {

   public:

      /**
      * Constructor.
      */
      MixAndMatchEnvTmpl();

      /**
      * Destructor.
      */
      ~MixAndMatchEnvTmpl();
      
      /**
      * Read parameters from input stream.
      * 
      * \param in  input stream
      */
      void readParameters(std::istream& in);

      /**
      * Checks if fields need to be (re)generated. If so, generates them. 
      */
      void generate();

      /**
      * Get the first FieldGenerator by const reference.
      */
      FG const & fieldGenerator1() const;

      /**
      * Get the second FieldGenerator (if any) by const reference.
      */
      FG const & fieldGenerator2() const;

      /**
      * Does a second FieldGenerator exist?
      */
      bool hasFieldGenerator2() const;

      /**
      * Get specialized sweep parameter types to add to a Sweep object.
      */
      GArray<ParameterType> getParameterTypes();

      /**
      * Set the value of a specialized sweep parameter.
      * 
      * \param name  name of the specialized parameter
      * \param ids  array of integer indices specifying value to set
      * \param value  value to which the parameter is set
      * \param success  boolean flag used to indicate if parameter was set
      */
      void setParameter(std::string name, DArray<int> ids, double value, 
                        bool& success);
      
      /**
      * Get the value of a specialized sweep parameter.
      * 
      * \param name  name of the specialized parameter
      * \param ids  array of integer indices specifying the value to get
      * \param success  boolean  Was the parameter found?
      */
      double getParameter(std::string name, DArray<int> ids, bool& success)
      const;

      using Env::reset;
      using Env::needsUpdate;
      using ParameterModifier::setParameter; // overloaded method
      using ParameterModifier::getParameter; // overloaded method
      using ParamComposite::setClassName;
      using ParamComposite::addParamComposite;

   protected:

      using Env::setNeedsUpdateFalse;
      using Env::setGenerateBools;

      /**
      * Create FieldGenerator objects for mask and/or external field.
      * 
      * This method must be defined by each subclass. In it, fieldGenPtr1_
      * and (optionally) fieldGenPtr2_ should be assigned to dynamically 
      * allocated FieldGenerator objects (created with "new" command). 
      * The actual class of each of these objects will be a subclass of 
      * FieldGenerator, the type of which will vary in each subclass of
      * MixAndMatchEnvTmpl.
      * 
      * If only fieldGenPtr1_ is assigned, then it can be either type of 
      * FieldGenerator (Mask or External). If both fieldGenPtr1_
      * and fieldGenPtr2_ are assigned, then one of them must be of type
      * Mask and the other must be of type External.
      * 
      * This is the only method that must be defined by subclasses of this
      * class. All other methods are fully defined in this parent class.
      */
      virtual void createGenerators() = 0;

      /**
      * Pointer to the first FieldGenerator object (required).
      * 
      * This FieldGenerator may generate a mask or a set of external fields.
      */
      FG* fieldGenPtr1_;

      /**
      * Pointer to the second FieldGenerator object (optional).
      * 
      * If fieldGenPtr1_ points to an object that generates a mask, then
      * fieldGenPtr2_ points to an object that generates a set of external
      * fields, and vice versa. 
      */
      FG* fieldGenPtr2_;

   };

}
#include "MixAndMatchEnvTmpl.tpp"
#endif

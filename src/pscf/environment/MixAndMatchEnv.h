#ifndef PSCF_MIX_AND_MATCH_ENV_H
#define PSCF_MIX_AND_MATCH_ENV_H

/*
* PSCF - Polymer Self-Consistent Field Theory
*
* Copyright 2015 - 2025, The Regents of the University of Minnesota
* Distributed under the terms of the GNU General Public License.
*/

#include "FieldGenerator.h"
#include "Environment.h"
#include <string>

namespace Pscf {

   using namespace Util;

   /**
   * Environment that can mix and match mask and external field generators.
   * 
   * An Environment in PSCF is an object that generates a mask and/or a 
   * set of external fields to impose upon the calculation. In some cases,
   * one may wish to design multiple types of Environment subclasses that
   * generate the same mask but different external fields, or vice versa.
   * The MixAndMatchEnv class is designed to allow for such cases without
   * the need for extensive duplicate code. In a MixAndMatchEnv object, the 
   * mask generation and external field generation are handled by separate 
   * FieldGenerator objects, which are owned by the MixAndMatchEnv object. 
   * This way, multiple subclasses of MixAndMatchEnv can use the same 
   * FieldGenerator for their masks or external fields, a "mix-and-match" 
   * approach.
   * 
   * Subclasses of MixAndMatchEnv should have names that correspond to
   * the type of Environment they generate (e.g., "FilmEnvironment" for
   * thin-film systems). The only method that needs to be implemented by 
   * subclasses is createGenerators, which should dynamically allocate 
   * the types of FieldGenerator objects needed by that Environment and 
   * store pointers to them. Otherwise, all behavior is defined in this
   * parent class.
   * 
   * \ingroup Pscf_Environment_Module
   */
   class MixAndMatchEnv : public Environment
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
      * Return the Environment's contribution to the stress.
      * 
      * The value returned is the sum of the stress contributions from
      * the two FieldGenerator objects owned by this object. 
      * 
      * \param paramId  index of the lattice parameter with this stress
      */
      double stress(int paramId) const;

      /**
      * Modify stress to minimize a property other than fHelmholtz. 
      * 
      * The property that should be minimized will vary depending on the
      * type of Environment. In cases where fHelmholtz is the desired 
      * property to be minimized, this method returns the stress that
      * it was provided, without modification. 
      * 
      * By default, this method calls fieldGenerator1().modifyStress(), 
      * then passes the result into fieldGenerator2().modifyStress(), 
      * and returns the result. However, the method is virtual so that
      * a subclass may define an alternative approach to modifying the
      * stress if needed.
      * 
      * \param stress  unmodified stress value
      * \param paramId  index of the lattice parameter with this stress
      */
      virtual double modifyStress(int paramId, double stress) const;

      /**
      * Get the first FieldGenerator by const reference.
      */
      FieldGenerator const & fieldGenerator1() const;

      /**
      * Get the second FieldGenerator (if any) by const reference.
      */
      FieldGenerator const & fieldGenerator2() const;

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

      using Environment::reset;
      using ParameterModifier::setParameter; // overloaded method
      using ParameterModifier::getParameter; // overloaded method

   protected:

      using Environment::needsUpdate_;

      /**
      * Create FieldGenerator objects for mask and/or external field.
      * 
      * This method must be defined by each subclass. In it, fieldGenPtr1_
      * and (optionally) fieldGenPtr2_ should be assigned to dynamically 
      * allocated FieldGenerator objects (created with "new" command). 
      * The actual class of each of these objects will be a subclass of 
      * FieldGenerator, the type of which will vary in each subclass of
      * MixAndMatchEnv.
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
      FieldGenerator* fieldGenPtr1_;

      /**
      * Pointer to the second FieldGenerator object (optional).
      * 
      * If fieldGenPtr1_ points to an object that generates a mask, then
      * fieldGenPtr2_ points to an object that generates a set of external
      * fields, and vice versa. 
      */
      FieldGenerator* fieldGenPtr2_;

   };

}
#endif

#ifndef PSCF_IMPOSED_FIELDS_TMPL_H
#define PSCF_IMPOSED_FIELDS_TMPL_H

/*
* PSCF - Polymer Self-Consistent Field Theory
*
* Copyright 2015 - 2025, The Regents of the University of Minnesota
* Distributed under the terms of the GNU General Public License.
*/

#include <pscf/iterator/FieldGenerator.h>
#include <util/param/ParamComposite.h>
#include <util/param/Factory.h>
#include <string>

namespace Pscf {

   using namespace Util;

   /**
   * Base class defining mask & external fields to impose on the calculation.
   * 
   * Subclasses, which should be named ImposedFieldsGenerator in their 
   * respective namespaces, need to define only the createGenerators method, 
   * which uses the input parameter "type" to determine which types of
   * FieldGenerator objects to create, and creates them. The types of
   * FieldGenerator subclasses will, in general, depend on the namespace,
   * which changes which strings are accepted values of the "type" parameter
   * of this class.
   * 
   * \ingroup Pscf_Iterator_Module
   */
   class ImposedFieldsTmpl : public ParamComposite, public ParameterModifier
   {

   public:

      /**
      * Constructor.
      */
      ImposedFieldsTmpl();

      /**
      * Destructor.
      */
      ~ImposedFieldsTmpl();
      
      /**
      * Read parameters from input stream.
      * 
      * \param in  input stream
      */
      void readParameters(std::istream& in);

      /**
      * Allocate, check compatibility, calculate, and store the field(s).
      */
      void setup();

      /**
      * Check whether system has changed, update the field(s) if necessary.
      */
      void update();

      /**
      * Correct the stress value if necessary.
      * 
      * There are two changes to the stress that may be necessary due to 
      * the presence of imposed fields. First, If the imposed fields 
      * change in a non-affine manner under changes in the lattice 
      * parameters, then the stress used to optimize the lattice 
      * parameters must contain additional contributions arising from the
      * imposed field(s). 
      * 
      * And second, it may be preferable with certain imposed fields to 
      * minimize a property other than fHelmholtz with respect to the 
      * lattice parameters. For instance, in a thin film it is useful to 
      * minimize the excess free energy per unit area, 
      * (fHelmholtz - fRef) * Delta, where fRef is a reference free 
      * energy and Delta is the film thickness. Therefore, a subclass 
      * of this class may modify the stress value accordingly via its 
      * modifyStress method.
      * 
      * This method first calls the stressTerm() methods of both
      * FieldGenerator objects and adds them to the original value of 
      * stress that was passed in as a parameter, correcting for the 
      * former effect described above. It then calls the 
      * modifyStress method of this object, which can correct for the
      * latter effect. 
      * 
      * The method should be called by Iterator classes that own this 
      * object and the return value should be used to compute error and 
      * optimize the lattice parameters. 
      * 
      * \param paramId  index of the lattice parameter with this stress
      * \param stress  stress value calculated by Mixture object
      */
      double correctedStress(int paramId, double stress) const;

      /**
      * Get the type string associated with this object.
      */
      std::string type() const;

      /**
       * Return const references to the FieldGenerator child objects
       */
      FieldGenerator const & fieldGenerator1() const;
      FieldGenerator const & fieldGenerator2() const;

      /**
      * Return specialized sweep parameter types to add to a Sweep object.
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
      * \param success  boolean flag used to indicate if parameter was gotten
      */
      double getParameter(std::string name, DArray<int> ids, bool& success)
      const;

      using ParameterModifier::setParameter; // overloaded method
      using ParameterModifier::getParameter; // overloaded method

   protected:

      /**
      * Create FieldGenerator objects for mask and/or external field.
      * 
      * This method must be defined in each subclass. In it, fieldGenPtr1_
      * and (optionally) fieldGenPtr2_ should be assigned to dynamically 
      * allocated FieldGenerator objects (created with "new" command). 
      * The actual class of each of these objects will be a subclass of 
      * FieldGenerator, and the type member of the FieldGenerator objects 
      * will depend on the "type" parameter that is read by this object. 
      * 
      * If only fieldGenPtr1_ is assigned, then it can be any type of 
      * FieldGenerator (Mask, External, or Both). If both fieldGenPtr1_
      * and fieldGenPtr2_ are assigned, then one of them must be of type
      * Mask and the other must be of type External.
      * 
      * This is the only method that must be defined by subclasses of this
      * class. All other methods are fully defined in this class.
      */
      virtual void createGenerators() = 0;

      /**
      * Modify the stress value if necessary.
      * 
      * It may be preferable with certain imposed fields to minimize a 
      * property other than fHelmholtz with respect to the lattice 
      * parameters. For instance, in a thin film it is useful to 
      * minimize the excess free energy per unit area, 
      * (fHelmholtz - fRef) * Delta, where fRef is a reference free 
      * energy and Delta is the film thickness. Therefore, a subclass 
      * of this class may modify the stress value accordingly via this
      * method, which is called by the method correctedStress after
      * adding in the stress contributions from non-affine distortions
      * of the imposed fields.
      * 
      * By default, this method will simply return the value of the 
      * stress that it was passed, without applying any modifications.
      * 
      * \param paramId  index of the lattice parameter with this stress
      * \param stress  stress value
      */
      virtual double modifyStress(int paramId, double stress) const;

      /**
      * Pointer to the first FieldGenerator object (required).
      * 
      * This FieldGenerator may generate a mask, a set of external fields,
      * or both.
      */
      FieldGenerator* fieldGenPtr1_;

      /**
      * Pointer to the second FieldGenerator object (optional).
      * 
      * If fieldGenPtr1_ points to an object that generates a mask, then
      * fieldGenPtr2_ points to an object that generates a set of external
      * fields, and vice versa. If fieldGenPtr1_ points to an object that 
      * generates both, then fieldGenPtr2_ must be a null pointer.
      */
      FieldGenerator* fieldGenPtr2_;

   private:

      /// String that defines the type of fields to impose.
      std::string type_;

   };
}
#endif

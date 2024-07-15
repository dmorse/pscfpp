#ifndef PSCF_IMPOSED_FIELDS_TMPL_H
#define PSCF_IMPOSED_FIELDS_TMPL_H

/*
* PSCF - Polymer Self-Consistent Field Theory
*
* Copyright 2016 - 2022, The Regents of the University of Minnesota
* Distributed under the terms of the GNU General Public License.
*/

#include <util/param/ParamComposite.h>
#include <util/param/Factory.h>
#include <pscf/iterator/FieldGenerator.h>
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
      * Constructor
      */
      ImposedFieldsTmpl();

      /**
      * Destructor
      */
      ~ImposedFieldsTmpl();
      
      /**
      * Read parameters from input stream
      * 
      * \param in  input stream
      */
      void readParameters(std::istream& in);

      /**
      * Allocate, check compatibility, calculate, and store the field(s)
      */
      void setup();

      /**
      * Return specialized sweep parameter types to add to the Sweep object
      */
      DArray<ParameterType> getParameterTypes();

      /**
      * Set the value of a specialized sweep parameter
      * 
      * \param name  name of the specialized parameter
      * \param ids  array of integer indices specifying the value to set
      * \param value  the value to which the parameter is set
      * \param success  boolean flag used to indicate if parameter was set
      */
      void setParameter(std::string name, DArray<int> ids, 
                                          double value, bool& success);
      
      /**
      * Get the value of a specialized sweep parameter
      * 
      * \param name  name of the specialized parameter
      * \param ids  array of integer indices specifying the value to get
      * \param success  boolean flag used to indicate if parameter was gotten
      */
      double getParameter(std::string name, DArray<int> ids, bool& success)
      const;

   protected:

      /**
      * Create FieldGenerator objects for the mask & external field
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
      *  and fieldGenPtr2_ are assigned, then one of them must be of type
      * Mask and the other must be of type External.
      * 
      * This is the only method that must be defined by subclasses of this
      * class. All other methods are fully defined in this class.
      */
      virtual void createGenerators() = 0;

      /**
      * Get the type string associated with this object
      */
      std::string type() const;

      /**
      * Pointer to the first FieldGenerator object (required).
      * 
      * This FieldGenerator may generate a mask, a set of external
      * fields, or both.
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

      /// String that defines the type of fields to impose
      std::string type_;

   };
}
#endif
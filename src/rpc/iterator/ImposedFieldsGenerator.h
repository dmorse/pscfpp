#ifndef RPC_IMPOSED_FIELDS_GENERATOR_H
#define RPC_IMPOSED_FIELDS_GENERATOR_H

/*
* PSCF - Polymer Self-Consistent Field Theory
*
* Copyright 2016 - 2022, The Regents of the University of Minnesota
* Distributed under the terms of the GNU General Public License.
*/

#include <rpc/iterator/Iterator.h>
#include <pscf/iterator/ImposedFieldsTmpl.h>

namespace Pscf {
namespace Rpc {

   using namespace Util;

   template <int D> 
   class System;

   /**
   * Class defining mask & external fields to impose on the calculation.
   * 
   * \ingroup Rpc_Iterator_Module
   */
   template <int D>
   class ImposedFieldsGenerator : public ImposedFieldsTmpl
   {

   public:

      /**
      * Constructor
      * 
      * \param sys  System parent object
      */
      ImposedFieldsGenerator(System<D>& sys);

      /**
      * Destructor
      */
      ~ImposedFieldsGenerator();

      using ParamComposite::isActive;

   private:

      /**
      * Create FieldGenerator objects for the mask & external field
      * 
      * This method dynamically allocates FieldGenerator objects and
      * assigns fieldGenPtr1_ and fieldGenPtr2_ to them, where the 
      * actual type of each of these objects will be a subclass of 
      * FieldGenerator, and the type will depend on the type_ parameter
      * that is read by this object.
      */
      void createGenerators();

      /// Pointer to the associated system object.
      System<D>* sysPtr_;

      using ImposedFieldsTmpl::fieldGenPtr1_;
      using ImposedFieldsTmpl::fieldGenPtr2_;
      using ImposedFieldsTmpl::type;
      using ParamComposite::setClassName;

   };
}
}
#endif

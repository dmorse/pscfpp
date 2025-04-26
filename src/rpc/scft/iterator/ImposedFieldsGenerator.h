#ifndef RPC_IMPOSED_FIELDS_GENERATOR_H
#define RPC_IMPOSED_FIELDS_GENERATOR_H

/*
* PSCF - Polymer Self-Consistent Field Theory
*
* Copyright 2015 - 2025, The Regents of the University of Minnesota
* Distributed under the terms of the GNU General Public License.
*/

#include <rpc/scft/iterator/Iterator.h>
#include <pscf/iterator/ImposedFieldsTmpl.h>

namespace Pscf {
namespace Rpc {

   using namespace Util;

   template <int D> 
   class System;

   /**
   * Class defining mask & external fields to impose on the calculation.
   * 
   * \ingroup Rpc_Scft_Iterator_Module
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

      using ImposedFieldsTmpl::type;
      using ParamComposite::isActive;

   protected:

      /**
      * Modify the stress value if necessary.
      * 
      * It may be preferable with certain imposed fields to minimize a 
      * property other than fHelmholtz with respect to the lattice 
      * parameters. For instance, in a thin film it is useful to 
      * minimize the excess free energy per unit area, 
      * (fHelmholtz - fRef) * Delta, where fRef is a reference free 
      * energy and Delta is the film thickness. Therefore, this method
      * modifies the stress value accordingly, depending on the "type"
      * parameter that is read from the parameter file. This method is 
      * called by the method correctedStress after adding in the stress 
      * contributions from non-affine distortions of the imposed fields.
      * 
      * In some cases, the modification may be performed directly by 
      * this class, and in others it requires parameters that are stored
      * in the FieldGenerator objects. In the latter case, the method
      * FieldGenerator::modifyStress is called to perform the actual
      * modification. The specific way that the modification is 
      * performed is determined based on the "type" parameter of this 
      * object.
      * 
      * \param paramId  index of the lattice parameter with this stress
      * \param stress  stress value
      */
      virtual double modifyStress(int paramId, double stress) const;

      using ImposedFieldsTmpl::fieldGenPtr1_;
      using ImposedFieldsTmpl::fieldGenPtr2_;
      using ParamComposite::setClassName;

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

   };

   #ifndef RPC_IMPOSED_FIELDS_GENERATOR_TPP
   extern template class ImposedFieldsGenerator<1>;
   extern template class ImposedFieldsGenerator<2>;
   extern template class ImposedFieldsGenerator<3>;
   #endif

} // namespace Rpc
} // namespace Pscf
#endif

#ifndef RPG_IMPOSED_FIELDS_GENERATOR_TPP
#define RPG_IMPOSED_FIELDS_GENERATOR_TPP

/*
* PSCF - Polymer Self-Consistent Field Theory
*
* Copyright 2015 - 2025, The Regents of the University of Minnesota
* Distributed under the terms of the GNU General Public License.
*/

#include "ImposedFieldsGenerator.h"
#include <rpg/scft/iterator/MaskGenFilm.h>
#include <rpg/scft/iterator/ExtGenFilm.h>

namespace Pscf {
namespace Rpg {

   // Constructor
   template <int D>
   ImposedFieldsGenerator<D>::ImposedFieldsGenerator(System<D>& sys)
    : ImposedFieldsTmpl::ImposedFieldsTmpl(),
      sysPtr_(&sys)
   {  setClassName("ImposedFieldsGenerator"); }

   // Destructor
   template <int D>
   ImposedFieldsGenerator<D>::~ImposedFieldsGenerator()
   {
      if (fieldGenPtr1_) {
         delete fieldGenPtr1_;
      }
      if (fieldGenPtr2_) {
         delete fieldGenPtr2_;
      }
   }

   // Modify the stress value if necessary.
   template <int D>
   double ImposedFieldsGenerator<D>::modifyStress(int paramId, double stress) 
   const
   {
      if (type() == "film") {
         return fieldGenPtr1_->modifyStress(paramId, stress);
      } else {
         return stress;
      }
   }

   // Create FieldGenerator objects for the mask & external field
   template <int D>
   void ImposedFieldsGenerator<D>::createGenerators()
   {
      if (type() == "film") {
         fieldGenPtr1_ = new MaskGenFilm<D>(*sysPtr_);
         fieldGenPtr2_ = new ExtGenFilm<D>(*sysPtr_);
      } else {
         UTIL_THROW(("Unrecognized type parameter: " + type()).c_str());
      }
   }
   
}
}
#endif
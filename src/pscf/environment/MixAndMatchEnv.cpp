/*
* PSCF - Polymer Self-Consistent Field Theory
*
* Copyright 2015 - 2025, The Regents of the University of Minnesota
* Distributed under the terms of the GNU General Public License.
*/

#include "MixAndMatchEnv.h"
#include <util/param/Label.h>

namespace Pscf {

   /*
   * Constructor.
   */
   MixAndMatchEnv::MixAndMatchEnv()
    : fieldGenPtr1_(0),
      fieldGenPtr2_(0)
   {  setClassName("MixAndMatchEnv"); }

   /*
   * Destructor.
   */
   MixAndMatchEnv::~MixAndMatchEnv()
   {
      if (fieldGenPtr1_) {
         delete fieldGenPtr1_;
      }
      if (fieldGenPtr2_) {
         delete fieldGenPtr2_;
      }
   }

   /*
   * Read parameters from input stream.
   *
   * Parameters for field generators are read at the same indentation level
   * as parameters of the parent MaxAndMatchEnv object, without enclosed curly
   * bracket delimiters around blocks for each generator.
   *
   * If the second generator is "dependent", it backs up and re-reads the
   * block of the first generator.
   */
   void MixAndMatchEnv::readParameters(std::istream& in)
   {
      // Before reading parameters, create FieldGenerator objects
      createGenerators();

      // Save current istream position (allows for possible rewind to here)
      std::streampos pos = in.tellg();

      // Read first FieldGenerator
      if (fieldGenPtr1_) {

         UTIL_CHECK(!fieldGenPtr1_->isDependent());

         // Make fieldGenPtr1_ a child ParamComponent of this object
         bool next = false; // Set indentation to be same as that of parent
         addParamComposite(*fieldGenPtr1_, next);

         // Read parameters for this FieldGenerator.
         // Reads body without indentation or curly bracket delimiters
         fieldGenPtr1_->readParameters(in);

      } else {

         UTIL_THROW("Object must contain at least one FieldGenerator.");

      }

      // Read second FieldGenerator (optional)
      if (fieldGenPtr2_) {

         // Check that one FieldGenerator is a Mask and other is External
         if (fieldGenPtr2_->type() == FieldGenerator::External) {
            UTIL_CHECK(fieldGenPtr1_->type() == FieldGenerator::Mask);
         } else if (fieldGenPtr2_->type() == FieldGenerator::Mask) {
            UTIL_CHECK(fieldGenPtr1_->type() ==
                                            FieldGenerator::External);
         } else {
            UTIL_THROW("fieldGenPtr2_ must have type Mask or External.");
         }

         // If this FieldGenerator is dependent, rewind the istream
         if (fieldGenPtr2_->isDependent()) {
            in.seekg(pos);
            Label::clear();
         }

         // Make this FieldGenerator a child ParamComponent of this object
         bool next = false; // Set indentation to be same as that of parent
         addParamComposite(*fieldGenPtr2_, next);

         // Read parameters for this FieldGenerator.
         // Reads body without indentation or curly bracket delimiters
         fieldGenPtr2_->readParameters(in);
      }
   }

   /*
   * Check if fields need to be (re)generated. If so, generates them.
   */
   void MixAndMatchEnv::generate()
   {
      if (!needsUpdate_) return;
      if (fieldGenPtr1_) fieldGenPtr1_->generate();
      if (fieldGenPtr2_) fieldGenPtr2_->generate();
      needsUpdate_ = false;
   }

   /*
   * Return the Environment's contribution to the stress.
   */
   double MixAndMatchEnv::stress(int paramId) const
   {
      UTIL_CHECK(!needsUpdate());

      double stress(0.0);
      if (fieldGenPtr1_) {
         stress += fieldGenPtr1_->stress(paramId);
      }
      if (fieldGenPtr2_) {
         stress += fieldGenPtr2_->stress(paramId);
      }
      return stress;
   }

   /*
   * Modify stress to minimize a property other than fHelmholtz.
   */
   double MixAndMatchEnv::modifyStress(int paramId, double stress) const
   {
      UTIL_CHECK(!needsUpdate());

      if (fieldGenPtr1_) {
         stress = fieldGenPtr1_->modifyStress(paramId, stress);
      }
      if (fieldGenPtr2_) {
         stress = fieldGenPtr2_->modifyStress(paramId, stress);
      }
      return stress;
   }

   /*
   * Return specialized sweep parameter types to add to the Sweep object.
   */
   GArray<ParameterType> MixAndMatchEnv::getParameterTypes()
   {
      GArray<ParameterType> a1, a2;

      if (fieldGenPtr1_) a1 = fieldGenPtr1_->getParameterTypes();
      if (fieldGenPtr2_) a2 = fieldGenPtr2_->getParameterTypes();

      for (int i = 0; i < a2.size(); i++) {
         a1.append(a2[i]);
      }

      return a1;
   }

   /*
   * Set the value of a specialized sweep parameter.
   */
   void MixAndMatchEnv::setParameter(std::string name, DArray<int> ids,
                                     double value, bool& success)
   {
      success = false;
      if (fieldGenPtr1_) {
         fieldGenPtr1_->setParameter(name, ids, value, success);
      }
      if ((!success) && (fieldGenPtr2_)) {
         fieldGenPtr2_->setParameter(name, ids, value, success);
      }
      if (success) reset();
   }

   /*
   * Get the value of a specialized sweep parameter.
   */
   double MixAndMatchEnv::getParameter(std::string name, DArray<int> ids,
                                       bool& success) const
   {
      double val(0);
      success = false;
      if (fieldGenPtr1_) {
         val = fieldGenPtr1_->getParameter(name, ids, success);
      }
      if ((!success) && (fieldGenPtr2_)) {
         val = fieldGenPtr2_->getParameter(name, ids, success);
      }
      return val;
   }

   /*
   * Get the first field generator by const reference.
   */
   FieldGenerator const & MixAndMatchEnv::fieldGenerator1() const
   {  
      UTIL_CHECK(fieldGenPtr1_);  
      return *fieldGenPtr1_; 
   }

   /*
   * Get the second field generator (if any).
   */
   FieldGenerator const & MixAndMatchEnv::fieldGenerator2() const
   { 
      UTIL_CHECK(fieldGenPtr2_);  
      return *fieldGenPtr2_; 
   }

   /*
   * Does a second FieldGenerator child object exist?
   */
   bool MixAndMatchEnv::hasFieldGenerator2() const
   {  return (bool) fieldGenPtr2_; }

}

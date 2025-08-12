#ifndef PSCF_MIX_AND_MATCH_ENV_TMPL_TPP
#define PSCF_MIX_AND_MATCH_ENV_TMPL_TPP

/*
* PSCF - Polymer Self-Consistent Field
*
* Copyright 2015 - 2025, The Regents of the University of Minnesota
* Distributed under the terms of the GNU General Public License.
*/

#include <util/param/Label.h>

namespace Pscf {

   /*
   * Constructor.
   */
   template <class Env, class FG>
   MixAndMatchEnvTmpl<Env,FG>::MixAndMatchEnvTmpl()
    : fieldGenPtr1_(nullptr),
      fieldGenPtr2_(nullptr)
   {  setClassName("MixAndMatchEnvTmpl"); }

   /*
   * Destructor.
   */
   template <class Env, class FG>
   MixAndMatchEnvTmpl<Env,FG>::~MixAndMatchEnvTmpl()
   {
      if (fieldGenPtr1_ != nullptr) {
         delete fieldGenPtr1_;
      }
      if (fieldGenPtr2_ != nullptr) {
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
   template <class Env, class FG>
   void MixAndMatchEnvTmpl<Env,FG>::readParameters(std::istream& in)
   {
      bool generatesMask(false), generatesExt(false);

      // Before reading parameters, create FieldGenerator objects
      createGenerators();

      // Save current istream position (allows for possible rewind to here)
      std::streampos pos = in.tellg();

      // Read first FieldGenerator
      if (fieldGenPtr1_ != nullptr) {

         UTIL_CHECK(!fieldGenPtr1_->isDependent());

         // Make fieldGenPtr1_ a child ParamComponent of this object
         bool next = false; // Set indentation to be same as that of parent
         addParamComposite(*fieldGenPtr1_, next);

         // Read parameters for this FieldGenerator.
         // Reads body without indentation or curly bracket delimiters
         fieldGenPtr1_->readParameters(in);

         if (fieldGenPtr1_->type() == FG::Mask) {
            generatesMask = true;
         } else if (fieldGenPtr1_->type() == FG::Mask) {
            generatesExt = true;
         } else {
            UTIL_THROW("fieldGenPtr1_ must have type Mask or External.");
         }

      } else {

         UTIL_THROW("Object must contain at least one FieldGenerator.");

      }

      // Read second FieldGenerator (optional)
      if (fieldGenPtr2_ != nullptr) {

         // Check that one FieldGenerator is a Mask and other is External
         if (fieldGenPtr2_->type() == FG::External) {
            generatesExt = true;
            UTIL_CHECK(fieldGenPtr1_->type() == FG::Mask);
         } else if (fieldGenPtr2_->type() == FG::Mask) {
            generatesMask = true;
            UTIL_CHECK(fieldGenPtr1_->type() == FG::External);
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

      setGenerateBools(generatesMask, generatesExt);
   }

   /*
   * Check if fields need to be (re)generated. If so, generates them.
   */
   template <class Env, class FG>
   void MixAndMatchEnvTmpl<Env,FG>::generate()
   {
      if (!needsUpdate()) return;

      if (fieldGenPtr1_ != nullptr) fieldGenPtr1_->generate();
      if (fieldGenPtr2_ != nullptr) fieldGenPtr2_->generate();
      
      setNeedsUpdateFalse();
   }

   /*
   * Return specialized sweep parameter types to add to the Sweep object.
   */
   template <class Env, class FG>
   GArray<ParameterType> MixAndMatchEnvTmpl<Env,FG>::getParameterTypes()
   {
      GArray<ParameterType> a1, a2;

      if (fieldGenPtr1_ != nullptr) a1 = fieldGenPtr1_->getParameterTypes();
      if (fieldGenPtr2_ != nullptr) a2 = fieldGenPtr2_->getParameterTypes();

      for (int i = 0; i < a2.size(); i++) {
         a1.append(a2[i]);
      }

      return a1;
   }

   /*
   * Set the value of a specialized sweep parameter.
   */
   template <class Env, class FG>
   void MixAndMatchEnvTmpl<Env,FG>::setParameter(std::string name, 
                                                 DArray<int> ids,
                                                 double value, 
                                                 bool& success)
   {
      success = false;
      if (fieldGenPtr1_ != nullptr) {
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
   template <class Env, class FG>
   double MixAndMatchEnvTmpl<Env,FG>::getParameter(std::string name, 
                                                   DArray<int> ids,
                                                   bool& success) const
   {
      double val(0.0);
      success = false;
      if (fieldGenPtr1_ != nullptr) {
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
   template <class Env, class FG>
   FG const & MixAndMatchEnvTmpl<Env,FG>::fieldGenerator1() 
   const
   {  
      UTIL_CHECK(fieldGenPtr1_);  
      return *fieldGenPtr1_; 
   }

   /*
   * Get the second field generator (if any).
   */
   template <class Env, class FG>
   FG const & MixAndMatchEnvTmpl<Env,FG>::fieldGenerator2() 
   const
   { 
      UTIL_CHECK(fieldGenPtr2_);  
      return *fieldGenPtr2_; 
   }

   /*
   * Does a second FieldGenerator child object exist?
   */
   template <class Env, class FG>
   bool MixAndMatchEnvTmpl<Env,FG>::hasFieldGenerator2() const
   {  return (bool) fieldGenPtr2_; }

}
#endif

#ifndef PSCF_MIX_AND_MATCH_ENV_TMPL_TPP
#define PSCF_MIX_AND_MATCH_ENV_TMPL_TPP

/*
* PSCF - Polymer Self-Consistent Field
*
* Copyright 2015 - 2026, David Morse (morse012@um.edu) and co-authors
* Distributed under the terms of the GNU General Public License.
*/

#include <util/param/Label.h>
#include <istream>
#include <streambuf>
#include <string>

namespace Pscf {

   namespace MixAndMatchEnvDetail {

      /*
      * Stream buffer that returns a given text, then the characters of
      * another stream buffer, and records the characters that it takes
      * from that other stream buffer.
      */
      class ReplayStreamBuf : public std::streambuf
      {
      public:

         /*
         * Constructor.
         *
         * \param source  stream buffer read after the text
         * \param text  text returned first
         */
         ReplayStreamBuf(std::streambuf* source, 
                         std::string const & text = "")
          : text_(text),
            sourcePtr_(source)
         {
            if (!text_.empty()) {
               char* begin = &text_[0];
               setg(begin, begin, begin + text_.size());
            }
         }

         /*
         * Get the characters taken from the source so far.
         */
         std::string const & record() const
         {  return record_; }

      protected:

         // Get the next character of the source, without extracting it
         int_type underflow()
         {  return sourcePtr_->sgetc(); }

         // Extract and record the next character of the source
         int_type uflow()
         {
            int_type c = sourcePtr_->sbumpc();
            if (!traits_type::eq_int_type(c, traits_type::eof())) {
               record_ += traits_type::to_char_type(c);
            }
            return c;
         }

      private:

         std::string text_;
         std::string record_;
         std::streambuf* sourcePtr_;

      };

   }

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
   * If the second generator is "dependent", it re-reads the block of the
   * first generator. The text read by the first generator is recorded and
   * given again to the second one, followed by the rest of the stream, so
   * that the stream is never rewound (a pipe cannot be).
   */
   template <class Env, class FG>
   void MixAndMatchEnvTmpl<Env,FG>::readParameters(std::istream& in)
   {
      bool generatesMask(false), generatesExt(false);

      // Before reading parameters, create FieldGenerator objects
      createGenerators();

      // If the second FieldGenerator is dependent, the first one reads
      // through a stream that records the text it reads
      bool isDependent2 = (fieldGenPtr2_ != nullptr)
                          && fieldGenPtr2_->isDependent();
      MixAndMatchEnvDetail::ReplayStreamBuf recordBuf(in.rdbuf());
      std::istream recordIn(&recordBuf);
      std::istream& in1 = isDependent2 ? recordIn : in;

      // Read first FieldGenerator
      if (fieldGenPtr1_ != nullptr) {

         UTIL_CHECK(!fieldGenPtr1_->isDependent());

         // Make fieldGenPtr1_ a child ParamComponent of this object
         bool next = false; // Set indentation to be same as that of parent
         addParamComposite(*fieldGenPtr1_, next);

         // Read parameters for this FieldGenerator.
         // Reads body without indentation or curly bracket delimiters
         fieldGenPtr1_->readParameters(in1);

         if (fieldGenPtr1_->type() == FG::Mask) {
            generatesMask = true;
         } else if (fieldGenPtr1_->type() == FG::External) {
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

         // Make this FieldGenerator a child ParamComponent of this object
         bool next = false; // Set indentation to be same as that of parent
         addParamComposite(*fieldGenPtr2_, next);

         // Read parameters for this FieldGenerator.
         // Reads body without indentation or curly bracket delimiters
         if (isDependent2) {
            // Read the text read by the first FieldGenerator, then the
            // rest of the stream. (A label left in the Label buffer by 
            // the first FieldGenerator is part of that text.)
            Label::clear();
            MixAndMatchEnvDetail::ReplayStreamBuf 
                                replayBuf(in.rdbuf(), recordBuf.record());
            std::istream replayIn(&replayBuf);
            fieldGenPtr2_->readParameters(replayIn);
         } else {
            fieldGenPtr2_->readParameters(in);
         }
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

      // Route these parameters through this Environment, so that a
      // change calls setParameter() below, which calls reset(), and
      // the fields are regenerated
      for (int i = 0; i < a1.size(); i++) {
         a1[i].modifierPtr = this;
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

/*
* PSCF - Polymer Self-Consistent Field Theory
*
* Copyright 2016 - 2022, The Regents of the University of Minnesota
* Distributed under the terms of the GNU General Public License.
*/

#include "ImposedFieldsTmpl.h"

namespace Pscf {

   // Constructor
   ImposedFieldsTmpl::ImposedFieldsTmpl()
    : fieldGenPtr1_(0),
      fieldGenPtr2_(0),
      type_()
   {}

   // Destructor
   ImposedFieldsTmpl::~ImposedFieldsTmpl()
   {}

   // Read parameters from input stream
   void ImposedFieldsTmpl::readParameters(std::istream& in)
   {
      // read type name and use it to create the generator objects
      read(in, "type", type_);
      createGenerators();

      // Read FieldGenerator for mask (optional)
      if (fieldGenPtr1_) {

         // Make fieldGenPtr1_ a child paramComponent of this object, so that 
         // it will be read/written correctly to/from param file with correct 
         // indentation
         setParent(*fieldGenPtr1_,false);
         addComponent(*fieldGenPtr1_,false);

         // Read parameters for this FieldGenerator
         fieldGenPtr1_->readParameters(in);
      }

      // Read second FieldGenerator (optional)
      if (fieldGenPtr2_) {

         // Make fieldGenPtr2_ a child paramComponent of this object, so that 
         // it will be read/written correctly to/from param file with correct 
         // indentation
         setParent(*fieldGenPtr2_,false);
         addComponent(*fieldGenPtr2_,false);

         // Check that one of the FieldGenerator objects has type Mask and
         // the other has type External
         if (fieldGenPtr2_->type() == FieldGenerator::External) {
            UTIL_CHECK(fieldGenPtr1_->type() == FieldGenerator::Mask);
         } else if (fieldGenPtr2_->type() == FieldGenerator::Mask) {
            UTIL_CHECK(fieldGenPtr1_->type() == 
                                            FieldGenerator::External);
         } else {
            UTIL_THROW("fieldGenPtr2_ must have type Mask or External.");
         }
         
         // Read parameters for external fields
         fieldGenPtr2_->readParameters(in);
      }
   }

   // Allocate, check compatibility, calculate, and store the field(s)
   void ImposedFieldsTmpl::setup()
   {
      if (fieldGenPtr1_) fieldGenPtr1_->setup();
      if (fieldGenPtr2_) fieldGenPtr2_->setup();
   }

   // Check whether system has changed and update the field(s) if necessary
   void ImposedFieldsTmpl::update()
   {
      if (fieldGenPtr1_) fieldGenPtr1_->update();
      if (fieldGenPtr2_) fieldGenPtr2_->update();
   }

   // Return specialized sweep parameter types to add to the Sweep object
   GArray<ParameterType> ImposedFieldsTmpl::getParameterTypes()
   {
      GArray<ParameterType> a1, a2;

      if (fieldGenPtr1_) a1 = fieldGenPtr1_->getParameterTypes();
      if (fieldGenPtr2_) a2 = fieldGenPtr2_->getParameterTypes();

      for (int i = 0; i < a2.size(); i++) {
         a1.append(a2[i]);
      }

      return a1;
   }

   // Set the value of a specialized sweep parameter
   void ImposedFieldsTmpl::setParameter(std::string name, 
                                                  DArray<int> ids, 
                                                  double value, bool& success)
   {
      success = false;
      if (fieldGenPtr1_) {
         fieldGenPtr1_->setParameter(name, ids, value, success);
      }
      if ((!success) && (fieldGenPtr2_)) {
         fieldGenPtr2_->setParameter(name, ids, value, success);
      }
   }

   // Get the value of a specialized sweep parameter
   double ImposedFieldsTmpl::getParameter(std::string name, 
                                                    DArray<int> ids, 
                                                    bool& success)
   const
   {
      double val;
      success = false;
      if (fieldGenPtr1_) {
         val = fieldGenPtr1_->getParameter(name, ids, success);
      }
      if ((!success) && (fieldGenPtr2_)) {
         val = fieldGenPtr2_->getParameter(name, ids, success);
      }
      return val;
   }

   // Get the type string associated with this object
   std::string ImposedFieldsTmpl::type() const
   {  return type_; }
}
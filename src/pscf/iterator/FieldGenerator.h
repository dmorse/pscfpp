#ifndef PSCF_IMPOSED_FIELD_GENERATOR_H
#define PSCF_IMPOSED_FIELD_GENERATOR_H

/*
* PSCF - Polymer Self-Consistent Field Theory
*
* Copyright 2016 - 2022, The Regents of the University of Minnesota
* Distributed under the terms of the GNU General Public License.
*/

#include <pscf/sweep/ParameterModifier.h> // base class
#include <util/param/ParamComposite.h>    // base class

namespace Pscf {

   using namespace Util;

   /**
   * Abstract base class for objects that generate fields for ImposedFields.
   * 
   * A FieldGenerator subclass will contain code to generate a mask (which
   * imposes geometric confinement), a set of external fields (one for each
   * monomer species), or both, based on a set of input parameters defined
   * by the subclass. The interface for the FieldGenerator to exchange this
   * information with the System object is contained within the class 
   * ImposedFieldsTmpl, which contains up to two FieldGenerator objects as
   * members.
   * 
   * \ingroup Pscf_Iterator_Module
   */
   class FieldGenerator : public ParamComposite, 
                                 public ParameterModifier
   {

   public:

      /**
      * Enum representing the type of field (mask, external field, or both).
      */
      enum Type {Mask, External, Both, None};
      
      /**
      * Constructor
      */
      FieldGenerator();

      /**
      * Destructor
      */
      ~FieldGenerator();
      
      /**
      * Allocate, check compatibility, calculate, and store the field(s)
      */
      void setup();

      /**
      * Check whether system has changed and update the field(s) if necessary
      */
      void update();

      /**
      * Check that the system is compatible with these fields
      * 
      * This may, for example, involve checking that the fields satisfy
      * the system's symmetry and lattice system, and ensuring that the
      * unit cell is only allowed to adjust its size in ways that agree
      * with the design of the fields.
      */
      virtual void checkCompatibility() = 0;

      /**
      * Check whether system has changed such that the field(s) need updating
      */
      virtual bool updateNeeded() const = 0;

      /**
      * Check whether the field(s) have been generated
      */
      virtual bool isGenerated() const = 0;

      /**
      * Return Type enumeration value (Mask, External, or None)
      *
      * This value should be initialized by subclasses during construction.
      */
      Type type() const
      {  return type_; }
   
   protected:

      /**
      * Allocate container(s) necessary to generate and store field(s)
      */ 
      virtual void allocate() = 0;

      /**
      * Generate the field(s) and store where the Iterator can access
      */
      virtual void generateMask() = 0;

      /**
      * Type of field (Mask, External, Both, or None)
      * 
      * This parameter should be a private member of subclasses, and should
      * be set in the constructor and left unchanged for the remainder of
      * the calculation.
      */
      Type type_;

   };

}
#endif
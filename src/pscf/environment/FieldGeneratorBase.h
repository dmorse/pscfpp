#ifndef PSCF_FIELD_GENERATOR_BASE_H
#define PSCF_FIELD_GENERATOR_BASE_H

/*
* PSCF - Polymer Self-Consistent Field
*
* Copyright 2015 - 2025, The Regents of the University of Minnesota
* Distributed under the terms of the GNU General Public License.
*/

#include <util/param/ParamComposite.h>    // base class
#include <pscf/sweep/ParameterModifier.h> // base class

namespace Pscf {

   using namespace Util;

   /**
   * Base class for objects that generate fields for a MixAndMatchEnv.
   * 
   * A subclass of FieldGeneratorBase will contain code to generate a 
   * mask (which imposes geometric confinement) or a set of external 
   * potential fields (one for each monomer type), based on a set of 
   * input parameters defined by the subclass. The interface for the 
   * subclass to exchange this information with the System object is 
   * contained within the class MixAndMatchEnvTmpl, which owns up to
   * two FieldGeneratorBase objects as members.
   * 
   * This class is a subclass of ParameterModifier so that, if desired,
   * subclasses may declare parameters that can be varied in a Sweep.
   * 
   * In their constructors, subclasses should set the variable type_
   * to the appropriate value.
   * 
   * \ingroup Pscf_Environment_Module
   */
   class FieldGeneratorBase : public ParamComposite, 
                              public ParameterModifier
   {

   public:

      /**
      * Enum representing the type of field (mask, external field, or none).
      */
      enum Type {Mask, External, None};
      
      /**
      * Constructor
      */
      FieldGeneratorBase();

      /**
      * Destructor
      */
      ~FieldGeneratorBase();
      
      /**
      * Checks if fields need to be (re)generated. If so, generates them. 
      */
      void generate();

      /**
      * Check whether system has changed such that the field(s) need updating.
      */
      virtual bool needsUpdate() const = 0;

      /**
      * Check that the system is compatible with these fields.
      * 
      * This may, for example, involve checking that the fields satisfy
      * the system's symmetry and lattice system, and ensuring that the
      * unit cell is only allowed to adjust its size in ways that agree
      * with the design of the fields.
      */
      virtual void checkCompatibility() = 0;

      /**
      * Return Type enumeration value (Mask, External, or None)
      *
      * This value should be initialized by subclasses during construction.
      */
      Type type() const;

      /**
      * Is this object dependent on parameters of another FieldGeneratorBase?
      * 
      * The parent MixAndMatchEnvTmpl object can contain up to two 
      * FieldGeneratorBase objects: a mask generator and an external field
      * generator. In some cases, the properties of one may be dependent on 
      * the properties of the other. isDependent allows a FieldGeneratorBase 
      * to indicate to the parent MixAndMatchEnvTmpl object that it needs to 
      * read the parameters of the other FieldGenerator in addition to its 
      * own. 
      * 
      * During readParameters, the parent object will first allow the 
      * independent FieldGeneratorBase to read its own parameters. It will 
      * then rewind the istream to allow the dependent FieldGeneratorBase
      * to read the parameters of the object on which it is dependent,
      * followed by its own parameters.
      * 
      * Therefore, a dependent FieldGeneratorBase must be the second of two
      * FieldGeneratorBases stored in a parent MixAndMatchEnvTmpl, and two
      * FieldGeneratorBases may not be dependent on each other. If the mask
      * and external field are dependent on each other, a different 
      * subclass of Environment should be used instead of a 
      * MixAndMatchEnvTmpl to generate these imposed fields.
      */
      bool isDependent() const;
   
   protected:

      /**
      * Compute the field(s) and store where the System can access.
      */
      virtual void compute() = 0;

      /**
      * Type of field (Mask, External, or None).
      * 
      * This parameter should be set in the constructor of subclasses 
      * and left unchanged for the remainder of the calculation.
      */
      Type type_;

      /// Is this object dependent on parameters of another FieldGeneratorBase?
      bool isDependent_;

   };

   // Inline member functions

   // Return Type enumeration value (Mask, External, or None)
   inline FieldGeneratorBase::Type FieldGeneratorBase::type() const
   {  return type_; }

   // Is this object dependent on parameters of another FieldGeneratorBase?
   inline bool FieldGeneratorBase::isDependent() const
   {  return isDependent_; }

}
#endif
#ifndef PSCF_FIELD_GENERATOR_H
#define PSCF_FIELD_GENERATOR_H

/*
* PSCF - Polymer Self-Consistent Field Theory
*
* Copyright 2015 - 2025, The Regents of the University of Minnesota
* Distributed under the terms of the GNU General Public License.
*/

#include <pscf/sweep/ParameterModifier.h> // base class
#include <util/param/ParamComposite.h>    // base class

namespace Pscf {

   using namespace Util;

   /**
   * Base class for objects that generate fields for MixAndMatchEnv.
   * 
   * A FieldGenerator subclass will contain code to generate a mask (which
   * imposes geometric confinement) or a set of external fields (one for 
   * each monomer species), based on a set of input parameters defined
   * by the subclass. The interface for the FieldGenerator to exchange this
   * information with the System object is contained within the class 
   * MixAndMatchEnv, which contains up to two FieldGenerator objects as
   * members.
   * 
   * \ingroup Pscf_Environment_Module
   */
   class FieldGenerator : public ParamComposite, 
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
      * Check whether system has changed, update the field(s) if necessary
      */
      void update();

      /**
      * Is this object dependent on the parameters of another FieldGenerator?
      * 
      * The parent MixAndMatchEnv object can contain up to two 
      * FieldGenerator objects: a mask and an external field. In some cases, 
      * the properties of one may be dependent on the properties of the 
      * other. isDependent allows a FieldGenerator to indicate to the 
      * parent object that it needs to read the parameters of the other
      * FieldGenerator in addition to its own. 
      * 
      * During readParameters, the parent object will first allow the 
      * independent FieldGenerator to read its own parameters. It will 
      * then rewind the istream to allow the dependent FieldGenerator
      * to read the parameters of the object on which it is dependent,
      * followed by its own parameters.
      * 
      * Therefore, a dependent FieldGenerator must be the second of two
      * FieldGenerators stored in a parent MixAndMatchEnv, and two
      * FieldGenerators may not be dependent on each other. In such a
      * circumstance, one should not use a MixAndMatchEnv, and should 
      * instead use a different subclass of Environment.
      */
      bool isDependent() const;

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
      * Get contribution to the stress from this imposed field.
      * 
      * If the imposed field(s) change in a non-affine manner under changes
      * in the lattice parameters, then the "stress" used to optimize the
      * lattice parameters must contain additional terms arising from the
      * imposed field(s). Thus, this method may return a nonzero value 
      * regardless of the Type of this object. A return value of zero 
      * indicates that the imposed field(s) stretch affinely under a change 
      * in the given lattice parameter. 
      * 
      * The default implementation of the method raises an error. If a
      * FieldGenerator subclass is designed to be used only in a rigid 
      * box, then the subclass does not need to redefine this method,
      * as it will not ever be called. But subclasses must define the 
      * method if they are to be used with non-rigid boxes, or else the
      * error message in this default implementation will be triggered.
      * 
      * \param paramId  index of the lattice parameter being varied
      */
      virtual double stress(int paramId) const;

      /**
      * Modify stress to minimize a property other than fHelmholtz. 
      * 
      * It may be preferable with certain imposed fields to minimize a  
      * property other than fHelmholtz with respect to the lattice  
      * parameters. For instance, in a thin film it is useful to minimize 
      * the excess free energy per unit area, (fHelmholtz - fRef) * Delta, 
      * where fRef is a reference free energy and Delta is the film 
      * thickness. The information needed to perform such a modification
      * is often contained within the FieldGenerator objects. Therefore,
      * this method allows subclasses of FieldGenerator to modify the 
      * stress.
      * 
      * The method is called by the MixAndMatchEnv object that owns this
      * object, and the return value should be used to optimize the 
      * lattice parameters. 
      * 
      * By default, this method will simply return the value of stress
      * that is provided as an input, without performing a modification,
      * which corresponds to a stress that will be used to minimize 
      * fHelmholtz.
      * 
      * \param paramId  index of the lattice parameter with this stress
      * \param stress  stress value calculated by Mixture object
      */
      virtual double modifyStress(int paramId, double stress) const;

      /**
      * Return Type enumeration value (Mask, External, or None)
      *
      * This value should be initialized by subclasses during construction.
      */
      Type type() const;
   
   protected:

      /**
      * Allocate container(s) necessary to generate and store field(s)
      */ 
      virtual void allocate() = 0;

      /**
      * Generate the field(s) and store where the System can access
      */
      virtual void generate() = 0;

      /**
      * Type of field (Mask, External, or None)
      * 
      * This parameter should be a private member of subclasses, and should
      * be set in the constructor and left unchanged for the remainder of
      * the calculation.
      */
      Type type_;

      /// Is this object dependent on the parameters of another FieldGenerator?
      bool isDependent_;

   };

   // Inline member functions

   // Get contribution to the stress from this imposed field.
   inline double FieldGenerator::stress(int paramId) const
   {  UTIL_THROW("Unimplemented stress() method called."); } 

   // Modify stress to minimize a property other than fHelmholtz. 
   inline double FieldGenerator::modifyStress(int paramId, double stress) 
   const
   {  return stress; }

   // Is this object dependent on the parameters of another FieldGenerator?
   inline bool FieldGenerator::isDependent() const
   {  return isDependent_; }

   // Return Type enumeration value (Mask, External, or None)
   inline FieldGenerator::Type FieldGenerator::type() const
   {  return type_; }

}
#endif
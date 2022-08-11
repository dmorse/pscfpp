#ifndef PSPC_ITERATOR_H
#define PSPC_ITERATOR_H

/*
* PSCF - Polymer Self-Consistent Field Theory
*
* Copyright 2016 - 2019, The Regents of the University of Minnesota
* Distributed under the terms of the GNU General Public License.
*/

#include <util/param/ParamComposite.h>    // base class
#include <util/containers/DArray.h>
#include <util/containers/FSArray.h>
#include <util/global.h>                  

namespace Pscf {
namespace Pspc
{

   template <int D>
   class System;

   using namespace Util;

   typedef DArray<double> FieldCPU;

   /**
   * Base class for iterative solvers for SCF equations.
   *
   * \ingroup Pspc_Iterator_Module
   */
   template <int D>
   class Iterator : public ParamComposite
   {

   public:

      /**
      * Default constructor.
      */
      Iterator();

      /**
      * Constructor.
      * 
      * \param system parent System object
      */
      Iterator(System<D>& system);

      /**
      * Destructor.
      */
      ~Iterator();

      /**
      * Setup iterator.
      */
      virtual void setup() = 0;

      /**
      * Iterate to solution.
      *
      * \return error code: 0 for success, 1 for failure.
      */
      virtual int solve() = 0;

      /**
      * Return const reference to parent system.
      */
      System<D> const & system() const
      {  return *sysPtr_;}

      /**
      * Return true if unit cell is flexible, false if rigid.
      */
      virtual bool isFlexible() const 
      {  return isFlexible_;}

      /**
      * Sets the flexibleParams_ array.
      * 
      * \param flexParams array of indices of each flexible lattice parameter
      */
      void setFlexibleParams(FSArray<int, 6> const & flexParams)
      {  
         flexibleParams_ = flexParams; 
         if (flexibleParams_.size() == 0) {
            isFlexible_ = false;
         } else {
            isFlexible_ = true;
         }
      }

      /**
      * Return an array containing indices of each flexible lattice 
      * parameter, e.g. [0,1,2] if there are 3 lattice parameters and
      * all 3 are flexible.
      */
      FSArray<int, 6> flexibleParams() const
      {  return flexibleParams_; }

      /**
      * Accepts a field array representing the concentration profile
      * of the mask that will be imposed upon the unit cell, and stores 
      * it for use during iteration.
      * 
      * \param field mask concentration profile, in basis format
      */
      void setMask(FieldCPU const & field);

      /**
      * Accepts an array of fields representing the external potential
      * field felt by each monomer species, and stores it for use 
      * during iteration.
      * 
      * \param fields external fields for each species, in basis format
      */
      void setExternalFields(DArray<FieldCPU> const & fields);

      /**
      * Returns the field that represents the mask imposed upon the 
      * unit cell, or returns an unallocated field if the iterator 
      * does not have a mask. Field is in symmetry-adapted basis format.
      */
      virtual FieldCPU const & maskField() const;

      /**
      * Returns the array of fields that represents the external potential
      * field felt by each monomer species, or returns an unallocated field
      * if the iterator does not have an external field. Fields are in 
      * symmetry-adapted basis format.
      */
      virtual DArray<FieldCPU> const & externalFields() const;

      /**
      * Returns the field that represents the external field imposed
      * on the monomer species monomerId, or returns an unallocated
      * field if the iterator does not have an external field. Field
      * is in symmetry-adapted basis format.
      * 
      * \param monomerId integer monomer type index
      */
      virtual FieldCPU const & externalField(int monomerId) const;

      /**
      * Return true if this iterator imposes a mask within the unit cell,
      * false if no mask is imposed.
      */
      virtual bool hasMask() const;

      /**
      * Return true if this iterator imposes an external field on any
      * monomer species, returns false if no external field is present.
      */
      virtual bool hasExternalField() const;

   protected:

      /**
      * Return reference to parent system.
      */
      System<D>& system() 
      {  return *sysPtr_;}

      /// Is the unit cell flexible during iteration?
      bool isFlexible_;

      /// Array of indices of the lattice parameters that are flexible
      FSArray<int, 6> flexibleParams_;

      /// Does this iterator have a mask?
      bool hasMask_;

      /// Does this iterator have an external field?
      bool hasExternalField_;

   private:

      /// Pointer to the associated system object.
      System<D>* sysPtr_;

      /// Field representing the mask imposed on the unit cell
      DArray<double> maskField_;

      /// Array containing the external fields felt by each species
      DArray<DArray<double> > externalFields_;
      
   };

   // Inline member functions

   // Default constructor
   template <int D>
   inline Iterator<D>::Iterator()
    : hasMask_(false),
      hasExternalField_(false),
      maskField_(),
      externalFields_()
   {  setClassName("Iterator"); }

   // Constructor
   template <int D>
   inline Iterator<D>::Iterator(System<D>& system)
    : hasMask_(false),
      hasExternalField_(false),
      sysPtr_(&system),
      maskField_(),
      externalFields_()
   {  setClassName("Iterator"); }

   // Destructor
   template <int D>
   inline Iterator<D>::~Iterator()
   {}

   // Store the concentration field for the mask imposed on the unit cell
   template <int D>
   inline void Iterator<D>::setMask(FieldCPU const & field)
   {  
      maskField_ = field; // copy field into maskField_
      hasMask_ = true;
   }

   // Store the external fields imposed on each monomer species
   template <int D> 
   inline void Iterator<D>::setExternalFields(DArray<FieldCPU> const & fields)
   {  
      externalFields_ = fields; // copy fields into externalFields_
      hasExternalField_ = true;
   }

   // Get the concentration field for the mask imposed on the unit cell
   template <int D>
   inline FieldCPU const & Iterator<D>::maskField() const
   {  return maskField_; }

   // Get array of all external fields felt by monomer species
   template <int D>
   inline DArray<FieldCPU> const & Iterator<D>::externalFields() 
   const
   {  return externalFields_; }

   // Get the external field felt by one monomer species
   template <int D>
   inline FieldCPU const & Iterator<D>::externalField(int monomerId) 
   const
   {  return externalFields_[monomerId]; }

   // Does this iterator have a mask?
   template <int D>
   inline bool Iterator<D>::hasMask() const
   {  return hasMask_; }

   // Does this iterator have any external field?
   template <int D>
   inline bool Iterator<D>::hasExternalField() const
   {  return hasExternalField_; }

} // namespace Pspc
} // namespace Pscf
#endif

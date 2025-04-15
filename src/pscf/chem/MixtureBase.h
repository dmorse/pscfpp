#ifndef PSCF_MIXTURE_BASE_H
#define PSCF_MIXTURE_BASE_H

/*
* PSCF - Polymer Self-Consistent Field Theory
*
* Copyright 2016 - 2022, The Regents of the University of Minnesota
* Distributed under the terms of the GNU General Public License.
*/

#include <pscf/chem/Monomer.h>
#include <util/containers/DArray.h>

namespace Pscf
{

   class PolymerSpecies;
   class SolventSpecies;

   using namespace Util;

   /**
   * Interface for a Mixture descriptor.
   *
   * \ingroup Pscf_Chem_Module
   */
   class MixtureBase 
   {

   public:

      /**
      * Constructor.
      */
      MixtureBase();

      /**
      * Destructor.
      */
      ~MixtureBase();

      /**
      * Set new vMonomer value.
      *
      * \param vMonomer  new vMonomer value
      */
      void setVmonomer(double vMonomer);

      /// \name Accessors (by value)
      ///@{

      /**
      * Get number of monomer types.
      */
      int nMonomer() const;

      /**
      * Get number of polymer species.
      */
      int nPolymer() const;

      /**
      * Get number of solvent (point particle) species.
      */
      int nSolvent() const;

      /** 
      * Get number of total blocks in the mixture across all polymers.
      */
      int nBlock() const;

      /**
      * Get monomer reference volume (set to 1.0 by default).
      */
      double vMonomer() const;

      ///@}
      /// \name Accessors (by const reference)
      ///@{

      /**
      * Get a Monomer type descriptor (const reference).
      *
      * \param id integer monomer type index (0 <= id < nMonomer)
      */
      Monomer const & monomer(int id) const;

      /**
      * Get a PolymerSpecies (const reference).
      *
      * \param id integer polymer species index (0 <= id < nPolymer)
      */
      virtual
      PolymerSpecies const & polymerSpecies(int id) const = 0;

      /**
      * Set a solvent solver object by const reference.
      *
      * \param id integer solvent species index (0 <= id < nSolvent)
      */
      virtual
      SolventSpecies const & solventSpecies(int id) const = 0;

      ///@}

   protected:

      /**
      * Array of monomer type descriptors.
      */
      DArray<Monomer> monomers_;

      /**
      * Number of monomer types.
      */
      int nMonomer_; 

      /**
      * Number of polymer species.
      */
      int nPolymer_;

      /**
      * Number of solvent species.
      */
      int nSolvent_;

      /**
      * Number of blocks total, across all polymers.
      */
      int nBlock_;

      /**
      * Monomer reference volume (set to 1.0 by default).
      */
      double vMonomer_;

      /**
      * Get a Monomer type descriptor (non-const reference).
      *
      * \param id integer monomer type index (0 <= id < nMonomer)
      */
      Monomer& monomer(int id);

   };

   // Inline member functions

   inline int MixtureBase::nMonomer() const
   {  return nMonomer_; }

   inline int MixtureBase::nPolymer() const
   {  return nPolymer_; }

   inline int MixtureBase::nSolvent() const
   {  return nSolvent_; }

   inline int MixtureBase::nBlock() const
   {  return nBlock_; }


   inline 
   Monomer const & MixtureBase::monomer(int id) const
   {  
      UTIL_CHECK(id < nMonomer_);
      return monomers_[id]; 
   }

   inline 
   Monomer& MixtureBase::monomer(int id)
   {  
      UTIL_CHECK(id < nMonomer_);
      return monomers_[id]; 
   }

   inline 
   double MixtureBase::vMonomer() const
   {  return vMonomer_; }

}
#endif

#ifndef PSCF_MIXTURE_BASE_H
#define PSCF_MIXTURE_BASE_H

/*
* PSCF - Polymer Self-Consistent Field Theory
*
* Copyright 2015 - 2025, The Regents of the University of Minnesota
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
   * Abstract descriptor for a mixture of polymer and solvent species.
   *
   * A MixtureBase has an array of Monomer objects and provides access to
   * PolymerSpecies and SolventSpecies objects describing all molecular
   * species in in a mixture.  The MixtureBase class does not have 
   * functions or data structures needed to solve the modified diffusion 
   * equation (MDE), and is thus a descriptor but not a solver for the 
   * mixture. 
   *
   * MixtureBase is an abstract base class for "solver" classes that can 
   * actually solve the single-molecule statistical mechanics problem for 
   * very species in a mixture. Each implementation level sub-namespace 
   * of Pscf (R1d, Rpc or Rpg) contains a concrete class named Mixture 
   * that is derived from Pscf::MixtureBase, and that acts as both a 
   * solver and descriptor for the mixture.  Each such subspace also 
   * defines an polymer MDE solver class named Polymer that is a subclass 
   * of PolymerSpecies and a solvent solver class named Solvent that is 
   * subclass of SolventSpecies.  
   * 
   * The Mixture class in each such implementation-level namespace is a 
   * subclass of an instantiation Pscf::PolymerTmpl<Polymer, Solvent> of 
   * the class template Pscf::PolymerTmpl. This template is derived 
   * directly from MixtureBase. A PolymerTmpl<Polymer, Solvent> object has 
   * two member private variables that are arrays of Polymer and Solvent 
   * solver objects.  The PolymerTmpl template defines implementations of 
   * functions polymerSpecies(int id) and solventSpecies(int id) that are
   * declared as pure virtual member functions of MixtureBase. These
   * functions return a single Polymer solver object as a reference to
   * a PolymerSpecies descriptor, or a Solvent solver object as a 
   * reference to a SolventSpecies descriptor, respectively.
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
      * Get total number blocks among all polymer species.
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
      * Get a Monomer type descriptor by const reference.
      *
      * \param id integer monomer type index (0 <= id < nMonomer)
      */
      Monomer const & monomer(int id) const;

      /**
      * Get a PolymerSpecies by const reference.
      *
      * \param id  integer polymer species index (0 <= id < nPolymer)
      */
      virtual
      PolymerSpecies const & polymerSpecies(int id) const = 0;

      /**
      * Set a solvent solver object by const reference.
      *
      * \param id  integer solvent species index (0 <= id < nSolvent)
      */
      virtual
      SolventSpecies const & solventSpecies(int id) const = 0;

      ///@}

   protected:

      /**
      * Get a Monomer type descriptor (non-const reference).
      *
      * \param id integer monomer type index (0 <= id < nMonomer)
      */
      Monomer& monomer(int id);

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

   };

   // Inline public member functions

   inline int MixtureBase::nMonomer() const
   {  return nMonomer_; }

   inline int MixtureBase::nPolymer() const
   {  return nPolymer_; }

   inline int MixtureBase::nSolvent() const
   {  return nSolvent_; }

   inline int MixtureBase::nBlock() const
   {  return nBlock_; }

   inline 
   double MixtureBase::vMonomer() const
   {  return vMonomer_; }

   inline 
   Monomer const & MixtureBase::monomer(int id) const
   {  
      UTIL_CHECK(id < nMonomer_);
      return monomers_[id]; 
   }

   // Inline protected member function

   inline 
   Monomer& MixtureBase::monomer(int id)
   {  
      UTIL_CHECK(id < nMonomer_);
      return monomers_[id]; 
   }

}
#endif

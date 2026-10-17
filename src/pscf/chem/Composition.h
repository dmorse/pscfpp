#ifndef PSCF_MIXTURE_BASE_H
#define PSCF_MIXTURE_BASE_H

/*
* PSCF - Polymer Self-Consistent Field
*
* Copyright 2015 - 2026, The Regents of the University of Minnesota
* Distributed under the terms of the GNU General Public License.
*/

#include <pscf/chem/Monomer.h>
#include <util/containers/DArray.h>

namespace Pscf
{

   template <typename WT> class PolymerSpecies;
   template <typename WT> class SolventSpecies;

   using namespace Util;

   /**
   * Descriptor for a mixture of polymer and solvent species.
   *
   * A Composition describes the contents of a a mixture that may 
   * contain polymer and solvent molecular species. An instance of this 
   * class contains an array of Monomer objects and provides access to 
   * PolymerSpecies and SolventSpecies objects that describe all molecular
   * species in the mixture. 
   *
   * Composition is an abstract base class for "solver" subclasses that 
   * can also solve the single-molecule statistical mechanics problem 
   * for every species in a mixture. Each program level sub-namespace 
   * of Pscf (R1d or Rp) contains a class or class template named Mixture 
   * that is derived from Pscf::Composition, and that acts as an MDE 
   * solver for the mixture, as well a descriptor. A reference to the 
   * Composition sub-object of such a Mixture solver object can be used
   * in generic algorithms that require a description of the contents of 
   * a mixture, but that do not actually need to solve the MDE. 
   * 
   * Within each program-level namespace, each Mixture class or class
   * template is derived from a specialization of the class template
   * Pscf::MixtureTmpl<class Polymer, class Solvent>. Specializations 
   * of this template are, in turn, derived from Pscf::Composition. The 
   * template parameters Polymer and Solvent of the MixtureTmpl template
   * are classes that are used to represent polymer and solvent species
   * in the relevant program level namespace, respectively.  The Polymer 
   * class or class template in each program level namespace is an MDE 
   * solver for a polymer species that is derived from the class 
   * Pscf::PolymerSpecies, which is a polymer species descriptor. 
   * Similarly, the Solvent class or class template in each program level
   * namespace is an MDE solver for a solvent species that is derived from
   * Pscf::SolventSpecies, which is a corresponding species descriptor. 
   *
   * The MixtureTmpl<Polymer, Solvent> template contains private member 
   * variables that are arrays of Polymer and Solvent solver objects.  
   * The MixtureTmpl template also defines implementations of the virtual
   * functions polymerSpecies(int id) and solventSpecies(int id), which 
   * are declared as pure virtual member functions of Composition. These 
   * functions return const references to the PolymerSpecies sub-object
   * of a Polymer solver object or to the SolventSpecies sub-object of 
   * of a Solvent solver object, respectively.
   * 
   * \ingroup Pscf_Chem_Module
   */
   template <typename WT>
   class Composition 
   {

   public:

      /**
      * Constructor.
      */
      Composition();

      /**
      * Destructor.
      */
      virtual ~Composition() = default;

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

      /**
      * Is this mixture being treated in canonical ensemble?
      *
      * Returns true iff a closed ensemble is used for every polymer
      * and solve species, by specifying a volume fraction phi rather
      * than a chemical potential mu for every species in the mixture.
      */
      bool isCanonical() const;

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
      * This function must be implemented by a specialization of the
      * Pscf::MixtureTmpl <Polymer, Solvent> class template, which is
      * derived from Composition.
      *
      * \param id  integer polymer species index (0 <= id < nPolymer)
      */
      virtual
      PolymerSpecies<WT> const & polymerSpecies(int id) const = 0;

      /**
      * Set a solvent solver object by const reference.
      *
      * This function must be implemented by a specialization of the
      * Pscf::MixtureTmpl <Polymer, Solvent> class template, which is
      * derived from Composition.
      *
      *
      * \param id  integer solvent species index (0 <= id < nSolvent)
      */
      virtual
      SolventSpecies<WT> const & solventSpecies(int id) const = 0;

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

   template <typename WT> inline 
   int Composition<WT>::nMonomer() const
   {  return nMonomer_; }

   template <typename WT> inline 
   int Composition<WT>::nPolymer() const
   {  return nPolymer_; }

   template <typename WT> inline 
   int Composition<WT>::nSolvent() const
   {  return nSolvent_; }

   template <typename WT> inline 
   int Composition<WT>::nBlock() const
   {  return nBlock_; }

   template <typename WT> inline 
   double Composition<WT>::vMonomer() const
   {  return vMonomer_; }

   template <typename WT> inline 
   Monomer const & Composition<WT>::monomer(int id) const
   {  
      UTIL_CHECK(id < nMonomer_);
      return monomers_[id]; 
   }

   // Inline protected member function

   template <typename WT> inline 
   Monomer& Composition<WT>::monomer(int id)
   {  
      UTIL_CHECK(id < nMonomer_);
      return monomers_[id]; 
   }

   extern template class Composition<double>;

}
#endif

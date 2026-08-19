#ifndef PSCF_SOLVENT_SPECIES_H
#define PSCF_SOLVENT_SPECIES_H

/*
* PSCF - Polymer Self-Consistent Field
*
* Copyright 2015 - 2026, The Regents of the University of Minnesota
* Distributed under the terms of the GNU General Public License.
*/

#include <pscf/chem/Species.h>             // base class

namespace Pscf {

   using namespace Util;

   /**
   * Descriptor for a solvent species.
   *
   * A SolventSpecies has values for monomer type index (monomerId) and
   * steric volume per molecule of a solvent species (size), in addition
   * to properties phi, mu, q, and ensemble inherited from the Species
   * base class. The size parameter is defined to be the ratio of solvent
   * molecule volume to the monomer reference volume.
   *
   * Each program-level sub-namespace of Pscf defines a subclass of
   * Pscf::SolventSpecies. The Solvent class in each such namespace defines
   * a function that can solve the single-particle statistical mechanics
   * problem for a solvent species.
   *
   * \ingroup Pscf_Chem_Module
   */
   template <typename WT = double>
   class SolventSpecies : public Species<WT>
   {

   public:

      /**
      * Constructor.
      */
      SolventSpecies();

      /**
      * Destructor.
      */
      virtual ~SolventSpecies() = default;

      /**
      * Read parameters and initialize.
      *
      * \param in  input parameter stream
      */
      virtual void readParameters(std::istream& in);

      /// \name Setters (set member data)
      ///@{

      /**
      * Set the monomer id for this solvent.
      *
      * \param monomerId  integer id of monomer type, in [0,nMonomer-1]
      */
      void setMonomerId(int monomerId);

      /**
      * Set the molecular volume of this solvent species.
      *
      * The ``size" is the ratio (solvent molecule volume) / vMonomer,
      * where vMonomer is the monomer reference volume, i.e., the
      * volume per monomer (or unit contour length) of any polymer.
      *
      * \param size  volume of solvent
      */
      void setSize(double size);

      ///@}
      /// \name Accessors (getters)
      ///@{

      /**
      * Get the monomer type id.
      */
      int monomerId() const;

      /**
      * Get the size (number of monomers) in this solvent.
      */
      double size() const;

      ///@}

   private:

      /// Identifier for the associated monomer type.
      int monomerId_;

      /// Size of this block = volume / monomer reference volume.
      double size_;

   };

   // Inline member functions

   /*
   * Get the monomer type id.
   */
   template <typename WT> inline 
   int SolventSpecies<WT>::monomerId() const
   {  return monomerId_; }

   /*
   * Get the size (number of monomers) in this block.
   */
   template <typename WT> inline 
   double SolventSpecies<WT>::size() const
   {  return size_; }

   // Explicit template instantiation declaration
   extern template class SolventSpecies<double>;
}
#endif

#ifndef PSCF_SOLVENT_DESCRIPTOR_H
#define PSCF_SOLVENT_DESCRIPTOR_H

/*
* PSCF - Polymer Self-Consistent Field Theory
*
* Copyright 2016 - 2022, The Regents of the University of Minnesota
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
   * \ingroup Pscf_Chem_Module
   */
   class SolventSpecies : public Species
   {

   public:

      /**
      * Constructor.
      */
      SolventSpecies();

      /**
      * Constructor.
      */
      ~SolventSpecies();

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

      // Inherited accessor functions
      using Pscf::Species::phi;
      using Pscf::Species::mu;
      using Pscf::Species::q;
      using Pscf::Species::ensemble;

   protected:

      /// Identifier for the associated monomer type.
      int monomerId_;

      /// Size of this block = volume / monomer reference volume.
      double size_;

   };

   // Inline member functions

   /*
   * Get the monomer type id.
   */
   inline int SolventSpecies::monomerId() const
   {  return monomerId_; }

   /*
   * Get the size (number of monomers) in this block.
   */
   inline double SolventSpecies::size() const
   {  return size_; }

}
#endif

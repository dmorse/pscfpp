#ifndef PSCF_FH_MOLECULE_H
#define PSCF_FH_MOLECULE_H

/*
* PSCF - Polymer Self-Consistent Field
*
* Copyright 2015 - 2026, The Regents of the University of Minnesota
* Distributed under the terms of the GNU General Public License.
*/

#include <util/param/ParamComposite.h>   // base class

#include <pscf/floryHuggins/FhClump.h>   // member
#include <util/containers/Pair.h>        // member
#include <util/containers/DArray.h>      // member

namespace Pscf {

   using namespace Util;

   /**
   * Molecular species in a homogeneous Flory-Huggins mixture.
   *
   * A FhMolecule has:
   *
   *  - An array of FhClump objects
   *  - An overall size (volume/monomer volume)
   *
   * Each FhClump has a monomer type id and a size. The size is the
   * total volume of monomers of that type in a molecule of this
   * species.
   *
   * \ingroup Pscf_FloryHuggins_Module
   */
   class FhMolecule : public Util::ParamComposite
   {

   public:

      /**
      * Constructor.
      */
      FhMolecule();

      /**
      * Destructor.
      */
      ~FhMolecule();

      /**
      * Read and initialize.
      *
      * Call either this or setNClump to initialize, not both.
      *
      * \param in input parameter stream
      */
      virtual void readParameters(std::istream& in);

      /**
      * Set the number of clumps, and allocate memory.
      *
      * Call either this or readParameters to initialize, not both.
      * If this is used to allocate memory, all clump properties must
      * be set using FhClump::setMonomerId() and FhClump::setSize().
      */
      void setNClump(int nClump);

      /**
      * Compute total molecule size by adding clump sizes.
      */
      void computeSize();

      /// \name Accessors
      //@{

      /**
      * Get a specified FhClump.
      *
      * \param id clump index, 0 <= id < nClump
      */
      FhClump& clump(int id);

      /**
      * Get a specified FhClump.
      *
      * \param id clump index, 0 <= id < nClump
      */
      FhClump const & clump(int id) const;

      /**
      * Number of monomer clumps (monomer types).
      */
      int nClump() const;

      /**
      * Total molecule size  = volume / reference volume.
      */
      double size() const;

      //@}

   private:

      /// Array of FhClump objects in this polymer.
      DArray<FhClump> clumps_;

      /// Number of clumps in this polymer
      int nClump_;

      /// Total size of all clumps (in units of reference size).
      double size_;

      /// Flag set when computeSize is called.
      bool hasSize_;

   };

   /*
   * Number of clumps.
   */
   inline 
   int FhMolecule::nClump() const
   {  return nClump_; }

   /*
   * Total size of all clumps = volume / reference volume
   */
   inline 
   double FhMolecule::size() const
   {
      UTIL_CHECK(hasSize_);
      return size_;
   }

   /*
   * Get a specified FhClump (non-constant reference)
   */
   inline 
   FhClump& FhMolecule::clump(int id)
   {  return clumps_[id]; }

   /*
   * Get a specified FhClump (constant reference)
   */
   inline 
   const FhClump& FhMolecule::clump(int id) const
   {  return clumps_[id]; }

}
#endif

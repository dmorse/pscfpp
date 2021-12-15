#ifndef PSCF_HOMOGENEOUS_MOLECULE_H
#define PSCF_HOMOGENEOUS_MOLECULE_H

/*
* PSCF - Polymer Self-Consistent Field Theory
*
* Copyright 2016 - 2019, The Regents of the University of Minnesota
* Distributed under the terms of the GNU General Public License.
*/

#include <util/param/ParamComposite.h>   // base class

#include <pscf/homogeneous/Clump.h>      // member template argument
#include <util/containers/Pair.h>        // member template
#include <util/containers/DArray.h>      // member template

#include <cmath>

namespace Pscf { 
   namespace Homogeneous { 
   
      using namespace Util;
   
      /**
      * Descriptor of a molecular species in a homogeneous mixture.
      *
      * A Homogeneous::Molecule has:
      *
      *  - An array of Homogeneous::Clump objects
      *  - An overall size (volume/monomer volume)
      *
      * Each Clump has a monomer type id and a size. The size is the
      * total volume of monomers of that type in a molecule of this
      * species.
      *
      * \ingroup Pscf_Homogeneous_Module
      */
      class Molecule : public Util::ParamComposite
      {
   
      public:
   
         /**
         * Constructor.
         */
         Molecule();
    
         /**
         * Destructor.
         */
         ~Molecule();
   
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
         * If this is used to allocate memory, all clump properties
         * must be set using Clump::setMonomerId() and Clump::setSize().
         */
         void setNClump(int nClump);
   
         /**
         * Compute total molecule size by adding clump sizes.
         */
         void computeSize();
   
         /// \name Accessors 
         //@{
   
         /**
         * Get a specified Clump.
         *
         * \param id clump index, 0 <= id < nClump
         */
         Clump& clump(int id);
   
         /**
         * Get a specified Clump.
         *
         * \param id clump index, 0 <= id < nClump
         */
         Clump const & clump(int id) const;
   
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
   
         /// Array of Clump objects in this polymer.
         DArray<Clump> clumps_;
   
         /// Number of clumps in this polymer
         int nClump_;
   
         /// Total size of all clumps (in units of reference size).
         double size_;

         /// Flag set when computeSize is called.
         bool hasSize_;
   
      };
   
   }

   /*
   * Number of clumps.
   */
   inline int Homogeneous::Molecule::nClump() const
   {  return nClump_; }

   /*
   * Total size of all clumps = volume / reference volume
   */
   inline double Homogeneous::Molecule::size() const
   {
      UTIL_CHECK(hasSize_);
      return size_; 
   }

   /*
   * Get a specified Clump (non-constant reference)
   */
   inline Homogeneous::Clump& Homogeneous::Molecule::clump(int id)
   {  return clumps_[id]; }
 
   /*
   * Get a specified Clump (constant reference)
   */
   inline const Homogeneous::Clump& Homogeneous::Molecule::clump(int id) const
   {  return clumps_[id]; }
 
}
#endif

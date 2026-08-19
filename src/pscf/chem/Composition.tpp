/*
* PSCF - Polymer Self-Consistent Field
*
* Copyright 2015 - 2026, The Regents of the University of Minnesota
* Distributed under the terms of the GNU General Public License.
*/

#include "Composition.h"
#include "PolymerSpecies.h"
#include "SolventSpecies.h"

namespace Pscf
{

   /*
   * Constructor.
   */
   template <typename WT>
   Composition<WT>::Composition()
    : monomers_(),
      nMonomer_(0),
      nPolymer_(0),
      nSolvent_(0),
      nBlock_(0),
      vMonomer_(1.0)
   {}

   /*
   * Set the monomer volume value.
   */
   template <typename WT>
   void Composition<WT>::setVmonomer(double vMonomer)
   {
      UTIL_CHECK(vMonomer > 0.0);  
      vMonomer_ = vMonomer; 
   }

   /*
   * Is the ensemble closed for all species?
   */
   template <typename WT>
   bool Composition<WT>::isCanonical() const
   {
      // Check ensemble of all polymer species
      for (int i = 0; i < nPolymer(); ++i) {
         if (polymerSpecies(i).ensemble() == Ensemble::Open) {
            return false;
         }
      }

      // Check ensemble of all solvent species
      for (int i = 0; i < nSolvent(); ++i) {
         if (solventSpecies(i).ensemble() == Ensemble::Open) {
            return false;
         }
      }

      // Returns true if false was not returned earlier
      return true;
   }

}

/*
* PSCF - Polymer Self-Consistent Field
*
* Copyright 2015 - 2025, The Regents of the University of Minnesota
* Distributed under the terms of the GNU General Public License.
*/

#include "MixtureBase.h"
#include "PolymerSpecies.h"
#include "SolventSpecies.h"

namespace Pscf
{

   /*
   * Constructor.
   */
   MixtureBase::MixtureBase()
    : monomers_(),
      nMonomer_(0),
      nPolymer_(0),
      nSolvent_(0),
      nBlock_(0),
      vMonomer_(1.0)
   {}

   /*
   * Destructor.
   */
   MixtureBase::~MixtureBase()
   {}

   void MixtureBase::setVmonomer(double vMonomer)
   {
      UTIL_CHECK(vMonomer > 0.0);  
      vMonomer_ = vMonomer; 
   }

   /*
   * Is the ensemble closed for every species.
   */
   bool MixtureBase::isCanonical() const
   {
      // Check ensemble of all polymer species
      for (int i = 0; i < nPolymer(); ++i) {
         if (polymerSpecies(i).ensemble() == Species::Open) {
            return false;
         }
      }

      // Check ensemble of all solvent species
      for (int i = 0; i < nSolvent(); ++i) {
         if (solventSpecies(i).ensemble() == Species::Open) {
            return false;
         }
      }

      // Returns true if false was not returned earlier
      return true;
   }

}

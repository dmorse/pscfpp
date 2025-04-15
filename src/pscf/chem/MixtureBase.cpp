/*
* PSCF - Polymer Self-Consistent Field Theory
*
* Copyright 2016 - 2025, The Regents of the University of Minnesota
* Distributed under the terms of the GNU General Public License.
*/

#include "MixtureBase.h"

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

}

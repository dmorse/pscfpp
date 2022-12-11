/*
* PSCF - Polymer Self-Consistent Field Theory
*
* Copyright 2016 - 2022, The Regents of the University of Minnesota
* Distributed under the terms of the GNU General Public License.
*/

#include "AmbdInteraction.h"
#include <pscf/math/LuSolver.h>

namespace Pscf 
{

   using namespace Util;

   /*
   * Constructor.
   */
   AmbdInteraction::AmbdInteraction()
    : nMonomer_(0),
      isAllocated_(false)
   {}

   /*
   * Destructor.
   */
   AmbdInteraction::~AmbdInteraction()
   {}

   /*
   * Set the number of monomer types and allocate memory.
   */
   void AmbdInteraction::setNMonomer(int nMonomer)
   {  
      UTIL_CHECK(isAllocated_ == false);
      UTIL_CHECK(nMonomer_ == 0);
      UTIL_CHECK(nMonomer > 0);
      nMonomer_ = nMonomer; 
      chi_.allocate(nMonomer, nMonomer);
      chiInverse_.allocate(nMonomer, nMonomer);
      p_.allocate(nMonomer, nMonomer);
      isAllocated_ = true;
   }

   void AmbdInteraction::update(Interaction const & interaction)
   {

      // Set nMonomer and allocate memory if not done previously
      if (!isAllocated_) {
         setNMonomer(interaction.nMonomer());
      }
      UTIL_CHECK(nMonomer_ == interaction.nMonomer());

      // Copy all chi and chiInverse matrix values
      int i, j;
      for (i = 0; i < nMonomer_; ++i) {
         for (j = 0; j < nMonomer_; ++j) {
           chi_(i, j) = interaction.chi(i, j);
           chiInverse_(i, j) = interaction.chiInverse(i, j);
         }
      }

      // Compute p and sumChiInverse
      int k;
      double sum = 0.0;
      for (i = 0; i < nMonomer_; ++i) {
         p_(0,i) = 0.0;
         for (j = 0; j < nMonomer_; ++j) {
            p_(0,i) -= chiInverse_(j,i);
         }
         sum -= p_(0,i);
         for (k = 0; k < nMonomer_; ++k) { //row
            p_(k,i) = p_(0,i);
         }
      }

      for (i = 0; i < nMonomer_; ++i) { //row
         for (j = 0; j < nMonomer_; ++j) { //coloumn
            p_(i,j) /= sum;
         }
         p_(i,i) += 1.0 ;
      }
      sumChiInverse_ = sum;

   }

} // namespace Pscf

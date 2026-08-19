/*
* PSCF - Polymer Self-Consistent Field
*
* Copyright 2015 - 2026, The Regents of the University of Minnesota
* Distributed under the terms of the GNU General Public License.
*/

#include "AmbdInteraction.h"
#include <pscf/interaction/Interaction.h>
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
           //chiInverse_(i, j) = interaction.chiInverse(i, j);
         }
      }

      if (nMonomer() == 2) {
         double det = chi_(0,0)*chi_(1, 1) - chi_(0,1)*chi_(1,0);
         double norm = chi_(0,0)*chi_(0, 0) + chi_(1,1)*chi_(1,1)
                     + 2.0*chi_(0,1)*chi_(1,0);
         if (fabs(det/norm) < 1.0E-8) {
            UTIL_THROW("Singular chi matrix");
         }
         chiInverse_(0,1) = -chi_(0,1)/det;
         chiInverse_(1,0) = -chi_(1,0)/det;
         chiInverse_(1,1) = chi_(0,0)/det;
         chiInverse_(0,0) = chi_(1,1)/det;

      } else {
         LuSolver solver;
         solver.allocate(nMonomer());
         solver.computeLU(chi_);
         solver.inverse(chiInverse_);
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

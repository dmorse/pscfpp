/*
* PSCF - Polymer Self-Consistent Field Theory
*
* Copyright 2016 - 2022, The Regents of the University of Minnesota
* Distributed under the terms of the GNU General Public License.
*/

#include "AmbdResidual.h"
#include <pscf/math/LuSolver.h>

namespace Pscf 
{

   using namespace Util;

   /*
   * Constructor.
   */
   AmbdResidual::AmbdResidual()
    : nMonomer_(0),
      isAllocated_(false)
   {}

   /*
   * Destructor.
   */
   AmbdResidual::~AmbdResidual()
   {}

   /*
   * Set the number of monomer types and allocate memory.
   */
   void AmbdResidual::allocate(int nMonomer)
   {  
      UTIL_CHECK(isAllocated_ == false);
      nMonomer_ = nMonomer; 
      chi_.allocate(nMonomer, nMonomer);
      chiInverse_.allocate(nMonomer, nMonomer);
      idemp_.allocate(nMonomer, nMonomer);
      isAllocated_ = true;
   }

   void AmbdResidual::update(Interaction const & interaction)
   {

      // Allocate memory if not done previously
      if (!isAllocated_) {
         allocate(interaction.nMonomer());
      }

      int i, j, k;

      // Copy all chi matrix values
      for (i = 0; i < nMonomer_; ++i) {
         for (j = 0; j < nMonomer_; ++j) {
           chi_(i, j) = interaction.chi(i, j);
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

      double sum = 0;
      for (i = 0; i < nMonomer(); ++i) {
         idemp_(0,i) = 0;
         for (j = 0; j < nMonomer(); ++j) {
            idemp_(0,i) -= chiInverse_(j,i);
         }
         sum -= idemp_(0,i);
         for (k = 0; k < nMonomer(); ++k) { //row
            idemp_(k,i) = idemp_(0,i);
         }
      }

      for (i = 0; i < nMonomer(); ++i) { //row
         for (j = 0; j < nMonomer(); ++j) { //coloumn
            idemp_(i,j) /= sum;
         }
         idemp_(i,i) += 1.0 ;
      }
      
      sum_inv_ = sum;

   }

} // namespace Pscf

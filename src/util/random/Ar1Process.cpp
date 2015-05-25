/*
* Simpatico - Simulation Package for Polymeric and Molecular Liquids
*
* Copyright 2010 - 2014, The Regents of the University of Minnesota
* Distributed under the terms of the GNU General Public License.
*/

#include "Ar1Process.h"

namespace Util
{

   /*
   * Default constructor.
   */   
   Ar1Process::Ar1Process()
    : x_(0.0),
      B_(1.0),
      C_(0.0),
      randomPtr_(0),
      isInitialized_(false)
   {}

   /*
   * Constructor.
   */   
   Ar1Process::Ar1Process(Random& random)
    : x_(0.0),
      B_(1.0),
      C_(0.0),
      randomPtr_(&random),
      isInitialized_(false)
   {}

   /*
   * Set pointer to random number generator.
   */
   void Ar1Process::setRNG(Random& random)
   { randomPtr_ = &random; }

   /*
   * Initialize process
   *
   * \param tau decay time (in discrete steps)
   */
   void Ar1Process::init(double tau)
   { 
      // Precondition
      if (randomPtr_ == 0) {
         UTIL_THROW("Random number generator not yet set");
      }

      C_ = exp(-1.0/tau);
      B_ = sqrt(1.0 - C_*C_);
      x_ = B_*randomPtr_->gaussian();
      isInitialized_ = true;
   }

} 

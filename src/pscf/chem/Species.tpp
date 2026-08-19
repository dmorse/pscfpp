/*
* PSCF - Polymer Self-Consistent Field
*
* Copyright 2015 - 2026, The Regents of the University of Minnesota
* Distributed under the terms of the GNU General Public License.
*/

#include "Species.h"
#include <pscf/math/arithmetic.h>
#include <cmath>

namespace Pscf { 

   using namespace Util;

   /*
   * Default constructor.
   */
   template <typename WT>
   Species<WT>::Species()
    : phiMu_(0.0),
      ensemble_(Ensemble::Closed)
   { 
      assign(phi_, 0.0); 
      assign(mu_, 0.0); 
      assign(q_, 0.0); 
      setClassName("Species"); 
   }

   /*
   * Read phi or mu (but not both) and set ensemble.
   */
   template <typename WT>
   void Species<WT>::readParameters(std::istream& in)
   {
      // Read phi or mu (but not both)
      bool hasPhi = readOptional(in, "phi", phiMu_).isActive();
      if (hasPhi) {
         ensemble_ = Ensemble::Closed;
         assign(phi_, phiMu_);
      } else {
         ensemble_ = Ensemble::Open;
         read(in, "mu", phiMu_);
         assign(mu_, phiMu_);
      }
   }

   /*
   *  Set volume fraction (if ensemble is closed).
   */ 
   template <typename WT>
   void Species<WT>::setPhi(double phi)
   {
      UTIL_CHECK(ensemble() == Ensemble::Closed);  
      UTIL_CHECK(phi >= 0.0);  
      UTIL_CHECK(phi <= 1.0);  
      assign(phiMu_, phi);
      assign(phi_, phi);
   }

   /*
   *  Set chemical potential (if ensemble is open).
   */ 
   template <typename WT>
   void Species<WT>::setMu(double mu)
   {
      UTIL_CHECK(ensemble() == Ensemble::Open);
      assign(phiMu_, mu);
      assign(mu_, mu);
   }

   /*
   * Set q and compute mu or phi (depending on ensemble).
   */
   template <typename WT>
   void Species<WT>::setQ(WT q)
   {
      assign(q_, q);
      if (ensemble() == Ensemble::Closed) {
         WT ratio;
         div(ratio, phi_, q_);
         assignLog(mu_, ratio);
         //mu_ = std::log(phi_/q_);
      } else
      if (ensemble() == Ensemble::Open) {
         WT lambda;
         assignExp(lambda, mu_);
         mul(phi_, lambda, q_);
         //phi_ = std::exp(mu_)*q_;
      }
   }

} // namespace Pscf

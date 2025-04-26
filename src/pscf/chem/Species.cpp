/*
* PSCF - Polymer Self-Consistent Field Theory
*
* Copyright 2015 - 2025, The Regents of the University of Minnesota
* Distributed under the terms of the GNU General Public License.
*/

#include "Species.h"

namespace Pscf
{ 

   using namespace Util;

   /*
   * Default constructor.
   */
   Species::Species()
    : phi_(0.0),
      mu_(0.0),
      q_(0.0),
      ensemble_(Species::Closed)
   {  setClassName("Species"); }

   /*
   * Read phi or mu (but not both) and set ensemble.
   */
   void Species::readParameters(std::istream& in)
   {
      // Read phi or mu (but not both)
      bool hasPhi = readOptional(in, "phi", phi_).isActive();
      if (hasPhi) {
         ensemble_ = Species::Closed;
      } else {
         ensemble_ = Species::Open;
         read(in, "mu", mu_);
      }
   }

   /*
   *  Set volume fraction (if ensemble is closed).
   */ 
   void Species::setPhi(double phi)
   {
      UTIL_CHECK(ensemble() == Species::Closed);  
      UTIL_CHECK(phi >= 0.0);  
      UTIL_CHECK(phi <= 1.0);  
      phi_ = phi;
   }

   /*
   *  Set chemical potential (if ensemble is open).
   */ 
   void Species::setMu(double mu)
   {
      UTIL_CHECK(ensemble() == Species::Open);  
      mu_ = mu;
   }

   /*
   * Set q and compute mu or phi (depending on ensemble).
   */
   void Species::setQ(double q)
   {
      q_ = q;
      if (ensemble() == Species::Closed) {
         mu_ = log(phi_/q_);
      } else
      if (ensemble() == Species::Open) {
         phi_ = exp(mu_)*q_;
      }
   }

   /* 
   * Extract a Species::Ensemble from an istream as a string.
   */
   std::istream& operator >> (std::istream& in, Species::Ensemble& policy)
   {
      std::string buffer;
      in >> buffer;
      if (buffer == "Closed" || buffer == "closed") {
         policy = Species::Closed;
      } else 
      if (buffer == "Open" || buffer == "open") {
         policy = Species::Open;
      } else {
         UTIL_THROW("Invalid Species::Ensemble string in operator >>");
      } 
      return in;
   }
   
   /* 
   * Insert a Species::Ensemble to an ostream as a string.
   */
   std::ostream& operator<<(std::ostream& out, Species::Ensemble policy) 
   {
      if (policy == Species::Closed) {
         out << "Closed";
      } else 
      if (policy == Species::Open) {
         out << "Open";
      } else 
      if (policy == Species::Unknown) {
         out << "Unknown";
      } else {
         std::cout << "Invalid Species::Ensemble value on input" 
                   << std::endl;
         UTIL_THROW("Unrecognized value for Species::Ensemble");
      } 
      return out; 
   }

} // namespace Pscf

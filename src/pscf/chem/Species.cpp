/*
* PSCF - Polymer Self-Consistent Field Theory
*
* Copyright 2016 - 2022, The Regents of the University of Minnesota
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
   {}

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
         std::cout << "Invalid Species::Ensemble value on input" << std::endl;
         UTIL_THROW("Unrecognized value for Species::Ensemble");
      } 
      return out; 
   }

} // namespace Pscf

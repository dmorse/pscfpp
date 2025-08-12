/*
* PSCF - Polymer Self-Consistent Field
*
* Copyright 2015 - 2025, The Regents of the University of Minnesota
* Distributed under the terms of the GNU General Public License.
*/

#include "Environment.h"

namespace Pscf {
namespace Prdc {

   // Constructor
   Environment::Environment()
    : hasStress_(false),
      nParam_(0)
   {  setClassName("Environment"); }

   // Destructor
   Environment::~Environment()
   {}

   // Sets needsUpdate() to true and hasStress() to false. 
   void Environment::reset()
   {
      EnvironmentBase::reset();
      hasStress_ = false;
      stress_.clear();
      stressIds_.clear();
   }

   // Modify stress to account for Environment for one lattice parameter
   double Environment::modifyStress(int paramId, double stress) const
   {  UTIL_THROW("Unimplemented modifyStress() method called."); }

   // Return the environment-modified SCFT stress for one parameter.
   double Environment::stress(int paramId) const
   {
      UTIL_CHECK(hasStress_);
      for (int i = 0; i < stress_.size(); i++) {
         if (stressIds_[i] == paramId) {
            return stress_[i];
         }
      }
      
      // If this point is reached, stress not found
      UTIL_THROW("Attempt to access stress value that was not calculated.");
   }

   // Write environment-modified stress to output stream.
   void Environment::writeStress(std::ostream& out) const
   {
      UTIL_CHECK(hasStress_);
      UTIL_CHECK(nParam_ > 0);

      // If no stresses have been set, do nothing and return
      if (stress_.size() == 0) return;

      // Write only the modified stresses to the ostream
      out << "environment-modified stress:" << std::endl;
      for (int i = 0; i < stress_.size(); i++) {
         out << Int(stressIds_[i], 5)
            << "  "
            << Dbl(stress_[i], 18, 11)
            << std::endl;
      }
      out << std::endl;
   }

   // Has the stress been calculated?
   bool Environment::hasStress() const
   {  return hasStress_; }

   // Set the number of lattice parameters
   void Environment::setNParams(int nParams)
   {
      UTIL_CHECK(nParam_ == 0);
      nParam_ = nParams;
   }

} // namespace Prdc
} // namespace Pscf
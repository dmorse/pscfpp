#ifndef PRDC_FILM_FIELD_GEN_EXT_BASE_TPP
#define PRDC_FILM_FIELD_GEN_EXT_BASE_TPP

/*
* PSCF - Polymer Self-Consistent Field
*
* Copyright 2015 - 2025, The Regents of the University of Minnesota
* Distributed under the terms of the GNU General Public License.
*/

#include "FilmFieldGenExtBase.h"
#include <prdc/crystal/SpaceGroup.h>
#include <pscf/sweep/ParameterType.h>
#include <util/param/ParamComponent.h>
#include <util/containers/GArray.h>

namespace Pscf {
namespace Prdc {

   using namespace Util;

   /*
   * Constructor
   */
   template <int D>
   FilmFieldGenExtBase<D>::FilmFieldGenExtBase()
    : FieldGenerator::FieldGenerator(),
      normalVecId_(-1),
      interfaceThickness_(-1.0),
      excludedThickness_(-1.0)
   {
      type_ = External; 
      isDependent_ = true;
   }

   /*
   * Destructor
   */
   template <int D>
   FilmFieldGenExtBase<D>::~FilmFieldGenExtBase()
   {}

   /*
   * Read and initialize.
   */
   template <int D>
   void FilmFieldGenExtBase<D>::readParameters(std::istream& in)
   {
      // First, read data defining the mask (quietly, with echo = false)
      bool echo = ParamComponent::echo();
      ParamComponent::setEcho(false);
      read(in, "normalVecId", normalVecId_);
      read(in, "interfaceThickness", interfaceThickness_);
      read(in, "excludedThickness", excludedThickness_);
      double tmp; 
      readOptional(in, "fBulk", tmp); // will not be used
      ParamComponent::setEcho(echo);

      // Remove all of these parameters from this ParamComposite, since
      // they are already managed by the FilmFieldGenMask ParamComposite
      ParamComposite::resetParam();

      // Make sure inputs are valid
      if (normalVecId_ > D || normalVecId_ < 0) {
         UTIL_THROW("bad value for normalVecId, must be in [0,D)");
      }
      if (interfaceThickness_ > excludedThickness_) {
         UTIL_THROW("excludedThickness must be larger than interfaceThickness");
      }
      if ((excludedThickness_ <= 0) || (interfaceThickness_ <= 0)) {
         UTIL_THROW("excludedThickness and interfaceThickness must be >0");
      }

      // Allocate chiBottom_ and chiTop_
      int nm = systemNMonomer();
      chiBottom_.allocate(nm);
      chiTop_.allocate(nm);

      // Read chiBottom_ and chiTop_ arrays
      readDArray(in, "chiBottom", chiBottom_, nm);
      readDArray(in, "chiTop", chiTop_, nm);
   }

   /*
   * Check whether system has changed such that the fields need updating
   */
   template <int D>
   bool FilmFieldGenExtBase<D>::needsUpdate() const
   {
      UTIL_CHECK(normalVecId_ >= 0); // Check that readParameters was called

      // If chiBottomCurrent_ and chiTopCurrent_ are unset, generate() has 
      // not been called. Therefore, needsUpdate is true
      if (!chiBottomCurrent_.isAllocated()) {
         UTIL_CHECK(!chiTopCurrent_.isAllocated());
         return true;
      }
      
      // Check if chiBottom or chiTop have been changed
      for (int i = 0; i < chiBottom_.capacity(); i++) {
         if ((chiBottom_[i] != chiBottomCurrent_[i]) || 
               (chiTop_[i] != chiTopCurrent_[i])) {
            return true;
         }
      }

      // If chiTop and chiBottom are unchanged and all zeros, no update needed
      if (isAthermal()) return false;

      // Check if system normalVec differ from normalVecCurrent_
      if (normalVecCurrent_ == systemLatticeVector(normalVecId_)) {
         return false;
      } else {
         return true;
      }
   }

   /*
   * Check that the system is compatible with this field
   */
   template <int D>
   void FilmFieldGenExtBase<D>::checkCompatibility()
   {
      // Ensure that space group symmetry is compatible with the fields

      // If chiBottom == chiTop, do nothing. All necessary checks have 
      // already been performed by the mask generator.
      if (hasSymmetricWalls()) return;

      // Otherwise, walls are asymmetric, so the space group must not 
      // have any operations that swap the two walls
      std::string groupName = systemSpaceGroup();
      SpaceGroup<D> group;
      std::ifstream in;

      // Open and read file containing space group's symmetry operations
      readGroup(groupName, group);

      // Make sure all symmetry operations are allowed
      std::string msg = "Space group contains forbidden symmetry operations";
      for (int i = 0; i < group.size(); i++) {
         for (int j = 0; j < D; j++) {
            int r = group[i].R(normalVecId_,j);
            if (j == normalVecId_) {
               if (r != 1) {
                  UTIL_THROW(msg.c_str());
               }
            } else { // j != normalVecId_
               if (r != 0) {
                  UTIL_THROW(msg.c_str());
               }
            }
         }
         if (group[i].t(normalVecId_) != 0) {
            UTIL_THROW(msg.c_str());
         }
      }
   }

   /*
   * Check whether or not the two walls are chemically identical using
   * the chi array.
   */
   template <int D>
   bool FilmFieldGenExtBase<D>::hasSymmetricWalls() const 
   {
      int nm = systemNMonomer();

      UTIL_CHECK(nm > 0);
      UTIL_CHECK(chiBottom_.capacity() == nm);
      UTIL_CHECK(chiTop_.capacity() == nm);

      for (int i = 0; i < nm; i++) {
         if (fabs(chiBottom_[i]-chiTop_[i]) > 1e-7) {
            return false;
         }
      }
      return true;
   }

   /*
   * Check whether or not the walls are athermal (only true if all values
   * in chi array are zero)
   */
   template <int D>
   bool FilmFieldGenExtBase<D>::isAthermal() const 
   {
      int nm = systemNMonomer();

      UTIL_CHECK(nm > 0);
      UTIL_CHECK(chiBottom_.capacity() == nm);
      UTIL_CHECK(chiTop_.capacity() == nm);

      for (int i = 0; i < nm; i++) {
         if ((fabs(chiBottom_[i]) >= 1e-7) || (fabs(chiTop_[i]) >= 1e-7)) {
            return false;
         }
      }
      return true;
   }

   /*
   * Return specialized sweep parameter types to add to the Sweep object.
   */
   template <int D>
   GArray<ParameterType> 
   FilmFieldGenExtBase<D>::getParameterTypes()
   {
      GArray<ParameterType> pTypes;
      pTypes.append(ParameterType("chi_top", 1, *this));
      pTypes.append(ParameterType("chi_bottom", 1, *this));
      return pTypes;
   }

   /*
   * Set the value of a specialized sweep parameter (chi_top or chi_bottom).
   */
   template <int D>
   void FilmFieldGenExtBase<D>::setParameter(std::string name, 
                                             DArray<int> ids, 
                                             double value,
                                             bool& success)
   {
      success = true;
      if (name == "chi_top") {
         chiTop_[ids[0]] = value;
      } else if (name == "chi_bottom") {
         chiBottom_[ids[0]] = value;
      } else {
         success = false;
      }
   }

   /*
   * Get the value of a specialized sweep parameter (chi_top or chi_bottom).
   */
   template <int D>
   double FilmFieldGenExtBase<D>::getParameter(std::string name, 
                                               DArray<int> ids, 
                                               bool& success) 
   const
   {
      success = true;
      if (name == "chi_top") {
         return chiTop_[ids[0]];
      } else if (name == "chi_bottom") {
         return chiBottom_[ids[0]];
      } else {
         success = false;
         return 0.0;
      }
   }

}
}
#endif

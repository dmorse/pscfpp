#ifndef PRDC_EXT_GEN_FILM_BASE_TPP
#define PRDC_EXT_GEN_FILM_BASE_TPP

/*
* PSCF - Polymer Self-Consistent Field Theory
*
* Copyright 2016 - 2022, The Regents of the University of Minnesota
* Distributed under the terms of the GNU General Public License.
*/

#include "ExtGenFilmBase.h"
#include "prdc/crystal/SpaceGroup.h"

namespace Pscf {
namespace Prdc
{

   using namespace Util;

   /*
   * Constructor
   */
   template <int D>
   ExtGenFilmBase<D>::ExtGenFilmBase()
    : FieldGenerator::FieldGenerator(),
      normalVecCurrent_(),
      chiBottomCurrent_(),
      chiTopCurrent_(),
      normalVecId_(-1),
      interfaceThickness_(-1.0),
      excludedThickness_(-1.0),
      chiBottom_(),
      chiTop_()
   {
      type_ = External; 
      isDependent_ = true;
   }

   /*
   * Destructor
   */
   template <int D>
   ExtGenFilmBase<D>::~ExtGenFilmBase()
   {}

   /*
   * Read and initialize.
   */
   template <int D>
   void ExtGenFilmBase<D>::readParameters(std::istream& in)
   {
      // First, read data defining the mask
      read(in, "normalVecId", normalVecId_);
      read(in, "interfaceThickness", interfaceThickness_);
      read(in, "excludedThickness", excludedThickness_);
      double fBulk; // will not be used
      readOptional(in, "fBulk", fBulk);

      // Allocate chiBottom_ and chiTop_
      int nm = systemNMonomer();
      chiBottom_.allocate(nm);
      chiTop_.allocate(nm);

      // Read chiBottom_ and chiTop_ arrays
      readDArray(in, "chiBottom", chiBottom_, nm);
      readDArray(in, "chiTop", chiTop_, nm);
   }

   /*
   * Check that the system is compatible with this field
   */
   template <int D>
   void ExtGenFilmBase<D>::checkCompatibility()
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
   * Check whether system has changed such that the field needs updating
   */
   template <int D>
   bool ExtGenFilmBase<D>::updateNeeded() const
   {
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
   * Check whether or not the two walls are chemically identical using
   * the chi array.
   */
   template <int D>
   bool ExtGenFilmBase<D>::hasSymmetricWalls() const 
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
   bool ExtGenFilmBase<D>::isAthermal() const 
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
   ExtGenFilmBase<D>::getParameterTypes()
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
   void ExtGenFilmBase<D>::setParameter(std::string name, 
                                        DArray<int> ids, 
                                        double value,
                                        bool& success)
   {
      success = true;
      bool wasAthermal = isAthermal();
      if (name == "chi_top") {
         chiTop_[ids[0]] = value;
      } else if (name == "chi_bottom") {
         chiBottom_[ids[0]] = value;
      } else {
         success = false;
      }
      
      // If walls were athermal but are not anymore, then the h fields
      // have not been allocated. Thus, we must call allocate() again.
      if (wasAthermal && success && (!isAthermal())) {
         allocate();
      }
   }

   /*
   * Get the value of a specialized sweep parameter (chi_top or chi_bottom).
   */
   template <int D>
   double ExtGenFilmBase<D>::getParameter(std::string name, 
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
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
      parameters_(),
      chiBottomCurrent_(),
      chiTopCurrent_(),
      normalVecId_(-1),
      chiBottom_(),
      chiTop_()
   {  type_ = External; }

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
      // Allocate chiBottom_ and chiTop_ and set to zero before 
      // reading them in
      int nm = getNMonomer();
      chiBottom_.allocate(nm);
      chiTop_.allocate(nm);

      // Read arrays
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
      std::string groupName = getSpaceGroup();
      SpaceGroup<D> group;
      std::ifstream in;

      // Open and read file containing space group's symmetry operations
      readGroup(groupName, group);

      // Get normalVecId from mask
      getNormalVecId();

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

      // Check if system lattice parameters are different than parameters_
      FSArray<double, 6> sysParams = getLatticeParameters();
      UTIL_CHECK(sysParams.size() == parameters_.size());
      for (int i = 0; i < parameters_.size(); i++) {
         if (fabs(sysParams[i] - parameters_[i]) > 1e-10) {
            return true;
         }
      }

      // If this point is reached, no update is needed
      return false;
   }

   /*
   * Check whether or not the two walls are chemically identical using
   * the chi array.
   */
   template <int D>
   bool ExtGenFilmBase<D>::hasSymmetricWalls() const 
   {
      int nm = getNMonomer();

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
      int nm = getNMonomer();

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
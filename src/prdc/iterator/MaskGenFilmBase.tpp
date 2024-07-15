#ifndef PRDC_MASK_GEN_FILM_BASE_TPP
#define PRDC_MASK_GEN_FILM_BASE_TPP

/*
* PSCF - Polymer Self-Consistent Field Theory
*
* Copyright 2016 - 2022, The Regents of the University of Minnesota
* Distributed under the terms of the GNU General Public License.
*/

#include "MaskGenFilmBase.h"
#include "prdc/crystal/SpaceGroup.h"
#include <cmath>

namespace Pscf {
namespace Prdc
{

   using namespace Util;

   /*
   * Constructor
   */
   template <int D>
   MaskGenFilmBase<D>::MaskGenFilmBase()
    : FieldGenerator::FieldGenerator(),
      parameters_(),
      normalVecId_(-1),
      interfaceThickness_(-1.0),
      excludedThickness_(-1.0)
   {  type_ = Mask; }

   /*
   * Destructor
   */
   template <int D>
   MaskGenFilmBase<D>::~MaskGenFilmBase()
   {}

   /*
   * Read and initialize.
   */
   template <int D>
   void MaskGenFilmBase<D>::readParameters(std::istream& in)
   {
      // Read required data defining the walls
      read(in, "normalVecId", normalVecId_);
      read(in, "interfaceThickness", interfaceThickness_);
      read(in, "excludedThickness", excludedThickness_);

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
   }

   /*
   * Check that the system is compatible with this field
   */
   template <int D>
   void MaskGenFilmBase<D>::checkCompatibility()
   {
      // If lattice parameters are flexible, determine which parameters
      // are allowed to vary, store them in this object, and pass them
      // into the Iterator. The flexibleParams_ member of the Iterator
      // should always be matched to that of this class.
      setFlexibleParams();

      // Ensure that space group symmetry is compatible with the wall
      checkSpaceGroup();

      // Ensure that unit cell is compatible with wall
      checkLatticeVectors();
   }

   /*
   * Check whether system has changed such that the field needs updating
   */
   template <int D>
   bool MaskGenFilmBase<D>::updateNeeded() const
   {
      // Check if system lattice parameters are different than parameters_
      FSArray<double, 6> sysParams = getLatticeParameters();
      UTIL_CHECK(sysParams.size() == parameters_.size());
      bool identical = true; // are the two arrays identical?
      for (int i = 0; i < parameters_.size(); i++) {
         if (fabs(sysParams[i] - parameters_[i]) > 1e-10) {
            identical = false;
            break;
         }
      }

      return identical;
   }

   /*
   * Check that space group is compatible with the mask.
   */
   template <int D>
   void MaskGenFilmBase<D>::checkSpaceGroup() const
   {
      // Setup
      std::string groupName = getSpaceGroup();
      SpaceGroup<D> group;
      std::ifstream in;

      // Open and read file containing space group's symmetry operations
      readGroup(groupName, group);

      // Make sure all symmetry operations are allowed
      int nv = normalVecId();
      std::string msg = "Space group contains forbidden symmetry operations";
      for (int i = 0; i < group.size(); i++) {
         for (int j = 0; j < D; j++) {
            int r = group[i].R(nv,j);
            if (j == nv) {
               if ((r != 1) && (r != -1)) {
                  UTIL_THROW(msg.c_str());
               }
            } else { // j != nv
               if (r != 0) {
                  UTIL_THROW(msg.c_str());
               }
            }
         }
         if (group[i].t(nv) != 0) {
            UTIL_THROW(msg.c_str());
         }
      }
   }

}
}
#endif
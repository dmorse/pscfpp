#ifndef PRDC_MASK_GEN_FILM_BASE_TPP
#define PRDC_MASK_GEN_FILM_BASE_TPP

/*
* PSCF - Polymer Self-Consistent Field Theory
*
* Copyright 2015 - 2025, The Regents of the University of Minnesota
* Distributed under the terms of the GNU General Public License.
*/

#include "MaskGenFilmBase.h"
#include "prdc/crystal/SpaceGroup.h"
#include "util/param/ScalarParam.h"
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
      normalVecCurrent_(),
      fBulk_(),
      normalVecId_(-1),
      interfaceThickness_(-1.0),
      excludedThickness_(-1.0),
      hasFBulk_(false)
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
      ScalarParam<double>& fBulkParam = readOptional(in, "fBulk", fBulk_);
      if (fBulkParam.isActive()) { // if we read fBulk
         hasFBulk_ = true;
      }

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
      UTIL_CHECK(isGenerated());
      UTIL_CHECK(normalVecId_ >= 0);
      
      // Check if system normalVec differ from normalVecCurrent_
      if (normalVecCurrent_ == systemLatticeVector(normalVecId_)) {
         return false;
      } else {
         return true;
      }
   }

   /*
   * Check that space group is compatible with the mask.
   */
   template <int D>
   void MaskGenFilmBase<D>::checkSpaceGroup() const
   {
      // Setup
      std::string groupName = systemSpaceGroup();
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

   // Explicit Specializations for checkLatticeVectors are in
   // MaskGenFilmBase.cpp
}
}
#endif
#ifndef PRDC_HAS_VARIABLE_ANGLE_H
#define PRDC_HAS_VARIABLE_ANGLE_H

#include <prdc/crystal/UnitCell.h>

/*
* PSCF - Polymer Self-Consistent Field Theory
*
* Copyright 2015 - 2025, The Regents of the University of Minnesota
* Distributed under the terms of the GNU General Public License.
*/

namespace Pscf { 
namespace Prdc { 

   /**
   * Return true if lattice type has variable angle parameters.
   *
   * This base template is not defined, but explicit specializations
   * are defined for D=1, 2, and 3.
   *
   * \param lattice  lattice system enumeration value
   * \return true iff lattice type has variable angles
   */
   template <int D>
   bool hasVariableAngle(typename UnitCell<D>::LatticeSystem lattice);

   // Explicit specializations

   template <>
   bool hasVariableAngle<1>(typename UnitCell<1>::LatticeSystem lattice);

   template <>
   bool hasVariableAngle<2>(typename UnitCell<2>::LatticeSystem lattice);

   template <>
   bool hasVariableAngle<3>(typename UnitCell<3>::LatticeSystem lattice);

}
}
#endif

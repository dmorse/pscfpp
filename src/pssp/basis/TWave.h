#ifndef PSSP_TWAVE_H
#define PSSP_TWAVE_H
/*
* PSCF - Polymer Self-Consistent Field Theory
*
* Copyright 2016, The Regents of the University of Minnesota
* Distributed under the terms of the GNU General Public License.
*/

#include <pscf/math/IntVec.h>

namespace Pscf {
namespace Pssp {

   using namespace Util;

   /**
   * Simple wave struct for use within Basis construction.
   *
   * \ingroup Pssp_Basis_Module
   */
   template <int D>
   struct TWave
   {
      double sqNorm;
      double phase;
      IntVec<D> indicesDft;
      IntVec<D> indicesBz;
   };

   /**
   * Comparator for TWave objects, based on TWave::sqNorm.
   *
   * Used to sort in ascending order of wavevector norm.
   *
   * \ingroup Pssp_Basis_Module
   */
   template <int D>
   struct TWaveNormComp {

      /**
      * Function (a, b) returns true iff a.sqNorm < b.sqNorm.
      */ 
      bool operator() (const TWave<D>& a, const TWave<D>& b) const
      {  return (a.sqNorm < b.sqNorm); }

   };

   /**
   * Comparator for TWave objects, based on TWave::indicesDft.
   *
   * Used to sort set of unique waves in ascending order of dft indices.
   *
   * \ingroup Pssp_Basis_Module
   */
   template <int D>
   struct TWaveDftComp {

      /**
      * Function (a, b) returns true iff a.indicesDft < b.indicesDft
      */ 
      bool operator() (const TWave<D>& a, const TWave<D>& b) const
      {  return (a.indicesDft < b.indicesDft); }

   };

   /**
   * Comparator for TWave objects, based on TWave::indicesBz.
   *
   * Used to sort in descending order of Bz (Brillouin zone) indices.
   *
   * \ingroup Pssp_Basis_Module
   */
   template <int D>
   struct TWaveBzComp {

      /**
      * Function (a, b) returns true iff a.indicesBz > b.indicesBz
      */ 
      bool operator() (const TWave<D>& a, const TWave<D>& b) const
      {  return (a.indicesBz > b.indicesBz); }

   };

} // namespace Pscf:Pssp
} // namespace Pscf
#endif

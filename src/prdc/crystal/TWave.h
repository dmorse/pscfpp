#ifndef PRDC_TWAVE_H
#define PRDC_TWAVE_H

/*
* PSCF - Polymer Self-Consistent Field
*
* Copyright 2015 - 2025, The Regents of the University of Minnesota
* Distributed under the terms of the GNU General Public License.
*/

#include <pscf/math/IntVec.h>

namespace Pscf {
namespace Prdc {

   using namespace Util;

   /**
   * Simple wave struct for use within Basis construction.
   *
   * \ingroup Prdc_Crystal_Module
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
   * \ingroup Prdc_Crystal_Module
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
   * \ingroup Prdc_Crystal_Module
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
   * \ingroup Prdc_Crystal_Module
   */
   template <int D>
   struct TWaveBzComp {

      /**
      * Function (a, b) returns true iff a.indicesBz > b.indicesBz
      */
      bool operator() (const TWave<D>& a, const TWave<D>& b) const
      {  return (a.indicesBz > b.indicesBz); }

   };

} // namespace Pscf::Prdc
} // namespace Pscf
#endif

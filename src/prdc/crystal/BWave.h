#ifndef PRDC_BWAVE_H
#define PRDC_BWAVE_H

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
   * Wave struct designed for use within Basis construction.
   *
   * \ingroup Prdc_Crystal_Module
   */
   template <int D>
   struct BWave
   {
      double sqNorm;
      double phase;
      IntVec<D> indicesStd;
      IntVec<D> indicesMin;
   };

   /**
   * Comparator for BWave objects, based on BWave::sqNorm.
   *
   * Used to sort in ascending order of wavevector norm.
   *
   * \ingroup Prdc_Crystal_Module
   */
   template <int D>
   struct BWaveNormComp {

      /**
      * Function (a, b) returns true iff a.sqNorm < b.sqNorm.
      */
      bool operator() (const BWave<D>& a, const BWave<D>& b) const
      {  return (a.sqNorm < b.sqNorm); }

   };

   /**
   * Comparator for BWave objects, based on BWave::indicesStd.
   *
   * Used to sort set of unique waves in ascending order of std indices.
   *
   * \ingroup Prdc_Crystal_Module
   */
   template <int D>
   struct BWaveStdComp {

      /**
      * Function (a, b) returns true iff a.indicesStd < b.indicesStd
      */
      bool operator() (const BWave<D>& a, const BWave<D>& b) const
      {  return (a.indicesStd < b.indicesStd); }

   };

   /**
   * Comparator for BWave objects, based on BWave::indicesMin.
   *
   * Used to sort in descending order of minimal indices.
   *
   * \ingroup Prdc_Crystal_Module
   */
   template <int D>
   struct BWaveMinComp {

      /**
      * Function (a, b) returns true iff a.indicesMin > b.indicesMin
      */
      bool operator() (const BWave<D>& a, const BWave<D>& b) const
      {  return (a.indicesMin > b.indicesMin); }

   };

} // namespace Pscf::Prdc
} // namespace Pscf
#endif

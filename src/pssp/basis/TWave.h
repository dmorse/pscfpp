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

   template <int D>
   struct TWaveNormComp {

      bool operator() (const TWave<D>& a, const TWave<D>& b) const
      {  return (a.sqNorm < b.sqNorm); }

   };

   template <int D>
   struct TWaveDftComp {

      bool operator() (const TWave<D>& a, const TWave<D>& b) const
      {  return (a.indicesDft < b.indicesDft); }

   };

   template <int D>
   struct TWaveBzComp {

      bool operator() (const TWave<D>& a, const TWave<D>& b) const
      {  return (a.indicesBz > b.indicesBz); }

   };

} // namespace Pscf:Pssp
} // namespace Pscf
#endif

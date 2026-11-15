/*
* PSCF - Polymer Self-Consistent Field
*
* Copyright 2015 - 2026, David Morse (morse012@um.edu) and co-authors
* Distributed under the terms of the GNU General Public License.
*/

#include "shiftToMinimum.h"
#include <pscf/math/RealVec.h>
#include <util/containers/FArray.h>
#include <util/global.h>
#include <cmath>

namespace Pscf {

   using namespace Util;

   namespace {

      /*
      * Scalar product of the reciprocal lattice vectors with integer
      * coordinates a and b.
      */
      template <int D>
      double kDot(IntVec<D> const & a, IntVec<D> const & b,
                  UnitCell<D> const & cell)
      {
         RealVec<D> ga(0.0);
         RealVec<D> gb(0.0);
         RealVec<D> p;
         for (int i = 0; i < D; ++i) {
            p.multiply(cell.kBasis(i), a[i]);
            ga += p;
            p.multiply(cell.kBasis(i), b[i]);
            gb += p;
         }
         double value = 0.0;
         for (int i = 0; i < D; ++i) {
            value += ga[i]*gb[i];
         }
         return value;
      }

      /*
      * Image shift vectors (see imageShiftNeighbors) of the last mesh
      * and unit cell used by this thread, so that the Selling reduction
      * is done once per mesh and cell, not once per wavevector.
      */
      template <int D>
      FSArray<IntVec<D>, 14> const &
      cachedNeighbors(IntVec<D> const & d, UnitCell<D> const & cell)
      {
         thread_local bool isValid = false;
         thread_local IntVec<D> dCache;
         thread_local double kBasisCache[D][D];
         thread_local FSArray<IntVec<D>, 14> neighbors;

         bool same = isValid && (dCache == d);
         for (int i = 0; same && i < D; ++i) {
            for (int j = 0; j < D; ++j) {
               if (kBasisCache[i][j] != cell.kBasis(i)[j]) same = false;
            }
         }
         if (!same) {
            isValid = false;
            imageShiftNeighbors<D>(d, cell, neighbors);
            dCache = d;
            for (int i = 0; i < D; ++i) {
               for (int j = 0; j < D; ++j) {
                  kBasisCache[i][j] = cell.kBasis(i)[j];
               }
            }
            isValid = true;
         }
         return neighbors;
      }

      /*
      * Minimum image of v (see shiftToMinimum).
      */
      template <int D>
      IntVec<D> minimumImage(IntVec<D> const & v, IntVec<D> const & d,
                             UnitCell<D> const & cell)
      {
         const double epsilon = 1.0E-8;
         FSArray<IntVec<D>, 14> const & nb = cachedNeighbors<D>(d, cell);
         const int n = nb.size();

         // Start from the image with components in (-d[i]/2, d[i]/2],
         // which is the same for all images v of a wavevector
         IntVec<D> r;
         for (int i = 0; i < D; ++i) {
            r[i] = v[i] % d[i];
            if (r[i] < 0) r[i] += d[i];
            if (2*r[i] > d[i]) r[i] -= d[i];
         }

         // Move r by the vector nb[j] that lowers |k|^2 most, for as
         // long as one lowers it by more than epsilon. No vector nb[j]
         // then lowers |k|^2, so r is a minimum image.
         double Gsq = cell.ksq(r);
         IntVec<D> u;
         for (;;) {
            int jBest = -1;
            double GsqBest = Gsq - epsilon;
            for (int j = 0; j < n; ++j) {
               u.add(r, nb[j]);
               double g = cell.ksq(u);
               if (g < GsqBest) {
                  jBest = j;
                  GsqBest = g;
               }
            }
            if (jBest < 0) break;
            r += nb[jBest];
            Gsq = GsqBest;
         }

         // Move r to a larger (in the order of IntVec) image r + nb[j]
         // of equal |k|^2, for as long as there is one. This ends at
         // the largest of the tied minimum images, which are connected
         // by the vectors nb[j].
         bool changed = true;
         while (changed) {
            changed = false;
            for (int j = 0; j < n; ++j) {
               u.add(r, nb[j]);
               if ((cell.ksq(u) < Gsq + epsilon) && (r < u)) {
                  r = u;
                  changed = true;
               }
            }
         }
         return r;
      }

   }

   template <> 
   IntVec<1> shiftToMinimum(IntVec<1> const & v, 
                            IntVec<1> const & d, 
                            UnitCell<1> const & cell)
   {
      IntVec<1> u;
      if (v[0] > d[0]/2) {
         u[0] = v[0] - d[0];
      } else 
      if (v[0] < -(d[0]/2) ) {
         u[0] = v[0] + d[0];
      } else {
         u[0] = v[0];
      }
      UTIL_ASSERT(abs(u[0]) <= d[0]/2);
      return u;
   }

   template <>
   IntVec<2> shiftToMinimum(IntVec<2> const & v, 
                            IntVec<2> const & d, 
                            UnitCell<2> const & cell)
   {  return minimumImage<2>(v, d, cell); }

   template <>
   IntVec<3> shiftToMinimum(IntVec<3> const & v, 
                            IntVec<3> const & d, 
                            UnitCell<3> const & cell)
   {  return minimumImage<3>(v, d, cell); }

   /*
   * Voronoi-relevant vectors of the lattice of image shifts.
   */
   template <int D>
   void imageShiftNeighbors(IntVec<D> const & d,
                            UnitCell<D> const & cell,
                            FSArray<IntVec<D>, 14>& neighbors)
   {
      neighbors.clear();
      if (D == 1) {
         neighbors.append(d);
         IntVec<D> w;
         w.negate(d);
         neighbors.append(w);
         return;
      }

      // Superbase: v[i+1] = d[i]*e_i, v[0] = -(v[1] + ... + v[D])
      FArray<IntVec<D>, D+1> v;
      v[0] = IntVec<D>(0);
      for (int i = 0; i < D; ++i) {
         v[i+1] = IntVec<D>(0);
         v[i+1][i] = d[i];
         v[0] -= v[i+1];
      }

      // Selling reduction: while some v_i.v_j > 0, replace v_i by -v_i
      // and add v_i to each other v_k (3D) or 2*v_i to the third (2D).
      // This keeps a superbase of the same lattice and lowers the sum
      // of the |v_k|^2 by 2 v_i.v_j (3D) or 4 v_i.v_j (2D). A product
      // v_i.v_j counts as positive if it exceeds rounding errors
      // relative to |v_i||v_j|, so that nearly orthogonal pairs (e.g.,
      // an angle of pi/2 rounded to 10 decimals) are also reduced.
      int iter = 0;
      for (;;) {
         double pMax = 0.0;
         int iMax = -1;
         int jMax = -1;
         for (int i = 0; i <= D; ++i) {
            for (int j = i + 1; j <= D; ++j) {
               double p = kDot<D>(v[i], v[j], cell);
               double s = std::sqrt(kDot<D>(v[i], v[i], cell)
                                    * kDot<D>(v[j], v[j], cell));
               if (p > 1.0E-14*s && p > pMax) {
                  pMax = p;
                  iMax = i;
                  jMax = j;
               }
            }
         }
         if (iMax < 0) break;
         for (int k = 0; k <= D; ++k) {
            if (k == iMax || k == jMax) continue;
            v[k] += v[iMax];
            if (D == 2) v[k] += v[iMax];
         }
         v[iMax].negate();
         ++iter;
         if (iter > 10000) {
            UTIL_THROW("Selling reduction of image shifts did not converge");
         }
      }

      // Neighbors: +-v_i, and in 3D +-(v_0 + v_i)
      IntVec<D> w;
      for (int i = 0; i <= D; ++i) {
         neighbors.append(v[i]);
         w.negate(v[i]);
         neighbors.append(w);
      }
      if (D == 3) {
         for (int i = 1; i <= D; ++i) {
            w.add(v[0], v[i]);
            neighbors.append(w);
            w.negate();
            neighbors.append(w);
         }
      }
   }

   template void imageShiftNeighbors<1>(IntVec<1> const &,
                                        UnitCell<1> const &,
                                        FSArray<IntVec<1>, 14>&);
   template void imageShiftNeighbors<2>(IntVec<2> const &,
                                        UnitCell<2> const &,
                                        FSArray<IntVec<2>, 14>&);
   template void imageShiftNeighbors<3>(IntVec<3> const &,
                                        UnitCell<3> const &,
                                        FSArray<IntVec<3>, 14>&);

}

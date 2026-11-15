#ifndef PSCF_MISC_TEST_H
#define PSCF_MISC_TEST_H

#include <test/UnitTest.h>
#include <test/UnitTestRunner.h>

#include <pscf/field/cpp/FFT.h>
#include <pscf/crystal/sortWaves.h>
#include <pscf/crystal/UnitCell.h>
#include <pscf/crystal/shiftToMinimum.h>
#include <pscf/mesh/Mesh.h>
#include <pscf/mesh/MeshIterator.h>
#include <pscf/math/Sort.h>
#include <util/math/Constants.h>
#include <util/containers/FSArray.h>
#include <util/format/Int.h>
#include <util/format/Dbl.h>

#include <vector>
#include <iostream>
#include <fstream>

using namespace Util;
using namespace Pscf;

class MiscTest : public UnitTest 
{

public:

   void setUp()
   {  setVerbose(0); }

   void tearDown()
   {  setVerbose(0); }
 
   void testSortWavesReal() 
   {
      printMethod(TEST_FUNC);

      UnitCell<3> cell;
      FSArray<double, 6> parameters;
      double a = 2.0;
      parameters.append(a);
      cell.set(UnitCell<3>::Cubic, parameters);

      IntVec<3> dimensions;
      dimensions[0] = 8;
      dimensions[1] = 8;
      dimensions[2] = 8;
      IntVec<3> kDimensions;
      int kSize;
      FFT<3,CPT>::computeKMesh(dimensions, kDimensions, kSize);
      Mesh<3> kMesh(kDimensions);

      std::vector< Sort::Item<double> > items;
      GArray< Sort::Bunch > bunches;
      double epsilon = 1.0E-8;
      bool isRealField = true;
      sortWaves(cell, dimensions, items, bunches, epsilon, isRealField);

      const double twoPi = 2.0 * Constants::Pi;
      const double bSqInv = (a*a)/(twoPi * twoPi);

      double oldVal, value;
      int nb = bunches.size();
      int begin, end, size, ib, iw, rank, itot;
      IntVec<3> position;
      itot = 0; 
      for (ib = 0; ib < nb; ++ib) {
         begin = bunches[ib][0];
         end = bunches[ib][1];
         if (ib > 0) {
           TEST_ASSERT(begin == bunches[ib-1][1]);
         }
         oldVal = items[begin].value;
         for (iw = begin; iw < end; ++iw) {
             value = items[iw].value;
             TEST_ASSERT( std::abs( value - oldVal) < epsilon);
             if (iw > 0) {
                TEST_ASSERT(value >= items[iw-1].value);
             } 
             ++itot;
         }
         if (verbose() > 0) {
            size = end - begin;
            std::cout << std::endl << Int(begin) << Int(end) << Int(size);
            for (iw = begin; iw < end; ++iw) {
               value = items[iw].value;
               value *= bSqInv;
               rank = items[iw].id;
               position = kMesh.position(rank);
               position = shiftToMinimum(position, dimensions, cell);
               std::cout << std::endl 
                         << Int(iw)
                         << Int(rank)
                         << Dbl(value);
               for (int k=0; k < 3; ++k) {
                  std::cout << Int(position[k]);
               }
            }
         }
      }
      TEST_ASSERT(itot == (int) items.size());
   }

   void testSortWavesComplex() 
   {
      printMethod(TEST_FUNC);

      UnitCell<3> cell;
      FSArray<double, 6> parameters;
      double a = 2.0;
      parameters.append(a);
      cell.set(UnitCell<3>::Cubic, parameters);

      IntVec<3> dimensions;
      dimensions[0] = 8;
      dimensions[1] = 8;
      dimensions[2] = 8;
      IntVec<3> kDimensions = dimensions;
      Mesh<3> kMesh(kDimensions);

      std::vector< Sort::Item<double> > items;
      GArray< Sort::Bunch > bunches;
      double epsilon = 1.0E-8;
      bool isRealField = false;
      sortWaves(cell, dimensions, items, bunches, epsilon, isRealField);

      const double twoPi = 2.0 * Constants::Pi;
      const double bSqInv = (a*a)/(twoPi * twoPi);

      double oldVal, value;
      int nb = bunches.size();
      int begin, end, size, ib, iw, rank, itot;
      IntVec<3> position;
      itot = 0; 
      for (ib = 0; ib < nb; ++ib) {
         begin = bunches[ib][0];
         end = bunches[ib][1];
         if (ib > 0) {
           TEST_ASSERT(begin == bunches[ib-1][1]);
         }
         oldVal = items[begin].value;
         for (iw = begin; iw < end; ++iw) {
             value = items[iw].value;
             TEST_ASSERT( std::abs( value - oldVal) < epsilon);
             if (iw > 0) {
                TEST_ASSERT(value >= items[iw-1].value);
             } 
             ++itot;
         }
         if (verbose() > 0) {
            size = end - begin;
            std::cout << std::endl << Int(begin) << Int(end) << Int(size);
            for (iw = begin; iw < end; ++iw) {
               value = items[iw].value;
               value *= bSqInv;
               rank = items[iw].id;
               position = kMesh.position(rank);
               position = shiftToMinimum(position, dimensions, cell);
               std::cout << std::endl 
                         << Int(iw)
                         << Int(rank)
                         << Dbl(value);
               for (int k=0; k < 3; ++k) {
                  std::cout << Int(position[k]);
               }
            }
         }
      }
      TEST_ASSERT(itot == (int) items.size());
   }

   /*
   * Skewed cells on meshes with very different grid spacings, for
   * which the minimum images of some wavevectors lie more than two
   * image shifts away from the wavevector.
   */
   void testShiftToMinimumSkewed()
   {
      printMethod(TEST_FUNC);

      FSArray<double, 6> parameters;
      parameters.append(2.0);
      parameters.append(2.0);
      parameters.append(2.4);
      UnitCell<2> cell2;
      cell2.set(UnitCell<2>::Oblique, parameters);
      IntVec<2> d2;
      d2[0] = 8;
      d2[1] = 64;
      checkMinimumImages<2>(cell2, d2, 8);

      parameters[2] = 2.0;
      parameters.append(2.8);
      UnitCell<3> cell3;
      cell3.set(UnitCell<3>::Monoclinic, parameters);
      IntVec<3> d3;
      d3[0] = 4;
      d3[1] = 2;
      d3[2] = 32;
      checkMinimumImages<3>(cell3, d3, 6);
   }

   /*
   * The image returned by shiftToMinimum must not depend on which
   * image of a wavevector is passed in. In this skewed cell, on a mesh
   * with very different grid spacings, tied minimum images lie several
   * image shifts apart.
   */
   void testShiftToMinimumRepresentative()
   {
      printMethod(TEST_FUNC);

      FSArray<double, 6> parameters;
      parameters.append(1.0);
      parameters.append(1.0);
      parameters.append(1.0);
      parameters.append(0.0);
      parameters.append(1.2);
      parameters.append(0.5*Constants::Pi);
      UnitCell<3> cell;
      cell.set(UnitCell<3>::Triclinic, parameters);
      IntVec<3> d;
      d[0] = 4;
      d[1] = 4;
      d[2] = 32;

      IntVec<3> v, r, u;
      MeshIterator<3> vItr(d);
      for (vItr.begin(); !vItr.atEnd(); ++vItr) {
         v = vItr.position();
         r = shiftToMinimum(v, d, cell);
         TEST_ASSERT(shiftToMinimum(r, d, cell) == r);
         for (int i = 0; i < 3; ++i) {
            u = v;
            u[i] += d[i];
            TEST_ASSERT(shiftToMinimum(u, d, cell) == r);
            u[i] -= 2*d[i];
            TEST_ASSERT(shiftToMinimum(u, d, cell) == r);
         }
      }
   }

private:

   /*
   * Check, for all mesh points v, that shiftToMinimum returns an image
   * of v whose |k|^2 is not larger than that of any image v + s, with
   * s[i] = m[i]*d[i] and all |m[i]| <= q.
   */
   template <int D>
   void checkMinimumImages(UnitCell<D> const & cell, IntVec<D> const & d,
                           int q)
   {
      IntVec<D> mDims;
      for (int i = 0; i < D; ++i) mDims[i] = 2*q + 1;
      MeshIterator<D> vItr(d);
      MeshIterator<D> mItr(mDims);
      IntVec<D> v, r, u;
      for (vItr.begin(); !vItr.atEnd(); ++vItr) {
         v = vItr.position();
         r = shiftToMinimum(v, d, cell);
         for (int i = 0; i < D; ++i) {
            TEST_ASSERT((r[i] - v[i]) % d[i] == 0);
         }
         double kSq = cell.ksq(r);
         for (mItr.begin(); !mItr.atEnd(); ++mItr) {
            for (int i = 0; i < D; ++i) {
               u[i] = v[i] + (mItr.position(i) - q)*d[i];
            }
            TEST_ASSERT(kSq <= cell.ksq(u)*(1.0 + 1.0E-10));
         }
      }
   }

};

TEST_BEGIN(MiscTest)
TEST_ADD(MiscTest, testSortWavesReal)
TEST_ADD(MiscTest, testSortWavesComplex)
TEST_ADD(MiscTest, testShiftToMinimumSkewed)
TEST_ADD(MiscTest, testShiftToMinimumRepresentative)
TEST_END(MiscTest)

#endif

#ifndef PRDC_MISC_TEST_H
#define PRDC_MISC_TEST_H

#include <test/UnitTest.h>
#include <test/UnitTestRunner.h>

#include <prdc/cpu/FFT.h>
#include <prdc/crystal/sortWaves.h>
#include <prdc/crystal/UnitCell.h>
#include <pscf/mesh/Mesh.h>
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
using namespace Pscf::Prdc;

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
      Cpu::FFT<3>::computeKMesh(dimensions, kDimensions, kSize);
      Mesh<3> kMesh(kDimensions);

      std::vector< Sort::Item<double> > items;
      std::vector< Sort::Bunch > bunches;
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
      TEST_ASSERT(itot == items.size());
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
      std::vector< Sort::Bunch > bunches;
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
      TEST_ASSERT(itot == items.size());
   }

};

TEST_BEGIN(MiscTest)
TEST_ADD(MiscTest, testSortWavesReal)
TEST_ADD(MiscTest, testSortWavesComplex)
TEST_END(MiscTest)

#endif

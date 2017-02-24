#ifndef PSCF_UNIT_CELL_TEST_H
#define PSCF_UNIT_CELL_TEST_H

#include <test/UnitTest.h>
#include <test/UnitTestRunner.h>

#include <pscf/crystal/UnitCell.h>
#include <util/math/Constants.h>

#include <iostream>
#include <fstream>

using namespace Util;
using namespace Pscf;

class UnitCellTest : public UnitTest 
{

public:

   void setUp()
   {}

   void tearDown()
   {}
 
   template <int D>
   bool isValidReciprocal(UnitCell<D> cell)
   {
      double sum;
      double twoPi = 2.0*Constants::Pi;
      int i, j, k;
      for (i=0; i < D; ++i ) {
         for (j=0; j < D; ++j ) {
            sum = 0.0;
            for (k=0; k < D; ++k ) {
               sum += cell.rBasis(i)[k]*cell.kBasis(j)[k];  
            }
            if (i == j) {
               sum -= twoPi;
            }
            if (std::abs(sum) > 1.0E-8) {
               return false;
            }
         }
      }
      return true;
   }

   void test1DLamellar() {
      printMethod(TEST_FUNC);
      printEndl();

      UnitCell<1> v;
      std::ifstream in;
      openInputFile("in/Lamellar", in);

      in >> v;
      std::cout.width(20);
      std::cout.precision(6);
      std::cout << v << std::endl ;

      TEST_ASSERT(isValidReciprocal(v));
   }

   void test2DSquare() {
      printMethod(TEST_FUNC);
      printEndl();

      UnitCell<2> v;
      std::ifstream in;
      openInputFile("in/Square", in);

      in >> v;
      std::cout.width(20);
      std::cout.precision(6);
      std::cout << v << std::endl ;

      TEST_ASSERT(isValidReciprocal(v));
   }

};

TEST_BEGIN(UnitCellTest)
TEST_ADD(UnitCellTest, test1DLamellar)
TEST_ADD(UnitCellTest, test2DSquare)
TEST_END(UnitCellTest)

#endif

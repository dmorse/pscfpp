#ifndef PRDC_CPU_FIELD_COMPARISON_TEST_H
#define PRDC_CPU_FIELD_COMPARISON_TEST_H

#include <test/UnitTest.h>
#include <test/UnitTestRunner.h>

#include <prdc/field/cpp/RField.h>
#include <prdc/field/cpp/RFieldComparison.h>

#include <prdc/field/cpp/RFieldDft.h>
#include <prdc/field/cpp/RFieldDftComparison.h>

#include <prdc/field/cpp/CField.h>
#include <prdc/field/cpp/CFieldComparison.h>

#include <util/format/Dbl.h>

using namespace Util;
using namespace Pscf::Prdc;

class CpuFieldComparisonTest : public UnitTest 
{

public:

   void setUp()
   {  setVerbose(0);  }

   void tearDown() 
   {}

   void testRFieldComparison_1D()
   {
      printMethod(TEST_FUNC);

      RField<1,CPT> rf_0, rf_1;
      int n = 10;
      IntVec<1> dimensions;
      dimensions[0] = n;
      rf_0.allocate(dimensions);
      rf_1.allocate(dimensions);
      int size = rf_0.capacity();
      TEST_ASSERT(size == n);
      for (int i = 0; i < n; ++i) {
         rf_0[i] = 2.0;
         rf_1[i] = 2.001;
      }
      RFieldComparison<1,CPT> comparison;
      comparison.compare(rf_0,  rf_1);
      if (verbose() > 0) {
         std::cout << "\n";
         std::cout << "MaxDiff = " 
                   << Dbl(comparison.maxDiff(), 20, 12) << "\n";
         std::cout << "RmsDiff = " 
                   << Dbl(comparison.rmsDiff(), 20, 12) << "\n";
      }
      TEST_ASSERT(comparison.maxDiff() < 0.0011);
      TEST_ASSERT(comparison.maxDiff() > 0.0009);
      TEST_ASSERT(comparison.rmsDiff() < 0.0011);
      TEST_ASSERT(comparison.rmsDiff() > 0.0009);

   }

   void testRFieldComparison_2D()
   {
      printMethod(TEST_FUNC);

      RField<2,CPT> rf_0, rf_1;
      int m = 5;
      int n = 10;
      IntVec<2> dimensions;
      dimensions[0] = m;
      dimensions[1] = n;
      rf_0.allocate(dimensions);
      rf_1.allocate(dimensions);
      int size = rf_0.capacity();
      TEST_ASSERT(size == m*n);
      for (int i = 0; i < size; ++i) {
         rf_0[i] = 2.0;
         rf_1[i] = 2.001;
      }
      RFieldComparison<2,CPT> comparison;
      comparison.compare(rf_0,  rf_1);
      if (verbose() > 0) {
         std::cout << "\n";
         std::cout << "MaxDiff = " 
                   << Dbl(comparison.maxDiff(), 20, 12) << "\n";
         std::cout << "RmsDiff = " 
                   << Dbl(comparison.rmsDiff(), 20, 12) << "\n";
      }
      TEST_ASSERT(comparison.maxDiff() < 0.0011);
      TEST_ASSERT(comparison.maxDiff() > 0.0009);
      TEST_ASSERT(comparison.rmsDiff() < 0.0011);
      TEST_ASSERT(comparison.rmsDiff() > 0.0009);

   }

   void testRFieldDftComparison_2D()
   {
      printMethod(TEST_FUNC);
      //setVerbose(1);

      RFieldDft<2,CPT> f_0, f_1;
      int m = 5;
      int n = 10;
      IntVec<2> dimensions;
      dimensions[0] = m;
      dimensions[1] = n;
      f_0.allocate(dimensions);
      f_1.allocate(dimensions);
      int size = f_0.capacity();
      TEST_ASSERT(size == f_1.capacity());
      for (int i = 0; i < size; ++i) {
         f_0[i][0] = 2.0;
         f_0[i][1] = 2.0;
         f_1[i][0] = 1.999;
         f_1[i][1] = 2.001;
      }
      RFieldDftComparison<2,CPT> comparison;
      comparison.compare(f_0,  f_1);
      if (verbose() > 0) {
         std::cout << "\n";
         std::cout << "MaxDiff = " 
                   << Dbl(comparison.maxDiff(), 20, 12) << "\n";
         std::cout << "RmsDiff = " 
                   << Dbl(comparison.rmsDiff(), 20, 12) << "\n";
      }
      TEST_ASSERT(comparison.maxDiff() < 0.00142);
      TEST_ASSERT(comparison.maxDiff() > 0.00141);
      TEST_ASSERT(comparison.rmsDiff() < 0.00142);
      TEST_ASSERT(comparison.rmsDiff() > 0.00141);

   }

   void testCFieldComparison_2D()
   {
      printMethod(TEST_FUNC);
      //setVerbose(1);

      CField<2,CPT> f_0, f_1;
      int m = 5;
      int n = 10;
      IntVec<2> dimensions;
      dimensions[0] = m;
      dimensions[1] = n;
      f_0.allocate(dimensions);
      f_1.allocate(dimensions);
      int size = f_0.capacity();
      TEST_ASSERT(size == m*n);
      for (int i = 0; i < size; ++i) {
         f_0[i][0] = 2.0;
         f_0[i][1] = 2.0;
         f_1[i][0] = 1.999;
         f_1[i][1] = 2.001;
      }
      CFieldComparison<2,CPT> comparison;
      comparison.compare(f_0,  f_1);
      if (verbose() > 0) {
         std::cout << "\n";
         std::cout << "MaxDiff = " 
                   << Dbl(comparison.maxDiff(), 20, 12) << "\n";
         std::cout << "RmsDiff = " 
                   << Dbl(comparison.rmsDiff(), 20, 12) << "\n";
      }
      TEST_ASSERT(comparison.maxDiff() < 0.00142);
      TEST_ASSERT(comparison.maxDiff() > 0.00141);
      TEST_ASSERT(comparison.rmsDiff() < 0.00142);
      TEST_ASSERT(comparison.rmsDiff() > 0.00141);

   }

};


TEST_BEGIN(CpuFieldComparisonTest)
TEST_ADD(CpuFieldComparisonTest, testRFieldComparison_1D)
TEST_ADD(CpuFieldComparisonTest, testRFieldComparison_2D)
TEST_ADD(CpuFieldComparisonTest, testRFieldDftComparison_2D)
TEST_ADD(CpuFieldComparisonTest, testCFieldComparison_2D)
TEST_END(CpuFieldComparisonTest)

#endif

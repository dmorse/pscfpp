#ifndef PSPC_FIELD_COMPARISON_TEST_H
#define PSPC_FIELD_COMPARISON_TEST_H

#include <test/UnitTest.h>
#include <test/UnitTestRunner.h>

#include <pspc/field/RField.h>
#include <pspc/field/RFieldDft.h>
#include <pspc/field/RFieldComparison.h>
#include <pspc/field/KFieldComparison.h>

#include <util/format/Dbl.h>

using namespace Util;
using namespace Pscf::Pspc;

class FieldComparisonTest : public UnitTest 
{
public:

   void setUp()
   {  setVerbose(0);  }

   void tearDown() 
   {}

   void testRFieldComparison()
   {
      printMethod(TEST_FUNC);

      RField<1> rf_0, rf_1;
      int n = 10;
      rf_0.allocate(n);
      rf_1.allocate(n);
      for (int i = 0; i < n; ++i) {
         rf_0[i] = 2.0;
         rf_1[i] = 2.001;
      }
      RFieldComparison<1> comparison;
      comparison.compare(rf_0,  rf_1);
      setVerbose(1);
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

};


TEST_BEGIN(FieldComparisonTest)
TEST_ADD(FieldComparisonTest, testRFieldComparison)
TEST_END(FieldComparisonTest)

#endif

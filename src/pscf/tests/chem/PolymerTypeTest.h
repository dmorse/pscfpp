#ifndef POLYMER_TYPE_TEST_H
#define POLYMER_TYPE_TEST_H

#include <test/UnitTest.h>
#include <test/UnitTestRunner.h>

#include <pscf/chem/PolymerType.h>

#include <fstream>

using namespace Pscf;
//using namespace Util;

class PolymerTypeTest : public UnitTest 
{

public:

   void setUp()
   {}

   void tearDown()
   {}

   void testReadWrite() {
      printMethod(TEST_FUNC);
      printEndl();

      std::ifstream in;
      openInputFile("in/PolymerType", in);

      PolymerType::Enum b;
      in >> b;
      TEST_ASSERT(b == PolymerType::Branched);
      std::cout << b << "  ";
      in >> b;
      TEST_ASSERT(b == PolymerType::Linear);
      std::cout << b << std::endl ;

      // If uncommented out, this one fails to read "Thingy"
      //in >> b;
   }

};

TEST_BEGIN(PolymerTypeTest)
TEST_ADD(PolymerTypeTest, testReadWrite)
TEST_END(PolymerTypeTest)

#endif

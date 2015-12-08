#ifndef MONOMER_TEST_H
#define MONOMER_TEST_H

#include <test/UnitTest.h>
#include <test/UnitTestRunner.h>

#include <chem/Monomer.h>

#include <fstream>

using namespace Chem;
//using namespace Util;

class MonomerTest : public UnitTest 
{

public:

   void setUp()
   {}

   void tearDown()
   {}

  
   void testConstructor()
   {
      printMethod(TEST_FUNC);
      Monomer v;
   } 

   void testReadWrite() {
      printMethod(TEST_FUNC);
      printEndl();

      Monomer v;
      std::ifstream in;
      openInputFile("in/Monomer", in);

      in >> v;
      std::cout << v << std::endl ;
   }

};

TEST_BEGIN(MonomerTest)
TEST_ADD(MonomerTest, testConstructor)
TEST_ADD(MonomerTest, testReadWrite)
TEST_END(MonomerTest)

#endif

#ifndef MONOMER_TEST_H
#define MONOMER_TEST_H

#include <test/UnitTest.h>
#include <test/UnitTestRunner.h>

#include <pscf/chem/Monomer.h>

#include <fstream>
#include <sstream>

using namespace Pscf;
using namespace Util;

class MonomerTest : public UnitTest 
{

public:

   void setUp()
   {
      //setVerbose(1);
   }

   void tearDown()
   {}

  
   void testConstructor()
   {
      printMethod(TEST_FUNC);
      Monomer v;
   } 

   void testReadWrite() {
      printMethod(TEST_FUNC);

      Monomer v;
      std::ifstream in;
      openInputFile("in/Monomer", in);

      in >> v;
      TEST_ASSERT(eq(v.kuhn(), 5.0));
      if (verbose() > 0) {
         printEndl();
         std::cout << v << std::endl ;
      }
   }

   /*
   * Output to a stream set to fixed format: setf(scientific) without
   * the floatfield mask gave hexadecimal output, which operator >>
   * cannot read.
   */
   void testWriteFixedStream()
   {
      printMethod(TEST_FUNC);

      Monomer v, u;
      v.setKuhn(1.25);
      std::stringstream stream;
      stream << std::fixed << v;
      TEST_ASSERT(stream.str().find("0x") == std::string::npos);
      stream >> u;
      TEST_ASSERT(u.kuhn() == 1.25);
   }

};

TEST_BEGIN(MonomerTest)
TEST_ADD(MonomerTest, testConstructor)
TEST_ADD(MonomerTest, testReadWrite)
TEST_ADD(MonomerTest, testWriteFixedStream)
TEST_END(MonomerTest)

#endif

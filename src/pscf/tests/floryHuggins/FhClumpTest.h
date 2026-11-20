#ifndef PSCF_FH_CLUMP_TEST_H
#define PSCF_FH_CLUMP_TEST_H

#include <test/UnitTest.h>
#include <test/UnitTestRunner.h>

#include <pscf/floryHuggins/FhClump.h>
#include <util/misc/Log.h>

#include <fstream>
#include <sstream>

using namespace Pscf;
using namespace Util;

class FhClumpTest : public UnitTest 
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
      FhClump clump;
   } 

   void testSetters()
   {
      printMethod(TEST_FUNC);
      FhClump clump;
      clump.setMonomerId(0);
      clump.setSize(2.0);
      TEST_ASSERT(clump.monomerId() == 0);
      TEST_ASSERT(eq(clump.size(), 2.0));
   } 

   void testReadWrite() {
      printMethod(TEST_FUNC);

      FhClump clump;
      std::ifstream in;
      openInputFile("in/FhClump", in);

      in >> clump;
      TEST_ASSERT(clump.monomerId() == 0);
      TEST_ASSERT(eq(clump.size(), 2.0));
      if (verbose() > 0) {
         printEndl();
         Log::file() << clump << std::endl ;
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

      FhClump v, u;
      v.setMonomerId(1);
      v.setSize(2.5);
      std::stringstream stream;
      stream << std::fixed << v;
      TEST_ASSERT(stream.str().find("0x") == std::string::npos);
      stream >> u;
      TEST_ASSERT(u.monomerId() == 1);
      TEST_ASSERT(u.size() == 2.5);
   }

};

TEST_BEGIN(FhClumpTest)
TEST_ADD(FhClumpTest, testConstructor)
TEST_ADD(FhClumpTest, testSetters)
TEST_ADD(FhClumpTest, testReadWrite)
TEST_ADD(FhClumpTest, testWriteFixedStream)
TEST_END(FhClumpTest)

#endif

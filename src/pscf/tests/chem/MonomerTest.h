#ifndef MONOMER_TEST_H
#define MONOMER_TEST_H

#include <test/UnitTest.h>
#include <test/UnitTestRunner.h>

#include <pscf/chem/Monomer.h>
#include <util/archives/MemoryOArchive.h>
#include <util/archives/MemoryIArchive.h>

#include <fstream>

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

   /*
   * Monomer::serialize took the archive by value, so it read from a
   * copy of a MemoryIArchive, whose cursor then did not advance, and it
   * could not write to a MemoryOArchive (private copy constructor).
   */
   void testSerialize()
   {
      printMethod(TEST_FUNC);

      // Read from data written as primitive types
      int id = 3;
      double kuhn = 1.25;
      int i1 = 13;
      MemoryOArchive oArchive;
      oArchive.allocate(64);
      oArchive << id;
      oArchive << kuhn;
      oArchive << i1;

      Monomer u;
      int i2 = 0;
      MemoryIArchive iArchive;
      iArchive = oArchive;
      iArchive >> u;
      iArchive >> i2;
      TEST_ASSERT(u.id() == 3);
      TEST_ASSERT(u.kuhn() == 1.25);
      TEST_ASSERT(i2 == 13);
      iArchive.release();

      // Write
      Monomer v;
      v.setId(4);
      v.setKuhn(2.5);
      MemoryOArchive oArchive2;
      oArchive2.allocate(64);
      oArchive2 << v;
      TEST_ASSERT(oArchive2.cursor()
                  == oArchive2.begin() + sizeof(int) + sizeof(double));
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

};

TEST_BEGIN(MonomerTest)
TEST_ADD(MonomerTest, testConstructor)
TEST_ADD(MonomerTest, testSerialize)
TEST_ADD(MonomerTest, testReadWrite)
TEST_END(MonomerTest)

#endif

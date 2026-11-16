#ifndef RP_SWEEP_PARAMETER_TEST_H
#define RP_SWEEP_PARAMETER_TEST_H

#include <test/UnitTest.h>
#include <test/UnitTestRunner.h>

#include <pscf/backend/cpp/CPT.h>
#include <rp/scft/sweep/SweepParameter.h>
#include <rp/fts/ramp/RampParameter.h>
#include <util/archives/MemoryOArchive.h>
#include <util/archives/MemoryIArchive.h>

#include <cstring>

using namespace Util;
using namespace Pscf;

class SweepParameterTest : public UnitTest
{

public:

   void setUp()
   {}

   void tearDown()
   {}

   /*
   * Load a parameter of type P into a new object from the record of a
   * block length parameter (type 0) with ids 0 and 1, and save it again.
   */
   template <class P>
   void checkSerialize()
   {
      unsigned int type = 0;
      int nId = 2;
      int id0 = 0;
      int id1 = 1;
      double initial = 0.5;
      double change = 0.25;
      MemoryOArchive record;
      record.allocate(64);
      record << type << nId << id0 << id1 << initial << change;
      size_t size = record.cursor() - record.begin();

      P p;
      MemoryIArchive iArchive;
      iArchive = record;
      iArchive >> p;
      iArchive.release();
      TEST_ASSERT(p.type() == "block");
      TEST_ASSERT(p.id(0) == 0);
      TEST_ASSERT(p.id(1) == 1);
      TEST_ASSERT(p.initial() == 0.5);
      TEST_ASSERT(p.change() == 0.25);

      MemoryOArchive oArchive;
      oArchive.allocate(64);
      oArchive << p;
      TEST_ASSERT((size_t)(oArchive.cursor() - oArchive.begin()) == size);
      TEST_ASSERT(std::memcmp(oArchive.begin(), record.begin(), size) == 0);
   }

   void testSweepParameterSerialize()
   {
      printMethod(TEST_FUNC);
      checkSerialize< Rp::SweepParameter<1, CPT> >();
   }

   void testRampParameterSerialize()
   {
      printMethod(TEST_FUNC);
      checkSerialize< Rp::RampParameter<1, CPT> >();
   }

};

TEST_BEGIN(SweepParameterTest)
TEST_ADD(SweepParameterTest, testSweepParameterSerialize)
TEST_ADD(SweepParameterTest, testRampParameterSerialize)
TEST_END(SweepParameterTest)

#endif

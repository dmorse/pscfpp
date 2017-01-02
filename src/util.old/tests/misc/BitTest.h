#ifndef BIT_TEST_H
#define BIT_TEST_H

#include <util/misc/Bit.h>
#include <util/global.h>

#include <test/UnitTest.h>
#include <test/UnitTestRunner.h>

using namespace Util;

class BitTest : public UnitTest 
{

public:

   void setUp()
   {};

   void tearDown()
   {};

   void testConstructor1() 
   {
      printMethod(TEST_FUNC);

      Bit bit(3);
      TEST_ASSERT(bit.mask() == 8);
   }

   void testConstructor2() 
   {
      printMethod(TEST_FUNC);

      Bit bit;
      bit.setMask(3);
      TEST_ASSERT(bit.mask() == 8);
   }

   void testIsSet() 
   {
      printMethod(TEST_FUNC);

      Bit bit(3);
      TEST_ASSERT(!bit.isSet(3));
      TEST_ASSERT(!bit.isSet(5));
      TEST_ASSERT(bit.isSet(8));
      TEST_ASSERT(bit.isSet(9));
      TEST_ASSERT(bit.isSet(12));
      TEST_ASSERT(bit.isSet(15));
      TEST_ASSERT(!bit.isSet(16));
      TEST_ASSERT(!bit.isSet(23));
      TEST_ASSERT(bit.isSet(24));
      TEST_ASSERT(bit.isSet(27));
      TEST_ASSERT(bit.isSet(31));
      TEST_ASSERT(!bit.isSet(32));
      TEST_ASSERT(!bit.isSet(35));
      TEST_ASSERT(!bit.isSet(39));
      TEST_ASSERT(bit.isSet(40));
      TEST_ASSERT(bit.mask() == 8);
   }

   void testSet() 
   {
      printMethod(TEST_FUNC);
      Bit bit(3);
      unsigned int flags  = 39;
      TEST_ASSERT(!bit.isSet(flags));
      bit.set(flags);
      TEST_ASSERT(bit.isSet(flags));
      TEST_ASSERT(flags == 47);
      bit.set(flags);
      TEST_ASSERT(bit.isSet(flags));
      TEST_ASSERT(flags == 47);
      TEST_ASSERT(bit.mask() == 8);
   }

   void testClear() 
   {
      printMethod(TEST_FUNC);
      Bit bit(3);
      unsigned int flags  = 47;
      TEST_ASSERT(bit.isSet(flags));
      bit.clear(flags);
      TEST_ASSERT(!bit.isSet(flags));
      TEST_ASSERT(flags == 39);
      bit.set(flags);
      TEST_ASSERT(flags == 47);
      TEST_ASSERT(bit.isSet(flags));
      TEST_ASSERT(bit.mask() == 8);
   }

};

TEST_BEGIN(BitTest)
TEST_ADD(BitTest, testConstructor1)
TEST_ADD(BitTest, testConstructor2)
TEST_ADD(BitTest, testIsSet)
TEST_ADD(BitTest, testSet)
TEST_ADD(BitTest, testClear)
TEST_END(BitTest)

#endif

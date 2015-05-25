#ifndef SETABLE_TEST_H
#define SETABLE_TEST_H

#include <util/misc/Setable.h>
#include <util/global.h>

#ifdef UTIL_MPI
#ifndef TEST_MPI
#define TEST_MPI
#endif
#endif

#include <test/UnitTest.h>
#include <test/UnitTestRunner.h>

using namespace Util;

class SetableTest : public UnitTest 
{

public:

   void setUp()
   {};

   void tearDown()
   {};

   void testDefaultConstructor() 
   {
      printMethod(TEST_FUNC);

      Setable<int> setable;
      TEST_ASSERT(!setable.isSet());
   }

   void testConstructorFromValue() 
   {
      printMethod(TEST_FUNC);

      Setable<int> setable(3);
      TEST_ASSERT(setable.isSet());
      TEST_ASSERT(setable.value() == 3);
      setable.unset();
      TEST_ASSERT(!setable.isSet());
   }

   void testAssignment() 
   {
      printMethod(TEST_FUNC);

      Setable<int> setable1(3);
      TEST_ASSERT(setable1.isSet());
      TEST_ASSERT(setable1.value() == 3);

      Setable<int> setable2;
      TEST_ASSERT(!setable2.isSet());
        
      setable2 = setable1;
      TEST_ASSERT(setable2.isSet());
      TEST_ASSERT(setable2.value() == 3);
   }

   void testAssignmentFromValue() 
   {
      printMethod(TEST_FUNC);

      Setable<int> setable;
      TEST_ASSERT(!setable.isSet());

      setable = 3;
      TEST_ASSERT(setable.isSet());
      TEST_ASSERT(setable.value() == 3);

      setable.unset();
      TEST_ASSERT(!setable.isSet());
   }

   void testSet() 
   {
      printMethod(TEST_FUNC);

      Setable<int> setable;
      TEST_ASSERT(!setable.isSet());

      setable.set(3);
      TEST_ASSERT(setable.isSet());
      TEST_ASSERT(setable.value() == 3);

      setable.unset();
      TEST_ASSERT(!setable.isSet());
   }

   #ifdef UTIL_MPI
   void testIsValid() 
   {
      printMethod(TEST_FUNC);

      Setable<int> setable;
      TEST_ASSERT(!setable.isSet());

      setable.set(3);
      TEST_ASSERT(setable.isSet());
      TEST_ASSERT(setable.value() == 3);
      TEST_ASSERT(setable.isValid(MPI::COMM_WORLD));

      setable.unset();
      TEST_ASSERT(!setable.isSet());
   }
   #endif

};

TEST_BEGIN(SetableTest)
TEST_ADD(SetableTest, testDefaultConstructor)
TEST_ADD(SetableTest, testConstructorFromValue)
TEST_ADD(SetableTest, testAssignment)
TEST_ADD(SetableTest, testAssignmentFromValue)
TEST_ADD(SetableTest, testSet)
#ifdef UTIL_MPI
TEST_ADD(SetableTest, testIsValid)
#endif
TEST_END(SetableTest)

#endif

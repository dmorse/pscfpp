#include "UnitTest.h"
#include "UnitTestRunner.h"
#include "CompositeTestRunner.h"
#include <iostream>

/**
* This file demonstrates the usage of a composite test runner.
* 
* A CompositeTestRunner is a TestRunner that runs and accumulates
* statistics of the tests associated with several child TestRunner
* objects. The child runners are usually instances of UnitTestRunner.
*
* To demonstrate the usage, we define two trivial unit tests, TestA
* and TestB, and use macros to define associated UnitTestRunner 
* subclasses, TEST_RUNNER(TestA) and TEST_RUNNER(TestB). The preprocessor
* macros in the main program then define a class CompositeExample that 
* is derived from CompositeTestRunner, which contains instances of the
* TEST_RUNNER(TestA) and TEST_RUNNER(TestB).  Calling the run() method 
* of the CompositeExample then runs the all of the tests defined in 
* TestA and TestB.
*/

/**
* Trivial UnitTest A.
*/
class TestA : public UnitTest
{

public:

   void test1() 
   { 
      printMethod(TEST_FUNC);
      TEST_ASSERT(eq(1, 2));
   }

   void test2() 
   { 
      printMethod(TEST_FUNC);
      TEST_ASSERT(false);
   }

   void test3() 
   { 
      printMethod(TEST_FUNC);
      TEST_ASSERT(eq(2, 6/3));
   }

};

TEST_BEGIN(TestA)
TEST_ADD(TestA, test1)
TEST_ADD(TestA, test2)
TEST_ADD(TestA, test3)
TEST_END(TestA)

/**
* Trivial UnitTest B.
*/
class TestB : public UnitTest
{

public:

   void test1() 
   { 
      printMethod(TEST_FUNC);
      TEST_ASSERT(eq(2, 9/3));
   }

   void test2() 
   { 
      printMethod(TEST_FUNC);
      TEST_ASSERT(eq(2, 4/2));
   }

};

TEST_BEGIN(TestB)
TEST_ADD(TestB, test1)
TEST_ADD(TestB, test2)
TEST_END(TestB)

int main() 
{

   TEST_COMPOSITE_BEGIN(CompositeExample)
   TEST_COMPOSITE_ADD_UNIT(TestA)
   TEST_COMPOSITE_ADD_UNIT(TestB)
   TEST_COMPOSITE_END
   
   CompositeExample runner;
   runner.run();

}

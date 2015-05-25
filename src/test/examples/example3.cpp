#define TEST_MPI

#include "UnitTest.h"
#include "UnitTestRunner.h"
#include <iostream>


/**
* Trivial example of UnitTest use for parallel MPI job.
*/

/**
* Trivial subclass of UnitTest for an MPI job.
*/
class TestA : public UnitTest
{

public:

   void test1() { 
      printMethod(TEST_FUNC);
      if (mpiRank() == 0) {
         TEST_ASSERT(true);
      } else {
         TEST_ASSERT(false);
      }
   }

   void test2() { 
      printMethod(TEST_FUNC);
      if (mpiRank() == 0) {
         TEST_ASSERT(false);
      } else {
         TEST_ASSERT(true);
      }
   }

   void test3() { 
      printMethod(TEST_FUNC);
      if (mpiRank() == 0) {
         TEST_ASSERT(false);
      } else {
         TEST_ASSERT(false);
      }
   }

   void test4() { 
      printMethod(TEST_FUNC);
      if (mpiRank() == 0) {
         TEST_ASSERT(true);
      } else {
         TEST_ASSERT(true);
      }
   }

};

/*
* Preprocessor macros to define an associated UnitTestRunner.
*/
TEST_BEGIN(TestA)
TEST_ADD(TestA, test1)
TEST_ADD(TestA, test2)
TEST_ADD(TestA, test3)
TEST_ADD(TestA, test4)
TEST_END(TestA)

int main() 
{ 
   MPI::Init();
   TestA_Runner test;
   test.run();
   MPI::Finalize();
}

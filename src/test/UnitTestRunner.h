#ifndef UNIT_TEST_RUNNER_H
#define UNIT_TEST_RUNNER_H

/*
* Simpatico - Simulation Package for Polymeric and Molecular Liquids
*
* Copyright 2010 - 2014, The Regents of the University of Minnesota
* Distributed under the terms of the GNU General Public License.
*/

#include "TestRunner.h"
#include "TestException.h"
#include <vector>

/**
* Template for a TestRunner that runs test methods of an associated UnitTest.
*
* A instance of UnitTestRunner<MyTest> holds an array of pointers to all of
* the test methods of a class MyTest that is a subclass of UnitTest. Each such
* test method must return void and take zero parameters. The addTestMethod() 
* method is used to register a test method with the UnitTestRunner instantiation, 
* by adding a pointer to a test method to this array.  The run() method runs 
* all of the registered test methods in sequence.

* To run a set of unit tests one must: 
*
*  - Define a subclass of UnitTest,
*  - Define an associated subclass of UnitTestRunner, 
*  - Construct a UnitTestRunner object and call run().
*
* The boilerplate code required to define a UnitTestRunner class may be 
* simplified by using set preprocessor macros that are defined at the end of 
* this file.
*
* Here is an example of the code to to define a subclass of UnitTestRunner<MyTest>,
* associated with a subclass MyTest of UnitTest, and then run all of its test
* methods, written without using any preprocessor macros:
* \code
*
* // Define a UnitTest class
* class MyTest : public UnitTest {
* public:
*
*    test1()
*    { .... }
*
*    test2
*    { ...  }
*
* };
*
* // Define a UnitTestRunner associated with MyTest
* class MyTest_Runner : public UnitTestRunner<MyTest> {
* public:
*
*    MyTest_Runner(){
*       addTestMethod(&MyTest::test1);
*       addTestMethod(&MyTest::test2);
*    }
*
* }
*
* // Run the tests.
* MyTest_Runner runner;
* runner.run();
*
* \endcode
* Note that, by convention:
*
*   - We defined a subclass of UnitTestRunner<MyTest>, called MyTest_Runner.
*   - All test methods of MyTest are registered in the MyTest_Runner constructor.
*
* Calling the run() method of MyTest_Runner will then run all of the tests.
*
* The following series of preprocessor macros may be used to generate the 
* definition of the MyTest_Runner class in the above example, and to create 
* an instance of this class:
* \code
* 
* TEST_BEGIN(MyTest)
* TEST_ADD(MyTest, test1)
* TEST_ADD(MyTest, test2)
* TEST_END(MyTest)
*
* TEST_RUNNER(MyTest) runner;
* runner.run();
*
* \endcode
* The macro TEST_BEGIN(TestClass) generates the beginning of the class
* definition for subclass MyTest_Runner of UnitTestRunner<TestClass>. 
* The TEST_ADD(TestClass, Method) adds a specified method of the associated
* class TestClass to the constructor of the new UnitTestRunner class. The
* TEST_END macro closes both the constructor definition and the class 
* definition. After expansion, the resulting code is completely equivalent
* to that given in the previous example, after the definition of MyTest.
* 
* The name of the UnitTestRunner class created by these preprocessor
* macros is created by appending the standard suffix "_Runner" to the 
* name of the unit test class. Thus, in the above example, the TestRunner 
* subclass is named MyTest_Runner. This TestRunner subclass name may
* be referred directly, using this name, or by using the preprocessor 
* macro TEST_RUNNER(TestClass), which expands to the name of the test
* runner class, e.g., to TestClass_Runner. In the above example, this 
* macro is used as a class name to instantiate an instance of the of
* required test runner. 
*
* \ingroup Test_Module
*/
template <class UnitTestClass>
class UnitTestRunner : public TestRunner
{

public:

   using TestRunner::nFailure;
   using TestRunner::isIoProcessor;

   /**
   * Pointer to a test method of the associated UnitTest class.
   */
   typedef void (UnitTestClass::*MethodPtr)();

   /**
   * Constructor.
   */
   UnitTestRunner();

   /**
   * Destructor.
   */
   ~UnitTestRunner();

   // Use compiler generated destructor.

   /**
   * Register a test method of the associated unit test class.
   */
   void addTestMethod(MethodPtr methodPtr);

   /**
   * Return the number of registered test methods.
   */
   int nTestMethod();

   /**
   * Run test method number i.
   *
   * \param i index of test method
   */
   void method(unsigned int i);

   /**
   * Run all registered test methods in the order added.
   */
   virtual int run();

private:

   std::vector<MethodPtr> methodPtrs_;

   #ifdef TEST_MPI
   std::vector<int> results_;
   #endif 

};


/*
* Constructor.
*/
template <class UnitTestClass>
UnitTestRunner<UnitTestClass>::UnitTestRunner()
 : TestRunner()
{
   #ifdef TEST_MPI
   if (isIoProcessor()) {
      results_.reserve(mpiSize());
      for (int i=0; i < mpiSize(); ++i) {
         results_.push_back(false);
      }
   }
   #endif
}

/*
* Destructor.
*/
template <class UnitTestClass>
UnitTestRunner<UnitTestClass>::~UnitTestRunner()
{}

/*
* Register a test method of the associated unit test class.
*/
template <class UnitTestClass>
void UnitTestRunner<UnitTestClass>::addTestMethod(MethodPtr methodPtr)
{  methodPtrs_.push_back(methodPtr); }

/*
* Return the number of registered test methods.
*/
template <class UnitTestClass>
int UnitTestRunner<UnitTestClass>::nTestMethod()
{  return methodPtrs_.size(); } 

/*
* Run test method number i.
*
* \param i index of test method
*/
template <class UnitTestClass>
void UnitTestRunner<UnitTestClass>::method(unsigned int i)
{
   UnitTestClass testCase;
   #ifdef TEST_MPI
   TestException exception;
   int result;
   #endif

   testCase.setFilePrefix(filePrefix());

   // Run test method (try / catch)
   try {
      testCase.setUp();
      (testCase.*methodPtrs_[i])();
      #ifndef TEST_MPI
      std::cout << ".";
      recordSuccess();
      #else
      result = 1;
      #endif
   } catch (TestException &e) {
      #ifndef TEST_MPI
      std::cout << std::endl;
      std::cout << " Failure " << std::endl << std::endl;
      std::cout << e.message() << std::endl;
      std::cout << ".";
      recordFailure();
      #else
      result = 0;
      exception = e;
      #endif
   }
   testCase.tearDown();

   #ifdef TEST_MPI
   // Process MPI Tests from all processors
   MPI::COMM_WORLD.Barrier();
   if (mpiRank() == 0) {
      results_[0] = result; // Result on processor 0
      if (results_[0] == 0) {
         std::cout << std::endl;
         std::cout << " Failure  on Processor 0" 
                   << std::endl << std::endl;
         std::cout << exception.message() << std::endl;
         std::cout.flush();
      }
      for (int i=1; i < mpiSize(); ++i) {
         // Receive result results_[i] of test on processor i.
         MPI::COMM_WORLD.Recv(&(results_[i]), 1, MPI_INT, i, i);
         // If the test failed on processor i
         if (results_[i] == 0) {
            result = 0; // fails (==0) if failure on any processor
            // Send permission to print failure on processor i
            MPI::COMM_WORLD.Send(&(results_[i]), 1, MPI_INT, i, mpiSize() + i);
            // Receive confirmation that processor i completed printing
            MPI::COMM_WORLD.Recv(&(results_[i]), 1, MPI_INT, i, 2*mpiSize() + i);
         }
      }
      // Record success on master (rank == 0) iff success on all processors
      if (result) {
         recordSuccess();
      } else {
         recordFailure();
      }
      std::cout << ".";
      std::cout.flush();
   } else {   // Slave node
      // Send result of test on this processor to master, tag = mpiRank().
      MPI::COMM_WORLD.Send(&result, 1, MPI_INT, 0, mpiRank());
      // If test failed on this processor
      if (result == 0) {
         // Receive permission to print failure statement
         MPI::COMM_WORLD.Recv(&result, 1, MPI_INT, 0, 
                              mpiSize() + mpiRank());
         std::cout.flush();
         std::cout << std::endl;
         std::cout << " Failure  on Processor " << mpiRank() 
                   << std::endl << std::endl;
         std::cout << exception.message() << std::endl;
         std::cout.flush();
         // Confirm completion of print.
         MPI::COMM_WORLD.Send(&result, 1, MPI_INT, 0, 
                              2*mpiSize() + mpiRank());
      }
   }
   MPI::COMM_WORLD.Barrier();
   #endif // ifdef TEST_MPI

}

/*
* Run all registered test methods in the order added.
*/
template <class UnitTestClass>
int UnitTestRunner<UnitTestClass>::run()
{
   for (unsigned int i = 0; i < methodPtrs_.size(); ++i) {
      method(i);
   }
   report();
   return nFailure();
}

// Preprocessor Macros -------------------------------------------------

/**
* Macro for name of the UnitTestRunner class associated with UnitTestClass.
*/
#define TEST_RUNNER(UnitTestClass) UnitTestClass##_Runner

/**
* Begin definition of class TEST_RUNNER(UnitTestClass).
*
* This macro generates code to open the class definition, and to open
* the definition of a default constructor.
*/
#define TEST_BEGIN(UnitTestClass) \
   class TEST_RUNNER(UnitTestClass) \
    : public UnitTestRunner<UnitTestClass> \
   { public: TEST_RUNNER(UnitTestClass)() { 

/**
* Macro to add a test method to TEST_RUNNER(UnitTestClass).
*/
#define TEST_ADD(UnitTestClass, Method) \
   addTestMethod(&UnitTestClass::Method);

/**
* Macro to end definition of a class TEST_RUNNER(UnitTestClass).
*
* This macro ends both the constructor and class definition.
*/
#define TEST_END(UnitTestClass) } }; 

#endif

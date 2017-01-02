#ifndef TEST_RUNNER_H
#define TEST_RUNNER_H

/*
* Simpatico - Simulation Package for Polymeric and Molecular Liquids
*
* Copyright 2010 - 2014, The Regents of the University of Minnesota
* Distributed under the terms of the GNU General Public License.
*/

#include <string>
#include <iostream>

#ifdef TEST_MPI
#include <mpi.h>
#endif

/**
* Abstract base class for classes that run tests.
*
* TestRunner is an abstract base class with two types of subclass:
* The UnitTestRunner class template defines a TestRunner that runs
* the tests for an associated UnitTest. A CompositeTestRunner runs 
* the tests for a sequence of other TestRunner objects, each of 
* which can be a UnitTestRunner or another CompositeTestRunner.
*
* An implementation of the pure virtual run() method of must run 
* all of the associated test methods, and records the number 
* nSuccess() of tests that succeed and the number nFailure() that 
* fail.  A test fails if it throws a TestException. Test methods 
* use the TEST_ASSERT(expr) macro to assert the truth of a logical
* exprression expr, which throws a TestException if expr is false.
* The implementation of run() for a UnitTestRunner runs each unit 
* test method of the associated UnitTest in a try-catch block and
* catches any thrown TestExceptions. The implementation of run 
* for a TestComposite calls the run() method for each of its
* children. 
* 
* Each TestRunner may have a parent TestRunner. The parent, if
* any, is always TestComposite. A TestComposite can have any
* number of children.  The recordFailure() and recordSuccess() 
* methods, which can be called by the run method, increment the 
* nSuccess or nFailure counters, and also call the corresponding 
* method of the parent, if any. The nSuccess and nFailure counters 
* for a TestComposite thereby keep track of the total number of 
* successful and failed tests run by all descendants. 
*
* Each TestRunner has a filePrefix string. This string is 
* prepended to the names of any files opened by the functions
* openInputFile() and openOutputFile() of the UnitTest subclass.
* The filePrefix is initialized to an empty string, and may be 
* modified by the virtual addFilePrefix() function. The default
* implementation of this function prepends a string argument to 
* the existing filePrefix. The implementation by TestComposite
* calls the addFilePrefix method of each of its children, thus
* allowing a common prefix to be added to the file paths used
* by all of its children.
* 
* \ingroup Test_Module
*/
class TestRunner
{

public:

   /**
   * Constructor.
   */
   TestRunner();

   /**
   * Destructor.
   */
   virtual ~TestRunner();

   /**
   * Run all tests.
   *
   * \return number of failures.
   */
   virtual int run() = 0;

   /**
   * Increment counter for failed tests, and that of parent (if any).
   */
   void recordFailure();

   /**
   * Increment counter for successful tests, and that of parent (if any).
   */
   void recordSuccess();

   /**
   * Set another TestRunner as the parent.
   *
   * \param parent parent CompositeTestRunner object
   */
   void setParent(TestRunner& parent);

   /**
   * Return the parent object, if any.
   */
   TestRunner& parent();

   /**
   * Does this object have a parent?
   */
   bool hasParent() const;

   /**
   * Return number of successful tests run.
   */
   int nSuccess() const;

   /**
   * Return number of failed tests run.
   */
   int nFailure() const;

   /**
   * If this object has no parent, report success and failure counters.
   */
   void report() const;

   /**
   * Is this the IO processor of an MPI communicator?
   */
   bool isIoProcessor() const;

   #ifdef TEST_MPI
   /**
   * Return the MPI rank in the communicator.
   */
   int mpiRank() const;

   /**
   * Return the size (number of processors) of the communicator.
   */
   int mpiSize() const;
   #endif

   /**
   * Prepend argument prefix to existing filePrefix.
   */
   virtual void addFilePrefix(const std::string& prefix);

   /**
   * Return file prefix by const reference.
   */
   const std::string& filePrefix() const;

protected:

   /// Prefix added to file names
   std::string  filePrefix_;

private:

   /// Pointer to a parent TestRunner (if any).
   TestRunner* parentPtr_;

   /// Total number of successful tests run.
   int  nSuccess_;

   /// Total number of failed tests run.
   int  nFailure_;

   /// Can this processor input and output data?
   /// This is always true when TEST_MPI is not defined.
   bool  isIoProcessor_;

   #ifdef TEST_MPI
   /// Rank of this processor within an MPI job.
   int  mpiRank_;

   /// Size of associated MPI communicator.
   int  mpiSize_;
   #endif

};

// Inline methods

/*
* Set another TestRunner as the parent.
*/
inline void TestRunner::setParent(TestRunner& parent)
{  parentPtr_ = &parent; }

/*
* Return the parent object, if any.
*/
inline TestRunner& TestRunner::parent()
{  return *parentPtr_; }

/*
* Does this object have a parent?
*/
inline bool TestRunner::hasParent() const
{  return (parentPtr_ != 0); }

/*
* Return number of successful tests run.
*/
inline int TestRunner::nSuccess() const
{  return nSuccess_; }

/*
* Return number of failed tests run.
*/
inline int TestRunner::nFailure() const
{  return nFailure_; }

/*
* Return file prefix by const reference.
*/
inline 
const std::string& TestRunner::filePrefix() const
{  return filePrefix_; }

/*
* Is this an Io processor? (always true without MPI)
*/
inline bool TestRunner::isIoProcessor() const
{  return isIoProcessor_; }

#ifdef TEST_MPI
inline int TestRunner::mpiRank() const
{  return mpiRank_; }

inline int TestRunner::mpiSize() const
{  return mpiSize_; } 
#endif

// Non-inline methods

/*
* Constructor.
*/
TestRunner::TestRunner()
 : parentPtr_(0),
   nSuccess_(0),
   nFailure_(0),
   isIoProcessor_(0)
   #ifdef TEST_MPI
   ,mpiRank_(0),
   mpiSize_(0)
   #endif
{
    #ifndef TEST_MPI
    isIoProcessor_ = true; 
    #else
    mpiRank_ = MPI::COMM_WORLD.Get_rank();
    mpiSize_ = MPI::COMM_WORLD.Get_size();
    if (mpiRank_ == 0) {
       isIoProcessor_ = true; 
    } else {
       isIoProcessor_ = false; 
    }
    #endif
}

/*
* Destructor.
*/
TestRunner::~TestRunner()
{}

/*
* Increment counter for failed tests, and that of parent (if any).
*/
void TestRunner::recordFailure() 
{
   if (isIoProcessor()) {
      ++nFailure_;
      if (hasParent()) {
         parent().recordFailure();
      }
   }
}

/*
* Increment counter for successful tests, and that of parent (if any).
*/
void TestRunner::recordSuccess()
{
   if (isIoProcessor()) {
      ++nSuccess_;
      if (hasParent()) {
         parent().recordSuccess();
      }
   }
}

/*
* If this object has no parent, report success and failure counters.
*/
void TestRunner::report() const
{ 
   if (!hasParent() && isIoProcessor()) {
      std::cout << std::endl;
      std::cout << nSuccess_ << "  successful tests  " << std::endl;
      std::cout << nFailure_ << "  failed tests  "    << std::endl;
      std::cout << std::endl;
   }
}

/*
* Prepend argument prefix to existing filePrefix (virtual).
*/
void TestRunner::addFilePrefix(const std::string& prefix) 
{
   std::string newPrefix = prefix;
   newPrefix += filePrefix_;
   filePrefix_ = newPrefix;
}

#endif

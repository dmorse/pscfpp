#ifndef COMPOSITE_TEST_RUNNER_H
#define COMPOSITE_TEST_RUNNER_H

/*
* Simpatico - Simulation Package for Polymeric and Molecular Liquids
*
* Copyright 2010 - 2014, The Regents of the University of Minnesota
* Distributed under the terms of the GNU General Public License.
*/

#include "TestRunner.h"
#include <vector>

/**
* A TestRunner comprised of one or more child TestRunners.
*
* \ingroup Test_Module
*/
class CompositeTestRunner : public TestRunner
{

public:

   // Default constructor.

   /**
   * Destructor.
   */
   virtual ~CompositeTestRunner();

   /**
   * Add an existing TestRunner as a child.
   *
   * Children added by this method are not destroyed by the parent
   * CompositeTestRunner destructor.
   *
   * \param child enclosed TestRunner object
   */
   void addChild(TestRunner& child);

   /**
   * Add a TestRunner as a child, and accept ownership.
   *
   * Children added by this method are owned by the parent CompositeTestRunner,
   * and so are destroyed by its destructor.  
   *
   * \param childPtr pointer to child TestRunner
   */
   void addChild(TestRunner* childPtr);

   /**
   * Add a TestRunner as a child, accept ownership, and set file prefix.
   *
   * Children added by this method are owned by the parent CompositeTestRunner,
   * and so are destroyed by its destructor.  
   *
   * \param childPtr  pointer to child TestRunner
   * \param prefix  prefix to append to file names in all descendants
   */
   void addChild(TestRunner* childPtr, const std::string& prefix);

   /**
   * Prepend argument prefix to existing filePrefix.
   *
   * \param prefix string to prepend to existing filePrefix.
   */
   virtual void addFilePrefix(const std::string& prefix);

   /**
   * Run all children in sequence, using depth-first recursion. 
   */
   virtual int run(); 

private:

   /// Vector of pointers to child TestRunner objects.
   std::vector<TestRunner*> children_;

   /// Vector of pointers to child TestRunner objects owned by this object.
   std::vector<TestRunner*> ownedChildren_;

};

// Function definitions

/*
* Destructor.
*/
CompositeTestRunner::~CompositeTestRunner()
{
   unsigned int i; 
   for (i = 0; i < children_.size(); ++i) {
      delete ownedChildren_[i]; 
   }
}

/*
* Add an existing TestRunner as a child.
*/
void CompositeTestRunner::addChild(TestRunner& child)
{
   children_.push_back(&child); 
   child.setParent(*this);
}

/*
* Add a TestRunner as a child, and accept ownership.
*/
void CompositeTestRunner::addChild(TestRunner* childPtr)
{
   children_.push_back(childPtr);
   ownedChildren_.push_back(childPtr);
   childPtr->setParent(*this);
}

/*
* Add a TestRunner as a child, accept ownership, and set file prefix.
*/
void CompositeTestRunner::addChild(TestRunner* childPtr, 
                                   const std::string& prefix)
{
   addChild(childPtr);
   childPtr->addFilePrefix(prefix);
}

/*
* Prepend argument prefix to existing filePrefix.
*/
void CompositeTestRunner::addFilePrefix(const std::string& prefix) 
{
   TestRunner::addFilePrefix(prefix);
   for (unsigned int i = 0; i < children_.size(); ++i) {
      children_[i]->addFilePrefix(prefix); 
   }
}

/*
* Run all children in sequence, using depth-first recursion. 
*/
int CompositeTestRunner::run() 
{
   for (unsigned int i = 0; i < children_.size(); ++i) {
      children_[i]->run(); 
   }
   report();
   return nFailure();
}

// Preprocessor macros

/**
* Macro to open a TestComposite class definition.
*
* This macro opens both the class definition and a constructor.
*/
#define TEST_COMPOSITE_BEGIN(CompositeClass) \
   class CompositeClass : public CompositeTestRunner  { public: \
      CompositeClass () {

/**
* Macro to add a UnitTest subclass to a CompositeTestRunner constructor.
*
* The parameter UnitTestClass is the name of the UnitTest subclass.
* An instance of corresponding UnitTestRunner is instantiated internally.
*
*/
#define TEST_COMPOSITE_ADD_UNIT(UnitTestClass) \
         addChild(new TEST_RUNNER(UnitTestClass)); 

/**
* Macro to add a UnitTest subclass to a CompositeTestRunner constructor.
*
* The parameter UnitTestClass is the name of the UnitTest subclass.
* An instance of corresponding UnitTestRunner is instantiated internally.
*
*/
#define TEST_COMPOSITE_ADD_CHILD(TestRunner, Prefix) \
         addChild(new TestRunner, Prefix);

/**
* Macro to close a TestComposite class definition.
*
* This macro closes both the constructor and class definitions.
*/
#define TEST_COMPOSITE_END } };

#endif

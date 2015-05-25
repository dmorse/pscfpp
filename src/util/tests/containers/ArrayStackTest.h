#ifndef ARRAY_STACK_TEST_H
#define ARRAY_STACK_TEST_H

#include <test/UnitTest.h>
#include <test/UnitTestRunner.h>
#include <util/containers/DArray.h>
#include <util/containers/ArrayStack.h>

using namespace Util;

class ArrayStackTest : public UnitTest 
{

private:

   const static int capacity = 10;

   typedef int Data;

   DArray<Data>* arrayPtr; 
   ArrayStack<Data>* stackPtr;
   
   DArray<Data>& array(){ return *arrayPtr; }
   ArrayStack<Data>& stack(){ return *stackPtr; }

public:

   void setUp();
   void tearDown();
  
   void testPush();
   void testPop();

};

void ArrayStackTest::setUp()
{
   TEST_ASSERT(Memory::total() == 0);
   arrayPtr = new DArray<Data>;
   stackPtr = new ArrayStack<Data>;

   array().allocate(capacity);
   for (int i=0; i < capacity; i++ ) {
      array()[i] = (i+1)*10 + 1;
   }
   stack().allocate(capacity);
}

void ArrayStackTest::tearDown()
{}

  
void ArrayStackTest::testPush()
{
   printMethod(TEST_FUNC);
   stack().push(array()[8]);
   TEST_ASSERT(stack().isValid());
   TEST_ASSERT(stack().size() == 1);

   stack().push(array()[4]);
   stack().push(array()[3]);
   TEST_ASSERT(stack().isValid());
   TEST_ASSERT(stack().size() == 3);

   delete arrayPtr;
   delete stackPtr;
   TEST_ASSERT(Memory::total() == 0);
} 

void ArrayStackTest::testPop()
{
   Data top;

   printMethod(TEST_FUNC);
   stack().push(array()[8]);
   stack().push(array()[4]);
   stack().push(array()[3]);

   TEST_ASSERT(stack().isValid());
   TEST_ASSERT(stack().size() == 3);

   top = stack().pop();
   TEST_ASSERT(top == array()[3]);
   TEST_ASSERT(stack().isValid());
   TEST_ASSERT(stack().size() == 2);

   top = stack().pop();
   TEST_ASSERT(top == array()[4]);
   TEST_ASSERT(stack().isValid());
   TEST_ASSERT(stack().size() == 1);

   stack().push(array()[7]);
   TEST_ASSERT(stack().isValid());
   TEST_ASSERT(stack().size() == 2);

   top = stack().pop();
   TEST_ASSERT(top == array()[7]);
   TEST_ASSERT(stack().isValid());
   TEST_ASSERT(stack().size() == 1);

   top = stack().peek();
   TEST_ASSERT(top == array()[8]);

   top = stack().pop();
   TEST_ASSERT(top == array()[8]);
   TEST_ASSERT(stack().isValid());
   TEST_ASSERT(stack().size() == 0);

   delete arrayPtr;
   delete stackPtr;
   TEST_ASSERT(Memory::total() == 0);
}

TEST_BEGIN(ArrayStackTest)
TEST_ADD(ArrayStackTest, testPush)
TEST_ADD(ArrayStackTest, testPop)
TEST_END(ArrayStackTest)

#endif

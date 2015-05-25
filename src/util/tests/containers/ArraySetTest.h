#ifndef ARRAY_SET_TEST_H
#define ARRAY_SET_TEST_H

#include <test/UnitTest.h>
#include <test/UnitTestRunner.h>

#include <util/containers/DArray.h>
#include <util/containers/ArraySet.h>
#include <util/containers/PArrayIterator.h>

using namespace Util;

class ArraySetTest : public UnitTest 
{

private:

   typedef int Data;

   const static int    capacity = 10;

   DArray<Data>* arrayPtr; 
   ArraySet<Data>* setPtr;
   PArrayIterator<Data> iterator;

   DArray<Data>& array() { return *arrayPtr; }
   ArraySet<Data>& set() { return *setPtr; }

   int memory_;
   
public:

   void setUp();
   void tearDown();
  
   void testAppend();
   void testRemove();
   void testPop();
   void testIterator();

};

void ArraySetTest::setUp()
{
   arrayPtr = new DArray<Data>;
   setPtr = new ArraySet<Data>;

   memory_ = Memory::total();
   array().allocate(capacity);
   TEST_ASSERT(Memory::total() == memory_ + capacity*((int)sizeof(Data)) );
   for (int i = 0; i < capacity; ++i) {
      array()[i] = (i+1)*10 + 1;
   }
   set().allocate(array());
}

void ArraySetTest::tearDown()
{}

  
void ArraySetTest::testAppend()
{
   printMethod(TEST_FUNC);
   {
      set().append(array()[8]);
      TEST_ASSERT(set().isValid());
   
      set().append(array()[4]);
      set().append(array()[3]);
      //set().dump();
      TEST_ASSERT(set().isValid());
      TEST_ASSERT(set().size() == 3);
      TEST_ASSERT(set()[0] == array()[8]);
      TEST_ASSERT(set()[1] == array()[4]);
      TEST_ASSERT(set()[2] == array()[3]);
   
      delete arrayPtr;
      delete setPtr;
   }
   TEST_ASSERT(Memory::total() == memory_);
} 

void ArraySetTest::testRemove()
{
   printMethod(TEST_FUNC);
   {
      set().append(array()[8]);
      set().append(array()[4]);
      set().append(array()[3]);
      set().append(array()[5]);
      TEST_ASSERT(set()[0] == array()[8]);
      TEST_ASSERT(set()[1] == array()[4]);
      TEST_ASSERT(set()[2] == array()[3]);
      TEST_ASSERT(set()[3] == array()[5]);
   
      TEST_ASSERT(set().isValid());
      TEST_ASSERT(set().size() == 4);
   
      set().remove(array()[4]);
      //set().dump();
      TEST_ASSERT(set().isValid());
      TEST_ASSERT(set().size() == 3);
      TEST_ASSERT(set()[0] == array()[8]);
      TEST_ASSERT(set()[1] == array()[5]);
      TEST_ASSERT(set()[2] == array()[3]);
   
      set().remove(array()[3]);
      //set().dump();
      TEST_ASSERT(set().isValid());
      TEST_ASSERT(set().size() == 2);
      TEST_ASSERT(set()[0] == array()[8]);
      TEST_ASSERT(set()[1] == array()[5]);
   
      set().remove(array()[8]);
      //set().dump();
      TEST_ASSERT(set().isValid());
      TEST_ASSERT(set().size() == 1);
      TEST_ASSERT(set()[0] == array()[5]);
   
      delete arrayPtr;
      delete setPtr;
   }
   TEST_ASSERT(Memory::total() == memory_);
}

void ArraySetTest::testPop()
{
   printMethod(TEST_FUNC);
   {
      set().append(array()[8]);
      set().append(array()[4]);
      set().append(array()[2]);
      set().append(array()[3]);
      set().append(array()[5]);
   
      set().remove(array()[2]);
   
      TEST_ASSERT(set().isValid());
      TEST_ASSERT(set().size() == 4);
      TEST_ASSERT(set()[3] == array()[3]);
   
      Data* ptr; 
      ptr = &set().pop();
      TEST_ASSERT(set().isValid());
      TEST_ASSERT(set().size() == 3);
      TEST_ASSERT(ptr = &array()[3]);
      
      delete arrayPtr;
      delete setPtr;
   }
   TEST_ASSERT(Memory::total() == memory_);
}

void ArraySetTest::testIterator()
{
   printMethod(TEST_FUNC);
   {
      set().append(array()[8]);
      set().append(array()[4]);
      set().append(array()[2]);
      set().append(array()[3]);
      set().append(array()[5]);
      set().remove(array()[2]);
   
      TEST_ASSERT(set().isValid());
      TEST_ASSERT(set().size() == 4);
   
      set().begin(iterator);
   
      TEST_ASSERT(!iterator.isEnd());
      TEST_ASSERT(*iterator == array()[8]);
      ++iterator;
   
      TEST_ASSERT(!iterator.isEnd());
      TEST_ASSERT(*iterator == array()[4]);
      ++iterator;
   
      TEST_ASSERT(!iterator.isEnd());
      TEST_ASSERT(*iterator == array()[5]);
      ++iterator;
   
      TEST_ASSERT(!iterator.isEnd());
      TEST_ASSERT(*iterator == array()[3]);
      ++iterator;
   
      TEST_ASSERT(iterator.isEnd());
   
      delete arrayPtr;
      delete setPtr;
   }
   TEST_ASSERT(Memory::total() == memory_);
}

TEST_BEGIN(ArraySetTest)
TEST_ADD(ArraySetTest, testAppend)
TEST_ADD(ArraySetTest, testRemove)
TEST_ADD(ArraySetTest, testPop)
TEST_ADD(ArraySetTest, testIterator)
TEST_END(ArraySetTest)

#endif

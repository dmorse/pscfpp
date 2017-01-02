#ifndef S_SET_TEST_H
#define S_SET_TEST_H

#include <test/UnitTest.h>
#include <test/UnitTestRunner.h>

#include <util/containers/SSet.h>
#include <util/containers/PArrayIterator.h>
#include <util/containers/DArray.h>

using namespace Util;

class SSetTest : public UnitTest 
{

public:

   void setUp();
   void tearDown();
  
   void testAppend();
   void testRemove();
   void testIterator();

private:

   const static int Capacity = 10;

   typedef int Data;

   SSet<Data, Capacity> set;
   PArrayIterator<Data> iterator;
   DArray<Data>         array;
   
};

void SSetTest::setUp()
{
   array.allocate(Capacity);
   for (int i = 0; i < Capacity; ++i) {
      array[i] = (i+1)*10 + 1;
   }
}

void SSetTest::tearDown()
{}

  
void SSetTest::testAppend()
{
   printMethod(TEST_FUNC);
   set.append(array[8]);
   //TEST_ASSERT(set.isValid());

   set.append(array[4]);
   set.append(array[3]);
   //TEST_ASSERT(set.isValid());
   TEST_ASSERT(set.size() == 3);
   TEST_ASSERT(set[0] == array[8]);
   TEST_ASSERT(set[1] == array[4]);
   TEST_ASSERT(set[2] == array[3]);
   TEST_ASSERT(set.isElement(array[3]));
   TEST_ASSERT(set.isElement(array[4]));
   TEST_ASSERT(set.isElement(array[8]));
   TEST_ASSERT(!set.isElement(array[1]));
   TEST_ASSERT(!set.isElement(array[2]));
   TEST_ASSERT(!set.isElement(array[5]));
   TEST_ASSERT(!set.isElement(array[6]));
   TEST_ASSERT(!set.isElement(array[7]));

} 

void SSetTest::testRemove()
{
   printMethod(TEST_FUNC);
   set.append(array[8]);
   set.append(array[4]);
   set.append(array[3]);
   set.append(array[5]);
   TEST_ASSERT(set[0] == array[8]);
   TEST_ASSERT(set[1] == array[4]);
   TEST_ASSERT(set[2] == array[3]);
   TEST_ASSERT(set[3] == array[5]);

   //TEST_ASSERT(set.isValid());
   TEST_ASSERT(set.size() == 4);

   set.remove(array[4]);
   //TEST_ASSERT(set.isValid());
   TEST_ASSERT(set.size() == 3);
   TEST_ASSERT(set[0] == array[8]);
   TEST_ASSERT(set[1] == array[5]);
   TEST_ASSERT(set[2] == array[3]);

   set.remove(array[3]);
   //TEST_ASSERT(set.isValid());
   TEST_ASSERT(set.size() == 2);
   TEST_ASSERT(set[0] == array[8]);
   TEST_ASSERT(set[1] == array[5]);

   set.remove(array[8]);
   //TEST_ASSERT(set.isValid());
   TEST_ASSERT(set.size() == 1);
   TEST_ASSERT(set[0] == array[5]);

}

void SSetTest::testIterator()
{
   printMethod(TEST_FUNC);
   set.append(array[8]);
   set.append(array[4]);
   set.append(array[2]);
   set.append(array[3]);
   set.append(array[5]);
   set.remove(array[2]);

   //TEST_ASSERT(set.isValid());
   TEST_ASSERT(set.size() == 4);

   set.begin(iterator);

   TEST_ASSERT(iterator.notEnd());
   TEST_ASSERT(*iterator == array[8]);
   ++iterator;

   TEST_ASSERT(iterator.notEnd());
   TEST_ASSERT(*iterator == array[4]);
   ++iterator;

   TEST_ASSERT(iterator.notEnd());
   TEST_ASSERT(*iterator == array[5]);
   ++iterator;

   TEST_ASSERT(iterator.notEnd());
   TEST_ASSERT(*iterator == array[3]);
   ++iterator;

   TEST_ASSERT(iterator.isEnd());

}

TEST_BEGIN(SSetTest)
TEST_ADD(SSetTest, testAppend)
TEST_ADD(SSetTest, testRemove)
TEST_ADD(SSetTest, testIterator)
TEST_END(SSetTest)

#endif

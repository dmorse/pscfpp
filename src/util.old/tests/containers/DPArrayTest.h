#ifndef DP_ARRAY_TEST_H
#define DP_ARRAY_TEST_H

#include <test/UnitTest.h>
#include <test/UnitTestRunner.h>

#include <util/containers/DArray.h>
#include <util/containers/DPArray.h>
#include <util/containers/PArrayIterator.h>

using namespace Util;

class DPArrayTest : public UnitTest 
{

private:

   const static int    capacity = 10;

   typedef int Data;

   DArray<Data>         array; 
   DPArray<Data>        parray;
   PArrayIterator<Data> iterator;
   
public:

   void setUp();
   void tearDown();
  
   void testAppend();
   void testModify();
   void testIterator();

};

void DPArrayTest::setUp()
{
   TEST_ASSERT(Memory::total() == 0);
   array.allocate(capacity);
   parray.allocate(capacity);
   for (int i=0; i < capacity; i++ ) {
      array[i] = (i+1)*10 + 1;
   }
}

void DPArrayTest::tearDown()
{}

  
void DPArrayTest::testAppend()
{
   printMethod(TEST_FUNC);
   parray.append(array[8]);
   parray.append(array[4]);
   parray.append(array[3]);

   TEST_ASSERT(parray.size() == 3);
   TEST_ASSERT(parray[0] == array[8]);
   TEST_ASSERT(parray[1] == array[4]);
   TEST_ASSERT(parray[2] == array[3]);

} 

void DPArrayTest::testModify()
{
   printMethod(TEST_FUNC);
   parray.append(array[8]);
   parray.append(array[4]);
   parray.append(array[3]);
   TEST_ASSERT(parray.size() == 3);

   try {
      parray[1] = 13;
   } catch (Exception e) {
      TEST_ASSERT(0);
   }
   TEST_ASSERT(parray[1] == 13);
   TEST_ASSERT(array[4]  == 13);
   TEST_ASSERT(parray[2] == array[3]);

} 

void DPArrayTest::testIterator()
{
   printMethod(TEST_FUNC);
   parray.append(array[8]);
   parray.append(array[4]);
   parray.append(array[3]);
   parray.append(array[5]);
   TEST_ASSERT(parray.size() == 4);

   parray.begin(iterator);

   TEST_ASSERT(!iterator.isEnd());
   TEST_ASSERT(iterator.notEnd());
   TEST_ASSERT(*iterator == array[8]);
   ++iterator;

   TEST_ASSERT(!iterator.isEnd());
   TEST_ASSERT(iterator.notEnd());
   TEST_ASSERT(*iterator == array[4]);
   ++iterator;

   TEST_ASSERT(!iterator.isEnd());
   TEST_ASSERT(iterator.notEnd());
   TEST_ASSERT(*iterator == array[3]);
   ++iterator;

   TEST_ASSERT(!iterator.isEnd());
   TEST_ASSERT(iterator.notEnd());
   TEST_ASSERT(*iterator == array[5]);
   ++iterator;

   TEST_ASSERT(iterator.isEnd());
   TEST_ASSERT(!iterator.notEnd());

}

TEST_BEGIN(DPArrayTest)
TEST_ADD(DPArrayTest, testAppend)
TEST_ADD(DPArrayTest, testModify)
TEST_ADD(DPArrayTest, testIterator)
TEST_END(DPArrayTest)

#endif

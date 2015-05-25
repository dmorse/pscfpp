#ifndef FP_ARRAY_TEST_H
#define FP_ARRAY_TEST_H

#include <test/UnitTest.h>
#include <test/UnitTestRunner.h>

#include <util/containers/DArray.h>
#include <util/containers/FPArray.h>
#include <util/containers/PArrayIterator.h>

using namespace Util;

class FPArrayTest : public UnitTest 
{

private:

   const static int    Capacity = 10;

   typedef int Data;

   DArray<Data>                array; 
   FPArray<Data,Capacity>      parray;
   PArrayIterator<Data>        iterator;
   ConstPArrayIterator<Data>   constIterator;
   
public:

   void setUp();
   void tearDown();
  
   void testAppend();
   void testModify();
   void testIterator();
   void testConstIterator();

};

void FPArrayTest::setUp()
{
   array.allocate(Capacity);
   for (int i=0; i < Capacity; i++ ) {
      array[i] = (i+1)*10 + 1;
   }

}

void FPArrayTest::tearDown()
{}

  
void FPArrayTest::testAppend()
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

void FPArrayTest::testModify()
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

void FPArrayTest::testIterator()
{
   printMethod(TEST_FUNC);
   parray.append(array[8]);
   parray.append(array[4]);
   parray.append(array[3]);
   parray.append(array[5]);
   TEST_ASSERT(parray.size() == 4);

   parray.begin(iterator);

   TEST_ASSERT(iterator.notEnd());
   TEST_ASSERT(*iterator == array[8]);
   ++iterator;

   TEST_ASSERT(iterator.notEnd());
   TEST_ASSERT(*iterator == array[4]);
   ++iterator;

   TEST_ASSERT(iterator.notEnd());
   TEST_ASSERT(*iterator == array[3]);
   ++iterator;

   TEST_ASSERT(iterator.notEnd());
   TEST_ASSERT(*iterator == array[5]);
   ++iterator;

   TEST_ASSERT(iterator.isEnd());

}

void FPArrayTest::testConstIterator()
{
   printMethod(TEST_FUNC);
   parray.append(array[8]);
   parray.append(array[4]);
   parray.append(array[3]);
   parray.append(array[5]);
   TEST_ASSERT(parray.size() == 4);

   // Make const reference to parray
   const FPArray<Data, Capacity>& carray = parray;

   carray.begin(constIterator);

   TEST_ASSERT(constIterator.notEnd());
   TEST_ASSERT(*constIterator == array[8]);
   ++constIterator;

   TEST_ASSERT(constIterator.notEnd());
   TEST_ASSERT(*constIterator == array[4]);
   ++constIterator;

   TEST_ASSERT(constIterator.notEnd());
   TEST_ASSERT(*constIterator == array[3]);
   ++constIterator;

   TEST_ASSERT(constIterator.notEnd());
   TEST_ASSERT(*constIterator == array[5]);
   ++constIterator;

   TEST_ASSERT(constIterator.isEnd());

}

TEST_BEGIN(FPArrayTest)
TEST_ADD(FPArrayTest, testAppend)
TEST_ADD(FPArrayTest, testModify)
TEST_ADD(FPArrayTest, testIterator)
TEST_ADD(FPArrayTest, testConstIterator)
TEST_END(FPArrayTest)

#endif

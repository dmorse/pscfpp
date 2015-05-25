#ifndef LIST_ARRAY_TEST_H
#define LIST_ARRAY_TEST_H

#include <test/UnitTest.h>
#include <test/UnitTestRunner.h>

#include <util/containers/List.h>
#include <util/containers/ListArray.h>

using namespace Util;

class ListArrayTest : public UnitTest 
{

private:

   const static int    capacity = 20;
   const static int    nList    =  1;

   typedef int Data;

   ListArray<Data>  listArray;
   
public:

   void setUp() {}
   void tearDown() {}
  
   void testAllocate();
   //void testPushBack();

};



void ListArrayTest::testAllocate()
{
   printMethod(TEST_FUNC);
   listArray.allocate(capacity, nList);

   for (int i=0; i < capacity; i++) {
      TEST_ASSERT(listArray[i].id() == i);
      TEST_ASSERT(&(listArray.node(i).data()) == &(listArray[i]));
   }
   TEST_ASSERT(listArray.isValid() );

   List<Data> &list = listArray.list(0);
   list.isValid();

} 

/*
void ListArrayTest::testPushBack()
{
   printMethod(TEST_FUNC);
   list.initialize(nodes, n);

   TEST_ASSERT(list.size_ == 0 );
   TEST_ASSERT(list.front_ == 0 );
   TEST_ASSERT(list.back_  == 0 );
   TEST_ASSERT(int(list.lower_ - list.nodes_) == n);
   TEST_ASSERT(int(list.upper_ - list.nodes_) == -1);

   list.pushBack(nodes[5]);
   list.pushBack(nodes[8]);
   list.pushBack(nodes[3]);

   TEST_ASSERT(list.size_ == 3);
   TEST_ASSERT(int(list.front_ - list.nodes_) == 5);
   TEST_ASSERT(int(list.back_  - list.nodes_) == 3);
   TEST_ASSERT(int(list.lower_ - list.nodes_) == 3);
   TEST_ASSERT(int(list.upper_ - list.nodes_) == 8);
   TEST_ASSERT(list.isValid());

}
*/ 

TEST_BEGIN(ListArrayTest)
TEST_ADD(ListArrayTest, testAllocate)
   //TEST_ADD(ListArrayTest, testPushBack);
TEST_END(ListArrayTest)

#endif

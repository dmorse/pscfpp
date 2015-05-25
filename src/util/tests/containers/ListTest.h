#ifndef LIST_TEST_H
#define LIST_TEST_H

#include <test/UnitTest.h>
#include <test/UnitTestRunner.h>

#include <util/containers/Node.h>
#include <util/containers/List.h>

using namespace Util;

class ListTest : public UnitTest 
{

private:

   const static int      n = 10;

   typedef int Data;

   List<Data>  list;
   Node<Data>  nodes[n];

public:

   void setUp();
   void tearDown();
  
   void testInitialize();
   void testPushBack();
   void testPushBack2();
   void testPushFront();
   void testPopBack();
   void testInsertNext1();
   void testInsertNext2();
   void testInsertPrev1();
   void testInsertPrev2();
   void testRemove1();
   void testRemove2();
   void testRemove3();
   void testInsert1();
   void testInsert2();
   void testInsert3();
   void testInsert4();
   void testInsert5();
   void dumpList();

};

void ListTest::setUp()
{
   for (int i=1; i < n; i++) {
      nodes[i].data() = i*10 + 1;
   }
}

void ListTest::tearDown()
{}

  
void ListTest::testInitialize()
{
   printMethod(TEST_FUNC);
   list.initialize(nodes, n);
   TEST_ASSERT(list.isValid());
} 

void ListTest::testPushBack()
{
   printMethod(TEST_FUNC);
   list.initialize(nodes, n);

   TEST_ASSERT(list.size_ == 0);
   TEST_ASSERT(list.front_ == 0);
   TEST_ASSERT(list.back_  == 0);
   TEST_ASSERT(int(list.lower_ - list.nodes_) == n);
   TEST_ASSERT(int(list.upper_ - list.nodes_) == -1);

   list.pushBack(nodes[5]);
   list.pushBack(nodes[8]);
   list.pushBack(nodes[3]);
   TEST_ASSERT(list.isValid());

   TEST_ASSERT(list.size_ == 3);
   TEST_ASSERT(int(list.front_ - list.nodes_) == 5);
   TEST_ASSERT(int(list.back_  - list.nodes_) == 3);
   TEST_ASSERT(int(list.lower_ - list.nodes_) == 3);
   TEST_ASSERT(int(list.upper_ - list.nodes_) == 8);

   //printf("size_  = %10i\n",   list.size_);
   //printf("front_ = %10i\n", int(list.front_ - list.nodes_));
   //printf("back_  = %10i\n", int(list.back_  - list.nodes_));
   //printf("lower_ = %10i\n", int(list.lower_ - list.nodes_));
   //printf("upper_ = %10i\n", int(list.upper_ - list.nodes_));

   //Node<Data> node;
   //for (int i=0; i < n; i++) {
   //   node = list.nodes_[i];
   //}

} 


void ListTest::testPushFront()
{
   printMethod(TEST_FUNC);
   list.initialize(nodes, n);

   TEST_ASSERT(list.size_ == 0);
   TEST_ASSERT(list.front_ == 0);
   TEST_ASSERT(list.back_  == 0);
   TEST_ASSERT(int(list.lower_ - list.nodes_) == n);
   TEST_ASSERT(int(list.upper_ - list.nodes_) == -1);

   list.pushFront(nodes[5]);
   list.pushFront(nodes[8]);
   list.pushFront(nodes[3]);

   TEST_ASSERT(list.size_ == 3);
   TEST_ASSERT(int(list.front_ - list.nodes_) == 3);
   TEST_ASSERT(int(list.back_  - list.nodes_) == 5);
   TEST_ASSERT(int(list.lower_ - list.nodes_) == 3);
   TEST_ASSERT(int(list.upper_ - list.nodes_) == 8);
   TEST_ASSERT(list.isValid());

   //printf("size_  = %10i\n",   list.size_);
   //printf("front_ = %10i\n", int(list.front_ - list.nodes_));
   //printf("back_  = %10i\n", int(list.back_  - list.nodes_));
   //printf("lower_ = %10i\n", int(list.lower_ - list.nodes_));
   //printf("upper_ = %10i\n", int(list.upper_ - list.nodes_));

   //Node<Data> node;
   //for (int i=0; i < n; i++) {
   //   node = list.nodes_[i];
   //}

} 


void ListTest::testPopBack()
{
   printMethod(TEST_FUNC);
   list.initialize(nodes, n);

   list.pushBack(nodes[5]);
   list.pushBack(nodes[8]);
   list.pushBack(nodes[3]);

   list.popBack();

   TEST_ASSERT(list.size_ == 2);
   TEST_ASSERT(int(list.front_ - list.nodes_) == 5);
   TEST_ASSERT(int(list.back_  - list.nodes_) == 8);
   TEST_ASSERT(int(list.lower_ - list.nodes_) == 5);
   TEST_ASSERT(int(list.upper_ - list.nodes_) == 8);
   TEST_ASSERT(list.isValid());

   list.popBack();

   TEST_ASSERT(list.size_ == 1);
   TEST_ASSERT(int(list.front_ - list.nodes_) == 5);
   TEST_ASSERT(int(list.back_  - list.nodes_) == 5);
   TEST_ASSERT(int(list.lower_ - list.nodes_) == 5);
   TEST_ASSERT(int(list.upper_ - list.nodes_) == 5);
   TEST_ASSERT(list.isValid());

   list.popBack();

   TEST_ASSERT(list.size_ == 0);
   TEST_ASSERT(list.front_ == 0);
   TEST_ASSERT(list.back_  == 0);
   TEST_ASSERT(int(list.lower_ - list.nodes_) ==  n);
   TEST_ASSERT(int(list.upper_ - list.nodes_) == -1);
   TEST_ASSERT(list.isValid());

} 


void ListTest::testInsertNext1()
{
   printMethod(TEST_FUNC);
   list.initialize(nodes, n);

   list.pushBack(nodes[5]);
   list.pushBack(nodes[8]);
   list.pushBack(nodes[3]);

   list.insertNext(nodes[5], nodes[9]);

   TEST_ASSERT(list.size_ == 4);
   TEST_ASSERT(int(list.front_ - list.nodes_) == 5);
   TEST_ASSERT(int(list.back_  - list.nodes_) == 3);
   TEST_ASSERT(int(list.lower_ - list.nodes_) == 3);
   TEST_ASSERT(int(list.upper_ - list.nodes_) == 9);
   TEST_ASSERT(list.isValid());

}


void ListTest::testInsertNext2()
{
   printMethod(TEST_FUNC);
   list.initialize(nodes, n);

   list.pushBack(nodes[5]);
   list.pushBack(nodes[8]);
   list.pushBack(nodes[3]);

   list.insertNext(nodes[3], nodes[9]);

   TEST_ASSERT(list.size_ == 4);
   TEST_ASSERT(int(list.front_ - list.nodes_) == 5);
   TEST_ASSERT(int(list.back_  - list.nodes_) == 9);
   TEST_ASSERT(int(list.lower_ - list.nodes_) == 3);
   TEST_ASSERT(int(list.upper_ - list.nodes_) == 9);
   TEST_ASSERT(list.isValid());

}


void ListTest::testInsertPrev1()
{
   printMethod(TEST_FUNC);
   list.initialize(nodes, n);

   list.pushBack(nodes[5]);
   list.pushBack(nodes[8]);
   list.pushBack(nodes[3]);

   list.insertPrev(nodes[5], nodes[9]);

   TEST_ASSERT(list.size_ == 4);
   TEST_ASSERT(int(list.front_ - list.nodes_) == 9);
   TEST_ASSERT(int(list.back_  - list.nodes_) == 3);
   TEST_ASSERT(int(list.lower_ - list.nodes_) == 3);
   TEST_ASSERT(int(list.upper_ - list.nodes_) == 9);
   TEST_ASSERT(list.isValid());

}


void ListTest::testInsertPrev2()
{
   printMethod(TEST_FUNC);
   list.initialize(nodes, n);

   list.pushBack(nodes[5]);
   list.pushBack(nodes[8]);
   list.pushBack(nodes[3]);

   list.insertPrev(nodes[3], nodes[9]);

   TEST_ASSERT(list.size_ == 4);
   TEST_ASSERT(int(list.front_ - list.nodes_) == 5);
   TEST_ASSERT(int(list.back_  - list.nodes_) == 3);
   TEST_ASSERT(int(list.lower_ - list.nodes_) == 3);
   TEST_ASSERT(int(list.upper_ - list.nodes_) == 9);
   TEST_ASSERT(list.isValid());

}

void ListTest::testRemove1()
{
   printMethod(TEST_FUNC);
   list.initialize(nodes, n);

   list.pushBack(nodes[5]);
   list.pushBack(nodes[8]);
   list.pushBack(nodes[3]);

   list.remove(nodes[5]);

   TEST_ASSERT(list.isValid());
   TEST_ASSERT(list.size_ == 2);
   TEST_ASSERT(int(list.front_ - list.nodes_) == 8);
   TEST_ASSERT(int(list.back_  - list.nodes_) == 3);
   TEST_ASSERT(int(list.lower_ - list.nodes_) == 3);
   TEST_ASSERT(int(list.upper_ - list.nodes_) == 8);

}

void ListTest::testRemove2()
{
   printMethod(TEST_FUNC);
   list.initialize(nodes, n);

   list.pushBack(nodes[5]);
   list.pushBack(nodes[8]);
   list.pushBack(nodes[3]);

   list.remove(nodes[8]);

   TEST_ASSERT(list.isValid());
   TEST_ASSERT(list.size_ == 2);
   TEST_ASSERT(int(list.front_ - list.nodes_) == 5);
   TEST_ASSERT(int(list.back_  - list.nodes_) == 3);
   TEST_ASSERT(int(list.lower_ - list.nodes_) == 3);
   TEST_ASSERT(int(list.upper_ - list.nodes_) == 5);

}


void ListTest::testRemove3()
{
   printMethod(TEST_FUNC);
   list.initialize(nodes, n);

   list.pushBack(nodes[5]);
   list.pushBack(nodes[8]);
   list.pushBack(nodes[3]);

   list.remove(nodes[3]);

   TEST_ASSERT(list.isValid());
   TEST_ASSERT(list.size_ == 2);
   TEST_ASSERT(int(list.front_ - list.nodes_) == 5);
   TEST_ASSERT(int(list.back_  - list.nodes_) == 8);
   TEST_ASSERT(int(list.lower_ - list.nodes_) == 5);
   TEST_ASSERT(int(list.upper_ - list.nodes_) == 8);

}


void ListTest::testInsert1()
{
   printMethod(TEST_FUNC);
   list.initialize(nodes, n);

   list.pushBack(nodes[5]);
   list.pushBack(nodes[8]);
   list.pushBack(nodes[3]);

   list.insert(nodes[7]);

   TEST_ASSERT(list.size_ == 4);
   TEST_ASSERT(int(list.front_ - list.nodes_) == 5);
   TEST_ASSERT(int(list.back_  - list.nodes_) == 3);
   TEST_ASSERT(int(list.lower_ - list.nodes_) == 3);
   TEST_ASSERT(int(list.upper_ - list.nodes_) == 8);
   TEST_ASSERT(list.isValid());

}


void ListTest::testInsert2()
{
   printMethod(TEST_FUNC);
   list.initialize(nodes, n);

   list.pushBack(nodes[5]);
   list.pushBack(nodes[8]);
   list.pushBack(nodes[3]);

   list.insert(nodes[2]);

   dumpList();
   TEST_ASSERT(list.size_ == 4);
   TEST_ASSERT(int(list.front_ - list.nodes_) == 5);
   TEST_ASSERT(int(list.back_  - list.nodes_) == 3);
   TEST_ASSERT(int(list.lower_ - list.nodes_) == 2);
   TEST_ASSERT(int(list.upper_ - list.nodes_) == 8);
   TEST_ASSERT(list.isValid());

}


void ListTest::testInsert3()
{
   printMethod(TEST_FUNC);
   list.initialize(nodes, n);

   list.pushBack(nodes[3]);
   list.pushBack(nodes[5]);
   list.pushBack(nodes[8]);

   list.insert(nodes[2]);

   dumpList();
   TEST_ASSERT(list.size_ == 4);
   TEST_ASSERT(int(list.front_ - list.nodes_) == 2);
   TEST_ASSERT(int(list.back_  - list.nodes_) == 8);
   TEST_ASSERT(int(list.lower_ - list.nodes_) == 2);
   TEST_ASSERT(int(list.upper_ - list.nodes_) == 8);
   TEST_ASSERT(list.isValid());

}


void ListTest::testInsert4()
{
   printMethod(TEST_FUNC);
   list.initialize(nodes, n);

   list.pushBack(nodes[5]);
   list.pushBack(nodes[8]);
   list.pushBack(nodes[3]);

   list.insert(nodes[9]);

   dumpList();
   TEST_ASSERT(list.size_ == 4);
   TEST_ASSERT(int(list.front_ - list.nodes_) == 5);
   TEST_ASSERT(int(list.back_  - list.nodes_) == 3);
   TEST_ASSERT(int(list.lower_ - list.nodes_) == 3);
   TEST_ASSERT(int(list.upper_ - list.nodes_) == 9);
   TEST_ASSERT(list.isValid());

}


void ListTest::testInsert5()
{
   printMethod(TEST_FUNC);
   list.initialize(nodes, n);

   list.pushBack(nodes[3]);
   list.pushBack(nodes[5]);
   list.pushBack(nodes[7]);

   list.insert(nodes[9]);

   dumpList();
   TEST_ASSERT(list.size_ == 4);
   TEST_ASSERT(int(list.front_ - list.nodes_) == 3);
   TEST_ASSERT(int(list.back_  - list.nodes_) == 9);
   TEST_ASSERT(int(list.lower_ - list.nodes_) == 3);
   TEST_ASSERT(int(list.upper_ - list.nodes_) == 9);
   TEST_ASSERT(list.isValid());

}


/* 
* Print integer indices of nodes in a linked list.
*/
void ListTest::dumpList()
{
   if (verbose() > 1 && list.front_ != 0) {

      Node<Data>* node = list.front_;
      while (node->next() != 0) {
         std::cout << int(node - list.nodes_) << std::endl;
         node = node->next();
      }
      std::cout << int(node - list.nodes_) << std::endl;
   }
}

TEST_BEGIN(ListTest)
TEST_ADD(ListTest, testInitialize)
TEST_ADD(ListTest, testPushBack)
TEST_ADD(ListTest, testPushFront)
TEST_ADD(ListTest, testPopBack)
TEST_ADD(ListTest, testInsertNext1)
TEST_ADD(ListTest, testInsertNext2)
TEST_ADD(ListTest, testInsertPrev1)
TEST_ADD(ListTest, testInsertPrev2)
TEST_ADD(ListTest, testRemove1)
TEST_ADD(ListTest, testRemove2)
TEST_ADD(ListTest, testRemove3)
TEST_ADD(ListTest, testInsert1)
TEST_ADD(ListTest, testInsert2)
TEST_ADD(ListTest, testInsert3)
TEST_ADD(ListTest, testInsert4)
TEST_ADD(ListTest, testInsert5)
TEST_END(ListTest)


#endif

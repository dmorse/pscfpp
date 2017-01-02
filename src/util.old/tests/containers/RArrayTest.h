#ifndef RARRAY_TEST_H
#define RARRAY_TEST_H

#include <test/UnitTest.h>
#include <test/UnitTestRunner.h>

#include <util/containers/RArray.h>
#include <util/containers/DArray.h>
#include <util/containers/Array.h>

using namespace Util;

class RArrayTest : public UnitTest 
{

public:

   void setUp() {}

   void tearDown() {}
  
   void testAssociate();

   void testBaseClassReference();

};


void RArrayTest::testAssociate()
{
   printMethod(TEST_FUNC);
   DArray<int> v;
   RArray<int> u;

   v.allocate(3);
   u.associate(v);
   v[0] = 3;
   v[1] = 4;
   v[2] = 5;
   TEST_ASSERT(u[0] == 3 );
   TEST_ASSERT(u[1] == 4 );
   TEST_ASSERT(u[2] == 5 );
} 

void RArrayTest::testBaseClassReference()
{
   printMethod(TEST_FUNC);
   DArray<int> v;
   RArray<int> u;

   v.allocate(3);
   u.associate(v);
   v[0] = 3;
   v[1] = 4;
   v[2] = 5;

   Array<int>& w = u;
   TEST_ASSERT(w[0] == 3);
   TEST_ASSERT(w[1] == 4);
   TEST_ASSERT(w[2] == 5);
} 

TEST_BEGIN(RArrayTest)
TEST_ADD(RArrayTest, testAssociate)
TEST_ADD(RArrayTest, testBaseClassReference)
TEST_END(RArrayTest)

#endif

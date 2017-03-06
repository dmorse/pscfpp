#ifndef CYLN_DOMAIN_TEST_H
#define CYLN_DOMAIN_TEST_H

#include <test/UnitTest.h>
#include <test/UnitTestRunner.h>

#include <cyln/misc/Domain.h>

#include <util/containers/DArray.h>
#include <util/math/Constants.h>
#include <util/format/Dbl.h>

using namespace Util;
using namespace Pscf::Cyln;

class DomainTest : public UnitTest 
{
public:

   void setUp() 
   {  }

   void tearDown() {}

   void testConstructor();
   void testSetParameters();

};

void DomainTest::testConstructor()
{
   printMethod(TEST_FUNC);
   {
      Domain v;
      //TEST_ASSERT(v.capacity() == 0 );
      //TEST_ASSERT(!v.isAllocated() );
   }
} 

void DomainTest::testSetParameters() 
{
   printMethod(TEST_FUNC);
   printEndl();

   Domain v;
   v.setParameters(2.0, 3.0, 21, 31);
   TEST_ASSERT(eq(v.radius(), 2.0));
   TEST_ASSERT(eq(v.length(), 3.0));
   TEST_ASSERT(eq(v.nr(), 21));
   TEST_ASSERT(eq(v.nz(), 31));
   TEST_ASSERT(eq(v.dr(), 0.1));
   TEST_ASSERT(eq(v.dz(), 0.1));
   TEST_ASSERT(eq(v.volume(), 12.0*Constants::Pi));
   //std::cout << std::endl;
}


TEST_BEGIN(DomainTest)
TEST_ADD(DomainTest, testConstructor)
TEST_ADD(DomainTest, testSetParameters)
TEST_END(DomainTest)

#endif

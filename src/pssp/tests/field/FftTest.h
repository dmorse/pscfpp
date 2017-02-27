#ifndef PSSP_FFT_TEST_H
#define PSSP_FFT_TEST_H

#include <test/UnitTest.h>
#include <test/UnitTestRunner.h>

#include <pssp/field/FFT.h>
#include <pssp/field/RField.h>
#include <pssp/field/RFieldDFT.h>

using namespace Util;
using namespace Pssp;

class FftTest : public UnitTest 
{
public:

   void setUp() 
   {  }

   void tearDown() {}

   void testConstructor();

};


void FftTest::testConstructor()
{
   printMethod(TEST_FUNC);
   {
      FFT<1> v;
      //TEST_ASSERT(v.capacity() == 0 );
      //TEST_ASSERT(!v.isAllocated() );
   }
} 


TEST_BEGIN(FftTest)
TEST_ADD(FftTest, testConstructor)
TEST_END(FftTest)

#endif

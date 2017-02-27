#ifndef PSSP_FFT_TEST_H
#define PSSP_FFT_TEST_H

#include <test/UnitTest.h>
#include <test/UnitTestRunner.h>

#include <pssp/field/FFT.h>
#include <pssp/field/RField.h>
#include <pssp/field/RFieldDFT.h>

#include <util/math/Constants.h>
#include <util/format/Dbl.h>

using namespace Util;
using namespace Pssp;

class FftTest : public UnitTest 
{
public:

   void setUp() 
   {  }

   void tearDown() {}

   void testConstructor();
   void testTransform1D();

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

void FftTest::testTransform1D() {
   printMethod(TEST_FUNC);
   printEndl();

   RField<1> in;
   RFieldDFT<1> out;
   int n = 10;
   IntVec<1> d;
   d[0] = n;
   in.allocate(d);
   out.allocate(d);

   double x;
   double twoPi = 2.0*Constants::Pi;
   for (int i = 0; i < n; ++i) {
      x = twoPi*float(i)/float(n); 
      in[i] = cos(x);
      std::cout << Dbl(in[i]);
   }
   std::cout << std::endl;

   FFT<1> v;
   v.setup(in, out);
   
   #if 0 
   double factor = 1.0/double(n); 
   for (int i = 0; i < n/2 + 1; ++i) {
      out[i][0] *= factor;
      out[i][1] *= factor;
      //std::cout << out[i][0] << "  " << out[i][1] << std::endl;
   }

   for (int i = 0; i < n; ++i) {
      std::cout << Dbl(in[i]);
   }
   std::cout << std::endl;
   #endif
}


TEST_BEGIN(FftTest)
TEST_ADD(FftTest, testConstructor)
TEST_ADD(FftTest, testTransform1D)
TEST_END(FftTest)

#endif

#ifndef PSSP_FFT_TEST_H
#define PSSP_FFT_TEST_H

#include <test/UnitTest.h>
#include <test/UnitTestRunner.h>

#include <pssp/field/FFT.h>
#include <pssp/field/RField.h>
#include <pssp/field/RFieldDft.h>

#include <util/math/Constants.h>
#include <util/format/Dbl.h>

using namespace Util;
using namespace Pscf::Pssp;

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
   RFieldDft<1> out;
   int n = 10;
   IntVec<1> d;
   d[0] = n;
   in.allocate(d);
   out.allocate(d);

   // Initialize input data
   double x;
   double twoPi = 2.0*Constants::Pi;
   for (int i = 0; i < n; ++i) {
      x = twoPi*float(i)/float(n); 
      in[i] = cos(x);
      // std::cout << Dbl(in[i]);
   }
   // std::cout << std::endl;

   FFT<1> v;
   v.setup(in, out);
   v.forwardTransform(in, out);
   RField<1> inCopy;
   inCopy.allocate(d);
   v.inverseTransform(out, inCopy);

   for (int i = 0; i < n; ++i) {
      //std::cout << Dbl(inCopy[i]);
      TEST_ASSERT(eq(in[i], inCopy[i]));
   }
   //std::cout << std::endl;
}


TEST_BEGIN(FftTest)
TEST_ADD(FftTest, testConstructor)
TEST_ADD(FftTest, testTransform1D)
TEST_END(FftTest)

#endif

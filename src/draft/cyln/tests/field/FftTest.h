#ifndef CYLN_FFT_TEST_H
#define CYLN_FFT_TEST_H

#include <test/UnitTest.h>
#include <test/UnitTestRunner.h>

#include <cyln/field/FFT.h>

#include <util/containers/DArray.h>
#include <util/math/Constants.h>
#include <util/format/Dbl.h>

using namespace Util;
using namespace Pscf::Cyln;

class FftTest : public UnitTest 
{
public:

   void setUp() 
   {  }

   void tearDown() {}

   void testConstructor();
   void testTransform();

};

void FftTest::testConstructor()
{
   printMethod(TEST_FUNC);
   {
      FFT v;
      //TEST_ASSERT(v.capacity() == 0 );
      //TEST_ASSERT(!v.isAllocated() );
   }
} 

void FftTest::testTransform() 
{
   printMethod(TEST_FUNC);
   printEndl();

   DArray<double> in;
   DArray<fftw_complex> out;
   int n = 10;
   int m = n/2 + 1;
   in.allocate(n);
   out.allocate(m);

   // Initialize input data
   double x;
   double twoPi = 2.0*Constants::Pi;
   for (int i = 0; i < n; ++i) {
      x = twoPi*float(i)/float(n); 
      in[i] = cos(x);
      //std::cout << Dbl(in[i]);
   }
   //std::cout << std::endl;

   FFT v;
   v.setup(in, out);
   v.forwardTransform(in, out);
   DArray<double> inCopy;
   inCopy.allocate(n);
   v.inverseTransform(out, inCopy);

   for (int i = 0; i < n; ++i) {
      //std::cout << Dbl(inCopy[i]);
      TEST_ASSERT(eq(in[i], inCopy[i]));
   }
   //std::cout << std::endl;
}


TEST_BEGIN(FftTest)
TEST_ADD(FftTest, testConstructor)
TEST_ADD(FftTest, testTransform)
TEST_END(FftTest)

#endif

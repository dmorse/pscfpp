#ifndef PSPC_FFT_TEST_H
#define PSPC_FFT_TEST_H

#include <test/UnitTest.h>
#include <test/UnitTestRunner.h>

#include <pspc/field/FFT.h>
#include <pspc/field/RField.h>
#include <pspc/field/RFieldDft.h>
#include <pspc/field/RFieldComparison.h>

#include <util/math/Constants.h>
#include <util/format/Dbl.h>

using namespace Util;
using namespace Pscf::Pspc;

class FftTest : public UnitTest 
{
public:

   void setUp() {}
   void tearDown() {}

   void testConstructor();
   void testTransform1D();
   void testTransform2D();
   void testTransform3D();

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
   //printEndl();

   int n = 10;
   IntVec<1> d;
   d[0] = n;

   FFT<1> v;
   v.setup(d);

   RField<1> in;
   in.allocate(d);
   TEST_ASSERT(in.capacity() == n);

   // Initialize input data
   double x;
   double twoPi = 2.0*Constants::Pi;
   for (int i = 0; i < n; ++i) {
      x = twoPi*float(i)/float(n); 
      in[i] = cos(x);
   }

   // Transform in -> out
   RFieldDft<1> out;
   out.allocate(d);
   v.forwardTransform(in, out);

   #if 0
   // Save a copy of out (to check if input was overwritten)
   RFieldDft<1> outCopy(out);
   TEST_ASSERT(out.capacity() == outCopy.capacity());
   #endif

   // Inverse transform out -> inCopy
   RField<1> inCopy;
   inCopy.allocate(d);
   v.inverseTransform(out, inCopy);

   for (int i = 0; i < n; ++i) {
      TEST_ASSERT(eq(in[i], inCopy[i]));
   }

   #if 0
   // Test if input of inverse transform was overwritten
   for (int i = 0; i < out.capacity(); ++i) {
      TEST_ASSERT(eq(out[i][0], outCopy[i][0]));
      TEST_ASSERT(eq(out[i][1], outCopy[i][1]));
   }
   #endif

}

void FftTest::testTransform2D() 
{
   printMethod(TEST_FUNC);
   //printEndl();

   int n1 = 3;
   int n2 = 3;
   IntVec<2> d;
   d[0] = n1;
   d[1] = n2;

   FFT<2> v;
   v.setup(d);

   RField<2> in;
   in.allocate(d);
   int rank = 0;
   double x, y, cx, sy;
   double twoPi = 2.0*Constants::Pi;
   for (int i = 0; i < n1; i++) {
      x = twoPi*float(i)/float(n1); 
      cx = cos(x);
      for (int j = 0; j < n2; j++) {
         y = twoPi*float(j)/float(n2); 
         sy = sin(y);
         rank = j + (i * n2);
         in[rank] = 0.5 + 0.2*cx + 0.6*cx*cx - 0.1*sy + 0.3*cx*sy;
      }
   }

   // Forward transform in -> out
   RFieldDft<2> out;
   out.allocate(d);
   TEST_ASSERT(eq(in.capacity() / in.meshDimensions()[1],
                  out.capacity() / (out.meshDimensions()[1]/2 + 1)));
   v.forwardTransform(in, out);

   #if 1
   // Save a copy of out
   RFieldDft<2> outCopy(out);
   TEST_ASSERT(out.capacity() == outCopy.capacity());
   for (int i = 0; i < out.capacity(); ++i) {
      TEST_ASSERT(eq(out[i][0], outCopy[i][0]));
      TEST_ASSERT(eq(out[i][1], outCopy[i][1]));
   }
   #endif

   // Inverse transform out -> inCopy
   RField<2> inCopy;
   inCopy.allocate(d);
   v.inverseTransform(out, inCopy);

   // Check if in == inCopy
   for (int i = 0; i < n1; i++) {
      for (int j = 0; j < n2; j++) {
         rank = j + (i * n2);
         TEST_ASSERT(eq(in[rank], inCopy[rank]));
      }
   }

   #if 0
   // Check if out was not modified my inverseTransform
   for (int i = 0; i < out.capacity(); ++i) {
      TEST_ASSERT(eq(out[i][0], outCopy[i][0]));
      TEST_ASSERT(eq(out[i][1], outCopy[i][1]));
   }
   #endif

}

void FftTest::testTransform3D() {
   printMethod(TEST_FUNC);
   //printEndl();

   int n1 = 3;
   int n2 = 3;
   int n3 = 3;
   IntVec<3> d;
   d[0] = n1;
   d[1] = n2;
   d[2] = n3;

   FFT<3> v;
   v.setup(d);

   RField<3> in;
   RFieldDft<3> out;
   in.allocate(d);
   out.allocate(d);

   TEST_ASSERT(eq(in.capacity() / in.meshDimensions()[2],
                  out.capacity() / (out.meshDimensions()[2]/2 + 1)));

   int rank = 0;
   for (int i = 0; i < n1; i++) {
      for (int j = 0; j < n2; j++) {
         for (int k = 0; k < n3; k++){
            rank = k + ((j + (i * n2)) * n3);
            in[rank] = 1.0 + double(rank)/double(in.capacity());
         }
      }
   }

   v.forwardTransform(in, out);
   RField<3> inCopy;
   inCopy.allocate(d);
   v.inverseTransform(out, inCopy);

   for (int i = 0; i < n1; i++) {
      for (int j = 0; j < n2; j++) {
         for (int k = 0; k < n3; k++){
            rank = k + ((j + (i * n1)) * n3);
            TEST_ASSERT(eq(in[rank], inCopy[rank]));
         }
      }
   }

   RFieldComparison<3> comparison;
   comparison.compare(in, inCopy);
   //std::cout << std::endl;
   //std::cout << "maxDiff = " 
   //          << Dbl(comparison.maxDiff(), 20, 13)
   //          << std::endl;
   TEST_ASSERT(comparison.maxDiff() < 1.0E-12);

}

TEST_BEGIN(FftTest)
TEST_ADD(FftTest, testConstructor)
TEST_ADD(FftTest, testTransform1D)
TEST_ADD(FftTest, testTransform2D)
TEST_ADD(FftTest, testTransform3D)
TEST_END(FftTest)

#endif

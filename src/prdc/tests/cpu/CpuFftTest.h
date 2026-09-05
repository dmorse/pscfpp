#ifndef PRDC_CPU_FFT_TEST_H
#define PRDC_CPU_FFT_TEST_H

#include <test/UnitTest.h>
#include <test/UnitTestRunner.h>

#include <prdc/field/cpp/FFT.h>
#include <prdc/field/cpp/RField.h>
#include <prdc/field/cpp/RFieldDft.h>
#include <prdc/field/cpp/RFieldComparison.h>

#include <util/math/Constants.h>
#include <util/format/Dbl.h>

using namespace Util;
using namespace Pscf::Prdc;

class CpuFftTest : public UnitTest 
{

public:

   void setUp() {}
   void tearDown() {}

   void testConstructor();
   void testTransformReal1D();
   void testTransformReal2D();
   void testTransformReal3D();

   void testTransformComplex1D();
   void testTransformComplex2D();
   void testTransformComplex3D();

};

void CpuFftTest::testConstructor()
{
   printMethod(TEST_FUNC);
   {
      FFT<1,CPT> v;
   }
} 

void CpuFftTest::testTransformReal1D() 
{
   printMethod(TEST_FUNC);

   int n = 10;
   IntVec<1> d;
   d[0] = n;

   FFT<1,CPT> v;
   v.setup(d);

   RField<1,CPT> in;
   in.allocate(d);
   TEST_ASSERT(in.capacity() == n);

   // Initialize input data
   double x;
   double twoPi = 2.0*Constants::Pi;
   for (int i = 0; i < n; ++i) {
      x = twoPi*float(i)/float(n); 
      in[i] = cos(x);
   }

   // Save a copy of in (to ensure input to forwardTransform is preserved)
   RField<1,CPT> inCopy(in);

   // Transform in -> out
   RFieldDft<1,CPT> out;
   out.allocate(d);
   v.forwardTransform(in, out);

   // Save a copy of out (to ensure input to inverseTransformSafe is preserved)
   RFieldDft<1,CPT> outCopy(out);

   // Inverse transform out -> in2
   RField<1,CPT> in2;
   in2.allocate(d);
   v.inverseTransformSafe(out, in2);

   // Elementwise compare in, in2, and inCopy
   for (int i = 0; i < in.capacity(); i++) {
      TEST_ASSERT(eq(in[i], in2[i]));
      TEST_ASSERT(eq(in[i], inCopy[i]));
   }

   // Elementwise compare out and outCopy
   for (int i = 0; i < out.capacity(); i++) {
      TEST_ASSERT(eq(out[i][0], outCopy[i][0]));
      TEST_ASSERT(eq(out[i][1], outCopy[i][1]));
   }

}

void CpuFftTest::testTransformReal2D() 
{
   printMethod(TEST_FUNC);

   // Create mesh
   int n1 = 3;
   int n2 = 4;
   IntVec<2> d;
   d[0] = n1;
   d[1] = n2;

   // Instantiate and initialize FFT
   FFT<2,CPT> v;
   v.setup(d);

   // Initialize input data
   RField<2,CPT> in;
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

   // Save a copy of in (to ensure input to forwardTransform is preserved)
   RField<2,CPT> inCopy(in);

   // Forward transform in -> out
   RFieldDft<2,CPT> out;
   out.allocate(d);
   TEST_ASSERT(eq(in.capacity() / in.meshDimensions()[1],
                  out.capacity() / (out.meshDimensions()[1]/2 + 1)));
   v.forwardTransform(in, out);

   // Save a copy of out (to ensure input to inverseTransformSafe is preserved)
   RFieldDft<2,CPT> outCopy(out);

   // Inverse transform out -> in2
   RField<2,CPT> in2;
   in2.allocate(d);
   v.inverseTransformSafe(out, in2);

   // Elementwise compare in, in2, and inCopy
   for (int i = 0; i < in.capacity(); i++) {
      TEST_ASSERT(eq(in[i], in2[i]));
      TEST_ASSERT(eq(in[i], inCopy[i]));
   }

   // Elementwise compare out and outCopy
   for (int i = 0; i < out.capacity(); i++) {
      TEST_ASSERT(eq(out[i][0], outCopy[i][0]));
      TEST_ASSERT(eq(out[i][1], outCopy[i][1]));
   }

}

void CpuFftTest::testTransformReal3D() 
{
   printMethod(TEST_FUNC);

   // Create mesh
   int n1 = 3;
   int n2 = 3;
   int n3 = 3;
   IntVec<3> d;
   d[0] = n1;
   d[1] = n2;
   d[2] = n3;

   // Instantiate and initialize objects
   FFT<3,CPT> v;
   v.setup(d);

   RField<3,CPT> in;
   RFieldDft<3,CPT> out;
   in.allocate(d);
   out.allocate(d);

   TEST_ASSERT(eq(in.capacity() / in.meshDimensions()[2],
                  out.capacity() / (out.meshDimensions()[2]/2 + 1)));

   // Generate test data
   int rank = 0;
   for (int i = 0; i < n1; i++) {
      for (int j = 0; j < n2; j++) {
         for (int k = 0; k < n3; k++){
            rank = k + ((j + (i * n2)) * n3);
            in[rank] = 1.0 + double(rank)/double(in.capacity());
         }
      }
   }

   // Save a copy of in (to ensure input to forwardTransform is preserved)
   RField<3,CPT> inCopy(in);

   // Forward transform in -> out
   v.forwardTransform(in, out);

   // Save a copy of out (to ensure input to inverseTransformSafe is preserved)
   RFieldDft<3,CPT> outCopy(out);

   // Inverse transform out -> in2
   RField<3,CPT> in2;
   in2.allocate(d);
   v.inverseTransformSafe(out, in2);

   // Elementwise compare in, in2, and inCopy
   for (int i = 0; i < in.capacity(); i++) {
      TEST_ASSERT(eq(in[i], in2[i]));
      TEST_ASSERT(eq(in[i], inCopy[i]));
   }

   // Elementwise compare out and outCopy
   for (int i = 0; i < out.capacity(); i++) {
      TEST_ASSERT(eq(out[i][0], outCopy[i][0]));
      TEST_ASSERT(eq(out[i][1], outCopy[i][1]));
   }

}

void CpuFftTest::testTransformComplex1D() 
{
   printMethod(TEST_FUNC);

   int n = 10;
   IntVec<1> d;
   d[0] = n;

   FFT<1,CPT> v;
   v.setup(d);

   CField<1,CPT> in;
   in.allocate(d);
   TEST_ASSERT(in.capacity() == n);

   // Initialize input data
   double x, c, s;
   double twoPi = 2.0*Constants::Pi;
   for (int i = 0; i < n; ++i) {
      x = twoPi*float(i)/float(n); 
      c = cos(x);
      s = sin(x);
      in[i][0] = c + 0.5*c*c;
      in[i][1] = c + s + 0.5*s*s;
   }

   // Save a copy of in (to ensure input to forwardTransform is preserved)
   CField<1,CPT> inCopy(in);

   // Transform in -> out
   CField<1,CPT> out;
   out.allocate(d);
   v.forwardTransform(in, out);

   // Save a copy of out (to ensure input to inverseTransform is preserved)
   CField<1,CPT> outCopy(out);

   // Inverse transform out -> in2
   CField<1,CPT> in2;
   in2.allocate(d);
   v.inverseTransform(out, in2);

   // Elementwise compare in, in2, and inCopy
   for (int i = 0; i < in.capacity(); i++) {
      TEST_ASSERT(eq(in[i][0], in2[i][0]));
      TEST_ASSERT(eq(in[i][1], in2[i][1]));
      TEST_ASSERT(eq(in[i][0], inCopy[i][0]));
      TEST_ASSERT(eq(in[i][1], inCopy[i][1]));
   }

   // Elementwise compare out and outCopy
   for (int i = 0; i < out.capacity(); i++) {
      TEST_ASSERT(eq(out[i][0], outCopy[i][0]));
      TEST_ASSERT(eq(out[i][1], outCopy[i][1]));
   }

}

void CpuFftTest::testTransformComplex2D() 
{
   printMethod(TEST_FUNC);

   int n1 = 3;
   int n2 = 4;
   IntVec<2> d;
   d[0] = n1;
   d[1] = n2;

   FFT<2,CPT> v;
   v.setup(d);

   // Initialize test data
   CField<2,CPT> in;
   in.allocate(d);
   int rank = 0;
   double x, y, cx, sx, cy, sy;
   double twoPi = 2.0*Constants::Pi;
   for (int i = 0; i < n1; i++) {
      x = twoPi*float(i)/float(n1); 
      cx = cos(x);
      sx = sin(x);
      for (int j = 0; j < n2; j++) {
         y = twoPi*float(j)/float(n2); 
         cy = cos(y);
         sy = sin(y);
         rank = j + (i * n2);
         in[rank][0] = 0.5 + 0.2*cx + 0.6*cx*cx*sy - 0.1*sy + 0.3*cx*sy;
         in[rank][1] = -0.2 - 0.2*sy + 0.4*sy*cx*sy + 0.2*cx - 0.7*sx*cy;
      }
   }

   // Save a copy of in (to ensure input to forwardTransform is preserved)
   CField<2,CPT> inCopy(in);

   // Forward transform in -> out
   CField<2,CPT> out;
   out.allocate(d);
   v.forwardTransform(in, out);

   // Save a copy of out (to ensure input to inverseTransform is preserved)
   CField<2,CPT> outCopy(out);

   // Inverse transform out -> in2
   CField<2,CPT> in2;
   in2.allocate(d);

   v.inverseTransform(out, in2);

   // Elementwise compare in, in2, and inCopy
   for (int i = 0; i < in.capacity(); i++) {
      TEST_ASSERT(eq(in[i][0], in2[i][0]));
      TEST_ASSERT(eq(in[i][1], in2[i][1]));
      TEST_ASSERT(eq(in[i][0], inCopy[i][0]));
      TEST_ASSERT(eq(in[i][1], inCopy[i][1]));
   }

   // Elementwise compare out and outCopy
   for (int i = 0; i < out.capacity(); i++) {
      TEST_ASSERT(eq(out[i][0], outCopy[i][0]));
      TEST_ASSERT(eq(out[i][1], outCopy[i][1]));
   }

}

void CpuFftTest::testTransformComplex3D() 
{
   printMethod(TEST_FUNC);

   // Create mesh
   int n1 = 3;
   int n2 = 3;
   int n3 = 3;
   IntVec<3> d;
   d[0] = n1;
   d[1] = n2;
   d[2] = n3;

   // Instantiate and initialize objects
   FFT<3,CPT> v;
   v.setup(d);

   CField<3,CPT> in;
   CField<3,CPT> out;
   in.allocate(d);
   out.allocate(d);

   // Generate test data
   int rank = 0;
   double frac;
   for (int i = 0; i < n1; i++) {
      for (int j = 0; j < n2; j++) {
         for (int k = 0; k < n3; k++){
            rank = k + ((j + (i * n2)) * n3);
            frac = double(rank)/double(in.capacity());
            in[rank][0] =  0.3 + frac;
            in[rank][1] = -2.0 + frac*frac - 3.0*frac;
         }
      }
   }

   // Save a copy of in (to ensure input to forwardTransform is preserved)
   CField<3,CPT> inCopy(in);

   // Forward transform in -> out
   v.forwardTransform(in, out);

   // Save a copy of out (to ensure input to inverseTransform is preserved)
   CField<3,CPT> outCopy(out);

   // Inverse transform out -> in2
   CField<3,CPT> in2;
   in2.allocate(d);
   v.inverseTransform(out, in2);

   // Elementwise compare in, in2, and inCopy
   for (int i = 0; i < in.capacity(); i++) {
      TEST_ASSERT(eq(in[i][0], in2[i][0]));
      TEST_ASSERT(eq(in[i][1], in2[i][1]));
      TEST_ASSERT(eq(in[i][0], inCopy[i][0]));
      TEST_ASSERT(eq(in[i][1], inCopy[i][1]));
   }

   // Elementwise compare out and outCopy
   for (int i = 0; i < out.capacity(); i++) {
      TEST_ASSERT(eq(out[i][0], outCopy[i][0]));
      TEST_ASSERT(eq(out[i][1], outCopy[i][1]));
   }

}

TEST_BEGIN(CpuFftTest)
TEST_ADD(CpuFftTest, testConstructor)
TEST_ADD(CpuFftTest, testTransformReal1D)
TEST_ADD(CpuFftTest, testTransformReal2D)
TEST_ADD(CpuFftTest, testTransformReal3D)
TEST_ADD(CpuFftTest, testTransformComplex1D)
TEST_ADD(CpuFftTest, testTransformComplex2D)
TEST_ADD(CpuFftTest, testTransformComplex3D)
TEST_END(CpuFftTest)

#endif

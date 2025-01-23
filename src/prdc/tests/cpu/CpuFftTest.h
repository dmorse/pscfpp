#ifndef PRDC_CPU_FFT_TEST_H
#define PRDC_CPU_FFT_TEST_H

#include <test/UnitTest.h>
#include <test/UnitTestRunner.h>

#include <prdc/cpu/FFT.h>
#include <prdc/cpu/RField.h>
#include <prdc/cpu/RFieldDft.h>
#include <prdc/cpu/RFieldComparison.h>

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
      Cpu::FFT<1> v;
   }
} 

void CpuFftTest::testTransformReal1D() 
{
   printMethod(TEST_FUNC);

   int n = 10;
   IntVec<1> d;
   d[0] = n;

   Cpu::FFT<1> v;
   v.setup(d);

   Cpu::RField<1> in;
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
   Cpu::RField<1> inCopy(in);

   // Transform in -> out
   Cpu::RFieldDft<1> out;
   out.allocate(d);
   v.forwardTransform(in, out);

   // Save a copy of out (to ensure input to inverseTransformSafe is preserved)
   Cpu::RFieldDft<1> outCopy(out);

   // Inverse transform out -> in2
   Cpu::RField<1> in2;
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
   Cpu::FFT<2> v;
   v.setup(d);

   // Initialize input data
   Cpu::RField<2> in;
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
   Cpu::RField<2> inCopy(in);

   // Forward transform in -> out
   Cpu::RFieldDft<2> out;
   out.allocate(d);
   TEST_ASSERT(eq(in.capacity() / in.meshDimensions()[1],
                  out.capacity() / (out.meshDimensions()[1]/2 + 1)));
   v.forwardTransform(in, out);

   // Save a copy of out (to ensure input to inverseTransformSafe is preserved)
   Cpu::RFieldDft<2> outCopy(out);

   // Inverse transform out -> in2
   Cpu::RField<2> in2;
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
   Cpu::FFT<3> v;
   v.setup(d);

   Cpu::RField<3> in;
   Cpu::RFieldDft<3> out;
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
   Cpu::RField<3> inCopy(in);

   // Forward transform in -> out
   v.forwardTransform(in, out);

   // Save a copy of out (to ensure input to inverseTransformSafe is preserved)
   Cpu::RFieldDft<3> outCopy(out);

   // Inverse transform out -> in2
   Cpu::RField<3> in2;
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

   Cpu::FFT<1> v;
   v.setup(d);

   Cpu::CField<1> in;
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
   Cpu::CField<1> inCopy(in);

   // Transform in -> out
   Cpu::CField<1> out;
   out.allocate(d);
   v.forwardTransform(in, out);

   // Save a copy of out (to ensure input to inverseTransform is preserved)
   Cpu::CField<1> outCopy(out);

   // Inverse transform out -> in2
   Cpu::CField<1> in2;
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

   Cpu::FFT<2> v;
   v.setup(d);

   // Initialize test data
   Cpu::CField<2> in;
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
   Cpu::CField<2> inCopy(in);

   // Forward transform in -> out
   Cpu::CField<2> out;
   out.allocate(d);
   v.forwardTransform(in, out);

   // Save a copy of out (to ensure input to inverseTransform is preserved)
   Cpu::CField<2> outCopy(out);

   // Inverse transform out -> in2
   Cpu::CField<2> in2;
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
   Cpu::FFT<3> v;
   v.setup(d);

   Cpu::CField<3> in;
   Cpu::CField<3> out;
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
   Cpu::CField<3> inCopy(in);

   // Forward transform in -> out
   v.forwardTransform(in, out);

   // Save a copy of out (to ensure input to inverseTransform is preserved)
   Cpu::CField<3> outCopy(out);

   // Inverse transform out -> in2
   Cpu::CField<3> in2;
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

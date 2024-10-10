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
   //printEndl();

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

   // Transform in -> out
   Cpu::RFieldDft<1> out;
   out.allocate(d);
   v.forwardTransform(in, out);

   // Save a copy of out (to check if input was overwritten)
   Cpu::RFieldDft<1> outCopy(out);
   TEST_ASSERT(out.capacity() == outCopy.capacity());

   // Inverse transform out -> inCopy
   Cpu::RField<1> inCopy;
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

void CpuFftTest::testTransformReal2D() 
{
   printMethod(TEST_FUNC);
   //printEndl();

   int n1 = 3;
   int n2 = 4;
   IntVec<2> d;
   d[0] = n1;
   d[1] = n2;

   Cpu::FFT<2> v;
   v.setup(d);

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

   // Forward transform in -> out
   Cpu::RFieldDft<2> out;
   out.allocate(d);
   TEST_ASSERT(eq(in.capacity() / in.meshDimensions()[1],
                  out.capacity() / (out.meshDimensions()[1]/2 + 1)));
   v.forwardTransform(in, out);

   #if 1
   // Save a copy of out
   Cpu::RFieldDft<2> outCopy(out);
   TEST_ASSERT(out.capacity() == outCopy.capacity());
   for (int i = 0; i < out.capacity(); ++i) {
      TEST_ASSERT(eq(out[i][0], outCopy[i][0]));
      TEST_ASSERT(eq(out[i][1], outCopy[i][1]));
   }
   #endif

   // Inverse transform out -> inCopy
   Cpu::RField<2> inCopy;
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

void CpuFftTest::testTransformReal3D() 
{
   printMethod(TEST_FUNC);
   //printEndl();

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

   // Forward transform in -> out
   v.forwardTransform(in, out);

   // Inverse transform out -> inCopy
   Cpu::RField<3> inCopy;
   inCopy.allocate(d);
   v.inverseTransform(out, inCopy);

   // Elementwise compare in and inCopy
   for (int i = 0; i < n1; i++) {
      for (int j = 0; j < n2; j++) {
         for (int k = 0; k < n3; k++){
            rank = k + ((j + (i * n1)) * n3);
            TEST_ASSERT(eq(in[rank], inCopy[rank]));
         }
      }
   }

   // Use RFieldComparison to compare in and inCopy
   Cpu::RFieldComparison<3> comparison;
   comparison.compare(in, inCopy);
   //std::cout << std::endl;
   //std::cout << "maxDiff = " 
   //          << Dbl(comparison.maxDiff(), 20, 13)
   //          << std::endl;
   TEST_ASSERT(comparison.maxDiff() < 1.0E-12);

}

void CpuFftTest::testTransformComplex1D() 
{
   printMethod(TEST_FUNC);
   //printEndl();

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

   // Transform in -> out
   Cpu::CField<1> out;
   out.allocate(d);
   v.forwardTransform(in, out);

   // Save a copy of out (to check if input was overwritten)
   Cpu::CField<1> outCopy(out);
   TEST_ASSERT(out.capacity() == outCopy.capacity());

   // Inverse transform out -> inCopy
   Cpu::CField<1> inCopy;
   inCopy.allocate(d);
   v.inverseTransform(out, inCopy);

   for (int i = 0; i < n; ++i) {
      TEST_ASSERT(eq(in[i][0], inCopy[i][0]));
      TEST_ASSERT(eq(in[i][1], inCopy[i][1]));
   }

   // Test if input of inverse transform was overwritten
   for (int i = 0; i < out.capacity(); ++i) {
      TEST_ASSERT(eq(out[i][0], outCopy[i][0]));
      TEST_ASSERT(eq(out[i][1], outCopy[i][1]));
   }

}

void CpuFftTest::testTransformComplex2D() 
{
   printMethod(TEST_FUNC);
   //printEndl();

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

   // Forward transform in -> out
   Cpu::CField<2> out;
   out.allocate(d);
   v.forwardTransform(in, out);

   #if 0
   // Save a copy of out
   Cpu::CField<2> outCopy(out);
   //outCopy.allocate(d);
   TEST_ASSERT(out.capacity() == outCopy.capacity());
   for (int i = 0; i < out.capacity(); ++i) {
      TEST_ASSERT(eq(out[i][0], outCopy[i][0]));
      TEST_ASSERT(eq(out[i][1], outCopy[i][1]));
   }

   // Inverse transform out -> inCopy
   Cpu::CField<2> inCopy;
   inCopy.allocate(d);

   v.inverseTransform(out, inCopy);

   // Check if in == inCopy
   for (int i = 0; i < n1; i++) {
      for (int j = 0; j < n2; j++) {
         rank = j + (i * n2);
         TEST_ASSERT(eq(in[rank][0], inCopy[rank][0]));
         TEST_ASSERT(eq(in[rank][1], inCopy[rank][1]));
      }
   }

   // Check if out was not modified my inverseTransform
   for (int i = 0; i < out.capacity(); ++i) {
      TEST_ASSERT(eq(out[i][0], outCopy[i][0]));
      TEST_ASSERT(eq(out[i][1], outCopy[i][1]));
   }
   #endif

}

void CpuFftTest::testTransformComplex3D() 
{
   printMethod(TEST_FUNC);
   //printEndl();

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

   // Forward transform in -> out
   v.forwardTransform(in, out);

   // Inverse transform out -> inCopy
   Cpu::CField<3> inCopy;
   inCopy.allocate(d);
   v.inverseTransform(out, inCopy);

   // Elementwise compare in and inCopy
   for (int i = 0; i < n1; i++) {
      for (int j = 0; j < n2; j++) {
         for (int k = 0; k < n3; k++){
            rank = k + ((j + (i * n1)) * n3);
            TEST_ASSERT(eq(in[rank][0], inCopy[rank][0]));
            TEST_ASSERT(eq(in[rank][1], inCopy[rank][1]));
         }
      }
   }

   #if 0
   // Use CFieldComparison to compare in and inCopy
   Cpu::CFieldComparison<3> comparison;
   comparison.compare(in, inCopy);
   //std::cout << std::endl;
   //std::cout << "maxDiff = " 
   //          << Dbl(comparison.maxDiff(), 20, 13)
   //          << std::endl;
   TEST_ASSERT(comparison.maxDiff() < 1.0E-12);
   #endif

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

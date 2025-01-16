#ifndef PRDC_CUDA_FFT_TEST_H
#define PRDC_CUDA_FFT_TEST_H

#include <test/UnitTest.h>
#include <test/UnitTestRunner.h>

#include <prdc/cuda/FFT.h>
#include <prdc/cuda/RField.h>
#include <prdc/cuda/RFieldDft.h>

#include <pscf/cuda/HostDArray.h>
#include <pscf/math/IntVec.h>

#include <util/math/Constants.h>
#include <util/format/Dbl.h>

using namespace Util;
using namespace Pscf;
using namespace Pscf::Prdc;

class CudaFftTest : public UnitTest {

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

void CudaFftTest::testConstructor()
{
   printMethod(TEST_FUNC);
   {
      using namespace Pscf::Prdc::Cuda;
      Prdc::Cuda::FFT<1> v;
   }
} 

void CudaFftTest::testTransformReal1D() 
{
   printMethod(TEST_FUNC);

   // Data size
   int n = 100;
   int rSize = n;
   IntVec<1> d;
   d[0] = n;

   Prdc::Cuda::RField<1> d_rField;
   Prdc::Cuda::RFieldDft<1> d_kField;
   d_rField.allocate(d);
   d_kField.allocate(d);

   Prdc::Cuda::FFT<1> v;
   v.setup(d);

   TEST_ASSERT(d_rField.capacity() == n);

   // Initialize input data in a temporary array in host memory 
   HostDArray<cudaReal> in(rSize);

   // Initialize input data in a temporary array in host memory 
   double x;
   double twoPi = 2.0*Constants::Pi;
   for (int i = 0; i < n; ++i) {
      x = twoPi*float(i)/float(n); 
      in[i] = cos(x);
   }

   // Copy data to device
   d_rField = in;

   // Transform forward, r to k
   v.forwardTransform(d_rField, d_kField);

   // Inverse transform, k to r
   Prdc::Cuda::RField<1> d_rField_out;
   d_rField_out.allocate(d);
   v.inverseTransform(d_kField, d_rField_out);

   // Copy to host memory
   HostDArray<cudaReal> out(rSize);
   out = d_rField_out;

   // Test agreement after transform and round trip
   for (int i = 0; i < n; ++i) {
      TEST_ASSERT( std::abs(in[i] - out[i]) < 1E-10 );
   }

}

void CudaFftTest::testTransformReal2D() {
   printMethod(TEST_FUNC);
   
   int n1 = 3;
   int n2 = 3;
   int rSize= n1*n2;
   IntVec<2> d;
   d[0] = n1;
   d[1] = n2;

   Prdc::Cuda::RField<2> d_rField;
   Prdc::Cuda::RFieldDft<2> d_kField;
   d_rField.allocate(d);
   d_kField.allocate(d);
   TEST_ASSERT(d_rField.capacity() == rSize);

   Prdc::Cuda::FFT<2> v;
   v.setup(d);

   // Initialize input data in a temporary array in host memory 
   HostDArray<cudaReal> in(rSize);

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

   // Copy data from host to device
   d_rField = in;

   // Transform forward, r to k
   v.forwardTransform(d_rField, d_kField);

   // Inverse transform, k to r
   Prdc::Cuda::RField<2> d_rField_out;
   d_rField_out.allocate(d);
   v.inverseTransform(d_kField, d_rField_out);

   // Copy to host memory
   HostDArray<cudaReal> out(rSize);
   out = d_rField_out;

   // Test agreement
   for (int i = 0; i < n1*n2; ++i) {
      TEST_ASSERT( std::abs(in[i] - out[i]) < 1E-10 );
   }

}

void CudaFftTest::testTransformReal3D() {
   printMethod(TEST_FUNC);
   
   int n1 = 1;
   int n2 = 3;
   int n3 = 4;
   int rSize= n1*n2*n3;
   IntVec<3> d;
   d[0] = n1;
   d[1] = n2;
   d[2] = n3;

   Prdc::Cuda::RField<3> d_rField;
   Prdc::Cuda::RFieldDft<3> d_kField;
   d_rField.allocate(d);
   d_kField.allocate(d);

   Prdc::Cuda::FFT<3> v;
   v.setup(d);
   //v.setup(d_rField, d_kField);

   TEST_ASSERT(d_rField.capacity() == rSize);

   // Declare input array in host memory
   HostDArray<cudaReal> in(rSize);

   // Initialize input data in host memory 
   int rank = 0;
   for (int i = 0; i < n1; i++) {
      for (int j = 0; j < n2; j++) {
         for (int k = 0; k < n3; k++){
            rank = k + ((j + (i * n2)) * n3);
            in[rank] = 1.0 + double(rank)/double(n1*n2*n3);
         }
      }
   }

   // Copy data to device
   d_rField = in;

   // Transform forward, r to k
   v.forwardTransform(d_rField, d_kField);

   // Inverse transform, k to r
   Prdc::Cuda::RField<3> d_rField_out;
   d_rField_out.allocate(d);
   v.inverseTransform(d_kField, d_rField_out);
     
   // Copy to host memory
   HostDArray<cudaReal> out(rSize);
   out = d_rField_out;

   // Test round trip agreement
   for (int i = 0; i < n1*n2*n3; ++i) {
      TEST_ASSERT( std::abs(in[i] - out[i]) < 1E-10 );
   }

}

void CudaFftTest::testTransformComplex1D() 
{
   printMethod(TEST_FUNC);

   // Data size
   int n = 100;
   int rSize = n;
   IntVec<1> d;
   d[0] = n;

   Prdc::Cuda::CField<1> d_rField;
   Prdc::Cuda::CField<1> d_kField;
   d_rField.allocate(d);
   d_kField.allocate(d);

   Prdc::Cuda::FFT<1> v;
   v.setup(d);

   TEST_ASSERT(d_rField.capacity() == n);

   // Initialize input data in a temporary array in host memory 
   HostDArray<cudaComplex> in(rSize);

   // Initialize input data in a temporary array in host memory 
   double x, c, s;
   double twoPi = 2.0*Constants::Pi;
   for (int i = 0; i < n; ++i) {
      x = twoPi*float(i)/float(n); 
      c = cos(x);
      s = sin(x);
      in[i].x = c + 0.5*c*c;
      in[i].y = c + s + 0.5*c*c;
   }

   // Copy data to device
   d_rField = in;

   // Transform forward, r to k
   v.forwardTransform(d_rField, d_kField);

   // Inverse transform, k to r
   Prdc::Cuda::CField<1> d_rField_out;
   d_rField_out.allocate(d);
   v.inverseTransform(d_kField, d_rField_out);

   // Copy to host memory
   HostDArray<cudaComplex> out(rSize);
   out = d_rField_out;

   // Test agreement after transform and round trip
   for (int i = 0; i < n; ++i) {
      TEST_ASSERT( std::abs(in[i].x - out[i].x) < 1E-10 );
      TEST_ASSERT( std::abs(in[i].y - out[i].y) < 1E-10 );
   }

}

void CudaFftTest::testTransformComplex2D() {
   printMethod(TEST_FUNC);
   
   int n1 = 3;
   int n2 = 4;
   int rSize= n1*n2;
   IntVec<2> d;
   d[0] = n1;
   d[1] = n2;

   Prdc::Cuda::CField<2> d_rField;
   Prdc::Cuda::CField<2> d_kField;
   d_rField.allocate(d);
   d_kField.allocate(d);
   TEST_ASSERT(d_rField.capacity() == rSize);

   Prdc::Cuda::FFT<2> v;
   v.setup(d);

   // Initialize input data in a temporary array in host memory 
   HostDArray<cudaComplex> in(rSize);

   int rank = 0;
   double x, y, cx, cy, sx, sy;
   double twoPi = 2.0*Constants::Pi;
   for (int i = 0; i < n1; i++) {
      x = twoPi*float(i)/float(n1); 
      cx = cos(x);
      cy = sin(x);
      for (int j = 0; j < n2; j++) {
         y = twoPi*float(j)/float(n2); 
         cy = cos(y);
         sy = sin(y);
         rank = j + (i * n2);
         in[rank].x =  0.5 + 0.2*cx + 0.6*cx*cx*sy - 0.1*sy + 0.3*cx*sy;
         in[rank].y = -0.2 - 0.2*cy + 0.4*sy*cx*sy + 0.2*cx - 0.7*sx*cy;
      }
   }

   // Copy data from host to device
   d_rField = in;

   // Transform forward, r to k
   v.forwardTransform(d_rField, d_kField);

   // Inverse transform, k to r
   Prdc::Cuda::CField<2> d_rField_out;
   d_rField_out.allocate(d);
   v.inverseTransform(d_kField, d_rField_out);

   // Copy to host memory
   HostDArray<cudaComplex> out(rSize);
   out = d_rField_out;

   // Test agreement
   for (int i = 0; i < n1*n2; ++i) {
      TEST_ASSERT( std::abs(in[i].x - out[i].x) < 1E-10 );
      TEST_ASSERT( std::abs(in[i].y - out[i].y) < 1E-10 );
   }

}

void CudaFftTest::testTransformComplex3D() {
   printMethod(TEST_FUNC);
   
   int n1 = 1;
   int n2 = 3;
   int n3 = 4;
   int rSize= n1*n2*n3;
   IntVec<3> d;
   d[0] = n1;
   d[1] = n2;
   d[2] = n3;

   Prdc::Cuda::CField<3> d_rField;
   Prdc::Cuda::CField<3> d_kField;
   d_rField.allocate(d);
   d_kField.allocate(d);

   Prdc::Cuda::FFT<3> v;
   v.setup(d);

   TEST_ASSERT(d_rField.capacity() == rSize);

   // Declare input array in host memory
   HostDArray<cudaComplex> in(rSize);

   // Initialize input data in host memory 
   int rank = 0;
   double frac;
   for (int i = 0; i < n1; i++) {
      for (int j = 0; j < n2; j++) {
         for (int k = 0; k < n3; k++){
            rank = k + ((j + (i * n2)) * n3);
            frac = double(rank)/double(n1*n2*n3);
	    in[rank].x =  0.3 + frac;
	    in[rank].y = -2.0 + frac*frac - 3.0*frac;
         }
      }
   }

   // Copy data to device
   d_rField = in;

   // Transform forward, r to k
   v.forwardTransform(d_rField, d_kField);

   // Inverse transform, k to r
   Prdc::Cuda::CField<3> d_rField_out;
   d_rField_out.allocate(d);
   v.inverseTransform(d_kField, d_rField_out);

   // Copy to host memory
   HostDArray<cudaComplex> out(rSize);
   out = d_rField_out;

   // Test round trip agreement
   for (int i = 0; i < n1*n2*n3; ++i) {
      TEST_ASSERT( std::abs(in[i].x - out[i].x) < 1E-10 );
      TEST_ASSERT( std::abs(in[i].y - out[i].y) < 1E-10 );
   }

}
#if 0
#endif

TEST_BEGIN(CudaFftTest)
TEST_ADD(CudaFftTest, testConstructor)
TEST_ADD(CudaFftTest, testTransformReal1D)
TEST_ADD(CudaFftTest, testTransformReal2D)
TEST_ADD(CudaFftTest, testTransformReal3D)
TEST_ADD(CudaFftTest, testTransformComplex1D)
TEST_ADD(CudaFftTest, testTransformComplex2D)
TEST_ADD(CudaFftTest, testTransformComplex3D)
TEST_END(CudaFftTest)

#endif

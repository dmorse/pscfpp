#ifndef PRDC_CUDA_FFT_TEST_H
#define PRDC_CUDA_FFT_TEST_H

#include <test/UnitTest.h>
#include <test/UnitTestRunner.h>

//#include <pspg/field/RFieldComparison.h>

#include <prdc/cuda/FFT.h>
#include <prdc/cuda/RField.h>
#include <prdc/cuda/RFieldDft.h>

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
   void testTransform1D();
   void testTransform2D();
   void testTransform3D();

};

void CudaFftTest::testConstructor()
{
   printMethod(TEST_FUNC);
   {
      using namespace Pscf::Prdc::Cuda;
      Cuda::FFT<1> v;
   }
} 

void CudaFftTest::testTransform1D() 
{
   printMethod(TEST_FUNC);

   // Data size
   int n = 100;
   IntVec<1> d;
   d[0] = n;

   Cuda::RField<1> d_rField;
   Cuda::RFieldDft<1> d_kField;
   d_rField.allocate(d);
   d_kField.allocate(d);

   Cuda::FFT<1> v;
   v.setup(d_rField, d_kField);

   TEST_ASSERT(d_rField.capacity() == n);

   // Initialize input data in a temporary array in host memory 
   cudaReal* in = new cudaReal[n];
   double x;
   double twoPi = 2.0*Constants::Pi;
   for (int i = 0; i < n; ++i) {
      x = twoPi*float(i)/float(n); 
      in[i] = cos(x);
   }

   // Copy data to device
   cudaMemcpy(d_rField.cField(), in, 
            n*sizeof(cudaReal), cudaMemcpyHostToDevice);

   // Transform forward, r to k
   v.forwardTransform(d_rField, d_kField);

   // Inverse transform, k to r
   Cuda::RField<1> d_rField_out;
   d_rField_out.allocate(d);
   v.inverseTransform(d_kField, d_rField_out);

   // Copy to host memory
   cudaReal* out = new cudaReal[n];
   cudaMemcpy(out, d_rField_out.cField(), 
            n*sizeof(cudaReal), cudaMemcpyDeviceToHost);

   for (int i = 0; i < n; ++i) {
      TEST_ASSERT( std::abs(in[i] - out[i]) < 1E-10 );
   }

}

void CudaFftTest::testTransform2D() {
   printMethod(TEST_FUNC);
   
   int n1 = 3, n2 = 3;
   IntVec<2> d;
   d[0] = n1;
   d[1] = n2;

   Cuda::RField<2> d_rField;
   Cuda::RFieldDft<2> d_kField;
   d_rField.allocate(d);
   d_kField.allocate(d);

   Cuda::FFT<2> v;
   v.setup(d_rField, d_kField);

   TEST_ASSERT(d_rField.capacity() == n1*n2);

   // Initialize input data in a temporary array in host memory 
   cudaReal* in = new cudaReal[n1*n2];

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

   // Copy data to device
   cudaMemcpy(d_rField.cField(), in, 
            n1*n2*sizeof(cudaReal), cudaMemcpyHostToDevice);

   // Transform forward, r to k
   v.forwardTransform(d_rField, d_kField);

   // Inverse transform, k to r
   Cuda::RField<2> d_rField_out;
   d_rField_out.allocate(d);
   v.inverseTransform(d_kField, d_rField_out);

   // Copy to host memory
   cudaReal* out = new cudaReal[n1*n2];
   cudaMemcpy(out, d_rField_out.cField(), 
            n1*n2*sizeof(cudaReal), cudaMemcpyDeviceToHost);

   for (int i = 0; i < n1*n2; ++i) {
      TEST_ASSERT( std::abs(in[i] - out[i]) < 1E-10 );
   }

}

void CudaFftTest::testTransform3D() {
   printMethod(TEST_FUNC);
   
   int n1 = 3, n2 = 3, n3 = 3;
   IntVec<3> d;
   d[0] = n1;
   d[1] = n2;
   d[2] = n3;

   Cuda::RField<3> d_rField;
   Cuda::RFieldDft<3> d_kField;
   d_rField.allocate(d);
   d_kField.allocate(d);

   Cuda::FFT<3> v;
   v.setup(d_rField, d_kField);

   TEST_ASSERT(d_rField.capacity() == n1*n2*n3);

   // Initialize input data in a temporary array in host memory 
   cudaReal* in = new cudaReal[n1*n2*n3];

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
   cudaMemcpy(d_rField.cField(), in, 
            n1*n2*n3*sizeof(cudaReal), cudaMemcpyHostToDevice);

   // Transform forward, r to k
   v.forwardTransform(d_rField, d_kField);

   // Inverse transform, k to r
   Cuda::RField<3> d_rField_out;
   d_rField_out.allocate(d);
   v.inverseTransform(d_kField, d_rField_out);

   // Copy to host memory
   cudaReal* out = new cudaReal[n1*n2*n3];
   cudaMemcpy(out, d_rField_out.cField(), 
            n1*n2*n3*sizeof(cudaReal), cudaMemcpyDeviceToHost);

   for (int i = 0; i < n1*n2*n3; ++i) {
      TEST_ASSERT( std::abs(in[i] - out[i]) < 1E-10 );
   }

}

TEST_BEGIN(CudaFftTest)
TEST_ADD(CudaFftTest, testConstructor)
TEST_ADD(CudaFftTest, testTransform1D)
TEST_ADD(CudaFftTest, testTransform2D)
TEST_ADD(CudaFftTest, testTransform3D)
TEST_END(CudaFftTest)

#endif

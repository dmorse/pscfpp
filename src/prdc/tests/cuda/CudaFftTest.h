#ifndef PRDC_CUDA_FFT_TEST_H
#define PRDC_CUDA_FFT_TEST_H

#include <test/UnitTest.h>
#include <test/UnitTestRunner.h>

#include <prdc/cuda/FFT.h>
#include <prdc/cuda/FFTBatched.h>
#include <prdc/cuda/RField.h>
#include <prdc/cuda/RFieldDft.h>
#include <prdc/cuda/CField.h>

#include <pscf/cuda/HostDArray.h>
#include <pscf/math/IntVec.h>

#include <util/math/Constants.h>
#include <util/format/Dbl.h>

using namespace Util;
using namespace Pscf;
using namespace Pscf::Prdc;

class CudaFftTest : public UnitTest 
{


private:

   // Error tolerance for equality
   #ifdef SINGLE_PRECISION
   constexpr static cudaReal tolerance_ = 1E-5;
   #else
   #ifdef DOUBLE_PRECISION
   constexpr static cudaReal tolerance_ = 1E-10;
   #endif
   #endif

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

   void testBatchedTransformReal1D();
   void testBatchedTransformReal2D();
   void testBatchedTransformReal3D();

};

void CudaFftTest::testConstructor()
{
   printMethod(TEST_FUNC);
   Prdc::Cuda::FFT<3> v1;
   Prdc::Cuda::FFTBatched<3> v2;
} 

void CudaFftTest::testTransformReal1D() 
{
   printMethod(TEST_FUNC);

   // Create mesh
   int n = 100;
   int rSize = n;
   IntVec<1> d;
   d[0] = n;

   // Instantiate and allocate objects
   Prdc::Cuda::RField<1> rField(d);
   Prdc::Cuda::RFieldDft<1> kField(d);
   HostDArray<cudaReal> rField1_h(rSize), 
                        rField2_h(rSize), 
                        rField3_h(rSize);
   HostDArray<cudaComplex> kField1_h(kField.capacity()), 
                           kField2_h(kField.capacity());

   Prdc::Cuda::FFT<1> v;
   v.setup(d);

   // Initialize input data in host memory 
   double x;
   double twoPi = 2.0*Constants::Pi;
   for (int i = 0; i < n; ++i) {
      x = twoPi*float(i)/float(n); 
      rField1_h[i] = cos(x);
   }

   // Copy data to device
   rField = rField1_h;

   // Transform forward, r to k
   v.forwardTransform(rField, kField);

   // Save a copy of rField on host
   // (to ensure input to forwardTransform was preserved)
   rField2_h = rField;

   // Save a copy of kField on host
   // (to ensure input to inverseTransformSafe is preserved)
   kField1_h = kField;

   // Inverse transform, k to r
   v.inverseTransformSafe(kField, rField);

   // Copy to host memory
   rField3_h = rField;
   kField2_h = kField;

   // Elementwise compare rField arrays for equality
   for (int i = 0; i < rField.capacity(); i++) {
      TEST_ASSERT(fabs(rField1_h[i] - rField2_h[i]) < tolerance_);
      TEST_ASSERT(fabs(rField1_h[i] - rField3_h[i]) < tolerance_);
   }

   // Elementwise compare kField arrays for equality
   for (int i = 0; i < kField.capacity(); i++) {
      TEST_ASSERT(fabs(kField1_h[i].x - kField2_h[i].x) < tolerance_);
      TEST_ASSERT(fabs(kField1_h[i].y - kField2_h[i].y) < tolerance_);
   }

}

void CudaFftTest::testTransformReal2D() 
{
   printMethod(TEST_FUNC);
   
   // create mesh
   int n1 = 12;
   int n2 = 32;
   int rSize = n1*n2;
   IntVec<2> d;
   d[0] = n1;
   d[1] = n2;

   // Instantiate and allocate objects
   Prdc::Cuda::RField<2> rField(d);
   Prdc::Cuda::RFieldDft<2> kField(d);
   HostDArray<cudaReal> rField1_h(rSize), 
                        rField2_h(rSize), 
                        rField3_h(rSize);
   HostDArray<cudaComplex> kField1_h(kField.capacity()), 
                           kField2_h(kField.capacity());

   Prdc::Cuda::FFT<2> v;
   v.setup(d);

   // Initialize input data in host memory 
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
         rField1_h[rank] = 0.5 + 0.2*cx + 0.6*cx*cx - 0.1*sy + 0.3*cx*sy;
      }
   }

   // Copy data to device
   rField = rField1_h;

   // Transform forward, r to k
   v.forwardTransform(rField, kField);

   // Save a copy of rField on host
   // (to ensure input to forwardTransform was preserved)
   rField2_h = rField;

   // Save a copy of kField on host
   // (to ensure input to inverseTransformSafe is preserved)
   kField1_h = kField;

   // Inverse transform, k to r
   v.inverseTransformSafe(kField, rField);

   // Copy to host memory
   rField3_h = rField;
   kField2_h = kField;

   // Elementwise compare rField arrays for equality
   for (int i = 0; i < rField.capacity(); i++) {
      TEST_ASSERT(fabs(rField1_h[i] - rField2_h[i]) < tolerance_);
      TEST_ASSERT(fabs(rField1_h[i] - rField3_h[i]) < tolerance_);
   }

   // Elementwise compare kField arrays for equality
   for (int i = 0; i < kField.capacity(); i++) {
      TEST_ASSERT(fabs(kField1_h[i].x - kField2_h[i].x) < tolerance_);
      TEST_ASSERT(fabs(kField1_h[i].y - kField2_h[i].y) < tolerance_);
   }

}

void CudaFftTest::testTransformReal3D() 
{
   printMethod(TEST_FUNC);
   
   // create mesh
   int n1 = 12;
   int n2 = 32;
   int n3 = 16;
   int rSize= n1*n2*n3;
   IntVec<3> d;
   d[0] = n1;
   d[1] = n2;
   d[2] = n3;

   // Instantiate and allocate objects
   Prdc::Cuda::RField<3> rField(d);
   Prdc::Cuda::RFieldDft<3> kField(d);
   HostDArray<cudaReal> rField1_h(rSize), 
                        rField2_h(rSize), 
                        rField3_h(rSize);
   HostDArray<cudaComplex> kField1_h(kField.capacity()), 
                           kField2_h(kField.capacity());

   Prdc::Cuda::FFT<3> v;
   v.setup(d);

   // Initialize input data in host memory 
   int rank = 0;
   for (int i = 0; i < n1; i++) {
      for (int j = 0; j < n2; j++) {
         for (int k = 0; k < n3; k++){
            rank = k + ((j + (i * n2)) * n3);
            rField1_h[rank] = 1.0 + double(rank)/double(n1*n2*n3);
         }
      }
   }

   // Copy data to device
   rField = rField1_h;

   // Transform forward, r to k
   v.forwardTransform(rField, kField);

   // Save a copy of rField on host
   // (to ensure input to forwardTransform was preserved)
   rField2_h = rField;

   // Save a copy of kField on host
   // (to ensure input to inverseTransformSafe is preserved)
   kField1_h = kField;

   // Inverse transform, k to r
   v.inverseTransformSafe(kField, rField);

   // Copy to host memory
   rField3_h = rField;
   kField2_h = kField;

   // Elementwise compare rField arrays for equality
   for (int i = 0; i < rField.capacity(); i++) {
      TEST_ASSERT(fabs(rField1_h[i] - rField2_h[i]) < tolerance_);
      TEST_ASSERT(fabs(rField1_h[i] - rField3_h[i]) < tolerance_);
   }

   // Elementwise compare kField arrays for equality
   for (int i = 0; i < kField.capacity(); i++) {
      TEST_ASSERT(fabs(kField1_h[i].x - kField2_h[i].x) < tolerance_);
      TEST_ASSERT(fabs(kField1_h[i].y - kField2_h[i].y) < tolerance_);
   }

}

void CudaFftTest::testTransformComplex1D() 
{
   printMethod(TEST_FUNC);

   // create mesh
   int n = 100;
   int rSize = n;
   IntVec<1> d;
   d[0] = n;

   // Instantiate and allocate objects
   Prdc::Cuda::CField<1> rField(d), 
                         kField(d);
   HostDArray<cudaComplex> rField1_h(rSize), 
                           rField2_h(rSize), 
                           rField3_h(rSize), 
                           kField1_h(rSize), 
                           kField2_h(rSize);

   Prdc::Cuda::FFT<1> v;
   v.setup(d);

   // Initialize input data in a temporary array in host memory 
   double x, c, s;
   double twoPi = 2.0*Constants::Pi;
   for (int i = 0; i < n; ++i) {
      x = twoPi*float(i)/float(n); 
      c = cos(x);
      s = sin(x);
      rField1_h[i].x = c + 0.5*c*c;
      rField1_h[i].y = c + s + 0.5*c*c;
   }

   // Copy data to device
   rField = rField1_h;

   // Transform forward, r to k
   v.forwardTransform(rField, kField);

   // Save a copy of rField on host
   // (to ensure input to forwardTransform was preserved)
   rField2_h = rField;

   // Save a copy of kField on host
   // (to ensure input to inverseTransform is preserved)
   kField1_h = kField;

   // Inverse transform, k to r
   v.inverseTransform(kField, rField);

   // Copy to host memory
   rField3_h = rField;
   kField2_h = kField;

   // Elementwise compare rField arrays for equality
   for (int i = 0; i < rField.capacity(); i++) {
      TEST_ASSERT(fabs(rField1_h[i].x - rField2_h[i].x) < tolerance_);
      TEST_ASSERT(fabs(rField1_h[i].y - rField2_h[i].y) < tolerance_);
      TEST_ASSERT(fabs(rField1_h[i].x - rField3_h[i].x) < tolerance_);
      TEST_ASSERT(fabs(rField1_h[i].y - rField3_h[i].y) < tolerance_);
   }

   // Elementwise compare kField arrays for equality
   for (int i = 0; i < kField.capacity(); i++) {
      TEST_ASSERT(fabs(kField1_h[i].x - kField2_h[i].x) < tolerance_);
      TEST_ASSERT(fabs(kField1_h[i].y - kField2_h[i].y) < tolerance_);
   }

}

void CudaFftTest::testTransformComplex2D() {
   printMethod(TEST_FUNC);
   
   // create mesh
   int n1 = 12;
   int n2 = 32;
   int rSize= n1*n2;
   IntVec<2> d;
   d[0] = n1;
   d[1] = n2;

   // Instantiate and allocate objects
   Prdc::Cuda::CField<2> rField(d), 
                         kField(d);
   HostDArray<cudaComplex> rField1_h(rSize), 
                           rField2_h(rSize), 
                           rField3_h(rSize), 
                           kField1_h(rSize), 
                           kField2_h(rSize);

   Prdc::Cuda::FFT<2> v;
   v.setup(d);

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
         rField1_h[rank].x =  0.5 + 0.2*cx + 0.6*cx*cx*sy - 0.1*sy + 0.3*cx*sy;
         rField1_h[rank].y = -0.2 - 0.2*cy + 0.4*sy*cx*sy + 0.2*cx - 0.7*sx*cy;
      }
   }

   // Copy data to device
   rField = rField1_h;

   // Transform forward, r to k
   v.forwardTransform(rField, kField);

   // Save a copy of rField on host
   // (to ensure input to forwardTransform was preserved)
   rField2_h = rField;

   // Save a copy of kField on host
   // (to ensure input to inverseTransform is preserved)
   kField1_h = kField;

   // Inverse transform, k to r
   v.inverseTransform(kField, rField);

   // Copy to host memory
   rField3_h = rField;
   kField2_h = kField;

   // Elementwise compare rField arrays for equality
   for (int i = 0; i < rField.capacity(); i++) {
      TEST_ASSERT(fabs(rField1_h[i].x - rField2_h[i].x) < tolerance_);
      TEST_ASSERT(fabs(rField1_h[i].y - rField2_h[i].y) < tolerance_);
      TEST_ASSERT(fabs(rField1_h[i].x - rField3_h[i].x) < tolerance_);
      TEST_ASSERT(fabs(rField1_h[i].y - rField3_h[i].y) < tolerance_);
   }

   // Elementwise compare kField arrays for equality
   for (int i = 0; i < kField.capacity(); i++) {
      TEST_ASSERT(fabs(kField1_h[i].x - kField2_h[i].x) < tolerance_);
      TEST_ASSERT(fabs(kField1_h[i].y - kField2_h[i].y) < tolerance_);
   }

}

void CudaFftTest::testTransformComplex3D() {
   printMethod(TEST_FUNC);
   
   // create mesh
   int n1 = 12;
   int n2 = 32;
   int n3 = 16;
   int rSize= n1*n2*n3;
   IntVec<3> d;
   d[0] = n1;
   d[1] = n2;
   d[2] = n3;

   // Instantiate and allocate objects
   Prdc::Cuda::CField<3> rField(d), 
                         kField(d);
   HostDArray<cudaComplex> rField1_h(rSize), 
                           rField2_h(rSize), 
                           rField3_h(rSize), 
                           kField1_h(rSize), 
                           kField2_h(rSize);

   Prdc::Cuda::FFT<3> v;
   v.setup(d);

   // Initialize input data in host memory 
   int rank = 0;
   double frac;
   for (int i = 0; i < n1; i++) {
      for (int j = 0; j < n2; j++) {
         for (int k = 0; k < n3; k++) {
            rank = k + ((j + (i * n2)) * n3);
            frac = double(rank)/double(n1*n2*n3);
            rField1_h[rank].x =  0.3 + frac;
            rField1_h[rank].y = -2.0 + frac*frac - 3.0*frac;
         }
      }
   }

   // Copy data to device
   rField = rField1_h;

   // Transform forward, r to k
   v.forwardTransform(rField, kField);

   // Save a copy of rField on host
   // (to ensure input to forwardTransform was preserved)
   rField2_h = rField;

   // Save a copy of kField on host
   // (to ensure input to inverseTransform is preserved)
   kField1_h = kField;

   // Inverse transform, k to r
   v.inverseTransform(kField, rField);

   // Copy to host memory
   rField3_h = rField;
   kField2_h = kField;

   // Elementwise compare rField arrays for equality
   for (int i = 0; i < rField.capacity(); i++) {
      TEST_ASSERT(fabs(rField1_h[i].x - rField2_h[i].x) < tolerance_);
      TEST_ASSERT(fabs(rField1_h[i].y - rField2_h[i].y) < tolerance_);
      TEST_ASSERT(fabs(rField1_h[i].x - rField3_h[i].x) < tolerance_);
      TEST_ASSERT(fabs(rField1_h[i].y - rField3_h[i].y) < tolerance_);
   }

   // Elementwise compare kField arrays for equality
   for (int i = 0; i < kField.capacity(); i++) {
      TEST_ASSERT(fabs(kField1_h[i].x - kField2_h[i].x) < tolerance_);
      TEST_ASSERT(fabs(kField1_h[i].y - kField2_h[i].y) < tolerance_);
   }

}

void CudaFftTest::testBatchedTransformReal1D() 
{
   printMethod(TEST_FUNC);

   // Create mesh
   int n = 100;
   int rSize = n;
   int kSize = n / 2 + 1;
   IntVec<1> d;
   d[0] = n;

   // Number of FFTs in batch
   int batchSize = 3;

   // Instantiate and allocate objects
   DeviceArray<cudaReal> rField(batchSize * rSize);
   DeviceArray<cudaComplex> kField(batchSize * kSize);
   HostDArray<cudaReal> rField1_h(rField.capacity()), 
                        rField2_h(rField.capacity()), 
                        rField3_h(rField.capacity());
   HostDArray<cudaComplex> kField_h(kField.capacity());

   Prdc::Cuda::FFTBatched<1> v;
   v.setup(d, batchSize);

   // Initialize input data in host memory 
   double x;
   double twoPi = 2.0*Constants::Pi;
   for (int batch = 0; batch < batchSize; batch++) {
      int batchInd = batch * rSize;
      for (int i = 0; i < n; ++i) {
         x = twoPi*float(i)/float(n); 
         rField1_h[i + batchInd] = cos(x);
      }
   }

   // Copy data to device
   rField = rField1_h;

   // First, calculate FFT using FFT<1> object
   Prdc::Cuda::FFT<1> altFFT;
   altFFT.setup(d);
   Prdc::Cuda::RField<1> rFieldAlt;
   rFieldAlt.associate(rField, 0, d);
   Prdc::Cuda::RFieldDft<1> kFieldAlt(d);
   altFFT.forwardTransform(rFieldAlt, kFieldAlt);
   HostDArray<cudaComplex> kFieldAlt_h(kFieldAlt.capacity());
   kFieldAlt_h = kFieldAlt; 

   // Transform forward, r to k
   v.forwardTransform(rField, kField);

   // Save a copy of rField on host
   // (to ensure input to forwardTransform was preserved)
   rField2_h = rField;

   // Save a copy of kField on host before it is destroyed
   kField_h = kField;

   // Inverse transform, k to r
   v.inverseTransformUnsafe(kField, rField);

   // Copy rField to host memory
   rField3_h = rField;

   // Elementwise compare rField arrays for equality
   for (int i = 0; i < rField.capacity(); i++) {
      TEST_ASSERT(fabs(rField1_h[i] - rField2_h[i]) < tolerance_);
      TEST_ASSERT(fabs(rField1_h[i] - rField3_h[i]) < tolerance_);
   }

   // Elementwise compare kField arrays for equality
   for (int batch = 0; batch < batchSize; batch++) {
      int batchInd = batch * kSize;
      for (int i = 0; i < kSize; i++) {
         int k = batchInd + i;
         TEST_ASSERT(fabs(kField_h[k].x - kFieldAlt_h[i].x) < tolerance_);
         TEST_ASSERT(fabs(kField_h[k].y - kFieldAlt_h[i].y) < tolerance_);
      }
   }

}

void CudaFftTest::testBatchedTransformReal2D() 
{
   printMethod(TEST_FUNC);

   // create mesh
   int n1 = 12;
   int n2 = 32;
   int rSize = n1 * n2;
   int kSize = n1 * (n2/2 + 1);
   IntVec<2> d;
   d[0] = n1;
   d[1] = n2;

   // Number of FFTs in batch
   int batchSize = 3;

   // Instantiate and allocate objects
   DeviceArray<cudaReal> rField(batchSize * rSize);
   DeviceArray<cudaComplex> kField(batchSize * kSize);
   HostDArray<cudaReal> rField1_h(rField.capacity()), 
                        rField2_h(rField.capacity()), 
                        rField3_h(rField.capacity());
   HostDArray<cudaComplex> kField_h(kField.capacity());

   Prdc::Cuda::FFTBatched<2> v;
   v.setup(d, batchSize);

   // Initialize input data in host memory 
   int rank = 0;
   double x, y, cx, sy;
   double twoPi = 2.0*Constants::Pi;
   for (int batch = 0; batch < batchSize; batch++) {
      int batchInd = batch * rSize;
      for (int i = 0; i < n1; i++) {
         x = twoPi*float(i)/float(n1); 
         cx = cos(x);
         for (int j = 0; j < n2; j++) {
            y = twoPi*float(j)/float(n2); 
            sy = sin(y);
            rank = j + (i * n2) + batchInd;
            rField1_h[rank] = 0.5 + 0.2*cx + 0.6*cx*cx - 0.1*sy + 0.3*cx*sy;
         }
      }
   }

   // Copy data to device
   rField = rField1_h;

   // First, calculate FFT using FFT<2> object
   Prdc::Cuda::FFT<2> altFFT;
   altFFT.setup(d);
   Prdc::Cuda::RField<2> rFieldAlt;
   rFieldAlt.associate(rField, 0, d);
   Prdc::Cuda::RFieldDft<2> kFieldAlt(d);
   altFFT.forwardTransform(rFieldAlt, kFieldAlt);
   HostDArray<cudaComplex> kFieldAlt_h(kFieldAlt.capacity());
   kFieldAlt_h = kFieldAlt; 

   // Transform forward, r to k
   v.forwardTransform(rField, kField);

   // Save a copy of rField on host
   // (to ensure input to forwardTransform was preserved)
   rField2_h = rField;

   // Save a copy of kField on host before it is destroyed
   kField_h = kField;

   // Inverse transform, k to r
   v.inverseTransformUnsafe(kField, rField);

   // Copy rField to host memory
   rField3_h = rField;

   // Elementwise compare rField arrays for equality
   for (int i = 0; i < rField.capacity(); i++) {
      TEST_ASSERT(fabs(rField1_h[i] - rField2_h[i]) < tolerance_);
      TEST_ASSERT(fabs(rField1_h[i] - rField3_h[i]) < tolerance_);
   }

   // Elementwise compare kField arrays for equality
   for (int batch = 0; batch < batchSize; batch++) {
      int batchInd = batch * kSize;
      for (int i = 0; i < kSize; i++) {
         int k = batchInd + i;
         TEST_ASSERT(fabs(kField_h[k].x - kFieldAlt_h[i].x) < tolerance_);
         TEST_ASSERT(fabs(kField_h[k].y - kFieldAlt_h[i].y) < tolerance_);
      }
   }

}

void CudaFftTest::testBatchedTransformReal3D() 
{
   printMethod(TEST_FUNC);

   // create mesh
   int n1 = 12;
   int n2 = 32;
   int n3 = 16;
   int rSize= n1 * n2 * n3;
   int kSize= n1 * n2 * (n3/2 + 1);
   IntVec<3> d;
   d[0] = n1;
   d[1] = n2;
   d[2] = n3;

   // Number of FFTs in batch
   int batchSize = 3;

   // Instantiate and allocate objects
   DeviceArray<cudaReal> rField(batchSize * rSize);
   DeviceArray<cudaComplex> kField(batchSize * kSize);
   HostDArray<cudaReal> rField1_h(rField.capacity()), 
                        rField2_h(rField.capacity()), 
                        rField3_h(rField.capacity());
   HostDArray<cudaComplex> kField_h(kField.capacity());

   Prdc::Cuda::FFTBatched<3> v;
   v.setup(d, batchSize);

   // Initialize input data in host memory 
   for (int batch = 0; batch < batchSize; batch++) {
      int batchInd = batch * rSize;
      // Initialize input data in host memory 
      int rank;
      for (int i = 0; i < n1; i++) {
         for (int j = 0; j < n2; j++) {
            for (int k = 0; k < n3; k++){
               rank = k + ((j + (i * n2)) * n3);
               rField1_h[rank+batchInd] = 1.0 + double(rank)/double(n1*n2*n3);
            }
         }
      }
   }

   // Copy data to device
   rField = rField1_h;

   // First, calculate FFT using FFT<3> object
   Prdc::Cuda::FFT<3> altFFT;
   altFFT.setup(d);
   Prdc::Cuda::RField<3> rFieldAlt;
   rFieldAlt.associate(rField, 0, d);
   Prdc::Cuda::RFieldDft<3> kFieldAlt(d);
   altFFT.forwardTransform(rFieldAlt, kFieldAlt);
   HostDArray<cudaComplex> kFieldAlt_h(kFieldAlt.capacity());
   kFieldAlt_h = kFieldAlt; 

   // Transform forward, r to k
   v.forwardTransform(rField, kField);

   // Save a copy of rField on host
   // (to ensure input to forwardTransform was preserved)
   rField2_h = rField;

   // Save a copy of kField on host before it is destroyed
   kField_h = kField;

   // Inverse transform, k to r
   v.inverseTransformUnsafe(kField, rField);

   // Copy rField to host memory
   rField3_h = rField;

   // Elementwise compare rField arrays for equality
   for (int i = 0; i < rField.capacity(); i++) {
      TEST_ASSERT(fabs(rField1_h[i] - rField2_h[i]) < tolerance_);
      TEST_ASSERT(fabs(rField1_h[i] - rField3_h[i]) < tolerance_);
   }

   // Elementwise compare kField arrays for equality
   for (int batch = 0; batch < batchSize; batch++) {
      int batchInd = batch * kSize;
      for (int i = 0; i < kSize; i++) {
         int k = batchInd + i;
         if (fabs(kField_h[k].x - kFieldAlt_h[i].x) > tolerance_) {
            Log::file() << std::endl << kField_h[k].x << " " << kFieldAlt_h[i].x << " " << k << std::endl;
         }
         TEST_ASSERT(fabs(kField_h[k].x - kFieldAlt_h[i].x) < tolerance_);
         TEST_ASSERT(fabs(kField_h[k].y - kFieldAlt_h[i].y) < tolerance_);
      }
   }

}



TEST_BEGIN(CudaFftTest)
TEST_ADD(CudaFftTest, testConstructor)
TEST_ADD(CudaFftTest, testTransformReal1D)
TEST_ADD(CudaFftTest, testTransformReal2D)
TEST_ADD(CudaFftTest, testTransformReal3D)
TEST_ADD(CudaFftTest, testTransformComplex1D)
TEST_ADD(CudaFftTest, testTransformComplex2D)
TEST_ADD(CudaFftTest, testTransformComplex3D)
TEST_ADD(CudaFftTest, testBatchedTransformReal1D)
TEST_ADD(CudaFftTest, testBatchedTransformReal2D)
TEST_ADD(CudaFftTest, testBatchedTransformReal3D)
TEST_END(CudaFftTest)

#endif

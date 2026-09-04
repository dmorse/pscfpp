#ifndef PSCF_CUDA_REDUCE_TEST_H
#define PSCF_CUDA_REDUCE_TEST_H

#include <test/UnitTest.h>
#include <test/UnitTestRunner.h>

#include <pscf/backend/cuda/cudaTypes.h>
#include <pscf/backend/cuda/complex.h>
#include <pscf/backend/cuda/Reduce.h>
#include <pscf/backend/cuda/CudaVecRandom.h>
#include <pscf/backend/cuda/DeviceArray.h>
#include <pscf/backend/cuda/HostDArray.h>
#include <pscf/math/IntVec.h>

#include <util/format/Dbl.h>
#include <util/misc/Timer.h>

#include <cstdlib>
#include <cmath>

using namespace Util;
using namespace Pscf;

class CudaReduceTest : public UnitTest
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

   // Random number generator on the GPU
   CudaVecRandom rand_;

public:

   /// Setup before each test function.
   void setUp()
   {
      setVerbose(0);
      rand_.setSeed(0);
   }

   /// Clean up after each test function.
   void tearDown()
   {  
      // Reduce::freeWorkSpace(); 
   }

   /// Test sum for an array of real elements.
   void testSum()
   {
      printMethod(TEST_FUNC);

      IntVec<3> nVals;
      nVals[0] = 50022;     // small array
      nVals[1] = 1934896;   // medium array
      nVals[2] = 109857634; // large array

      for (int j = 0; j < 3; j++) {

         int n = nVals[j];

         if (verbose() > 0) {
            Log::file() << std::endl << "n = " << n << std::endl;
         }

         // Generate random test data,
         // normally distributed about 0.5 with stdev = 2
         DeviceArray<cudaReal> num(n);
         rand_.normal(num, (cudaReal)2.0, (cudaReal)0.5);

         // Copy test data to host
         HostDArray<cudaReal> num_h(n);
         num_h = num;

         // Determine highest power of 2 less than n
         int nReduced = (int)(pow(2.0,floor(log2(n))) + 0.5);
         // note: 0.5 added to make sure it casts to the correct int value

         // Find sum on host using a binary tree
         // (numerical round-off error should match that from GPU sum)
         Timer timerCPU;
         if (verbose() > 0) {
            timerCPU.start();
         }

         for (int i = 0; i < nReduced; i++) {
            if (i + nReduced < n) {
               num_h[i] += num_h[nReduced+i];
            }
         }
         nReduced /= 2;
         for ( ; nReduced >= 1; nReduced /= 2) {
            for (int i = 0; i < nReduced; i++) {
               num_h[i] += num_h[nReduced+i];
            }
         }
         cudaReal sumCPU = num_h[0];
         if (verbose() > 0) {
            timerCPU.stop();
            Log::file() << "CPU wall time: " << Dbl(timerCPU.time())
                        << std::endl;
         }

         // Call kernel wrapper to calculate sum on GPU
         Timer timerGPU;
         if (verbose() > 0) {
            timerGPU.start();
         }
         cudaReal sumGPU = Reduce::sum(num);

         // Check answer
         if (verbose() > 0) {
            timerGPU.stop();
            Log::file()
                 << "GPU wall time: " << Dbl(timerGPU.time()) << "\n"
                 << "Sum on CPU:    " << Dbl(sumCPU) << "\n"
                 << "Sum on GPU:    " << Dbl(sumGPU) << "\n"
                 << "Difference:    " << fabs(sumCPU - sumGPU) << "\n"
                 << std::endl;
         }

         // Check that error is at least 5 (10) orders of magnitude smaller
         // than the value of the sum for single (double) precision data
         TEST_ASSERT((fabs(sumCPU - sumGPU) / sumCPU) < tolerance_);
      }
   }

   /// Test sum for an array of real elements.
   void testSumSlice()
   {
      printMethod(TEST_FUNC);

      IntVec<3> nVals;
      nVals[0] = 538;    // small array
      nVals[1] = 2796;   // medium array

      for (int j = 0; j < 2; j++) {

         int n = nVals[j];
	 int begin = 7;
	 int end = n - 3;

         if (verbose() > 0) {
            Log::file() << std::endl << "n = " << n << std::endl;
         }

         // Generate random test data,
         // normally distributed about 0.5 with stdev = 2
         DeviceArray<cudaReal> num(n);
         rand_.normal(num, (cudaReal)2.0, (cudaReal)0.5);

         // Copy test data to host
         HostDArray<cudaReal> num_h(n);
         num_h = num;

         // Determine highest power of 2 less than n
         int nReduced = (int)(pow(2.0,floor(log2(n))) + 0.5);
         // note: 0.5 added to make sure it casts to the correct int value

         // Find sum on host using a binary tree
         // (numerical round-off error should match that from GPU sum)
         Timer timerCPU;
         if (verbose() > 0) {
            timerCPU.start();
         }

	 // Evaluate sum on CPU (simple loop)
	 cudaReal sumCPU = 0.0;
	 for (int i = begin; i < end; ++i) {
            sumCPU += num_h[i];
         }		  

         if (verbose() > 0) {
            timerCPU.stop();
            Log::file() << "CPU wall time: " << Dbl(timerCPU.time())
                        << std::endl;
         }

         // Call kernel wrapper to calculate sum on GPU
         Timer timerGPU;
         if (verbose() > 0) {
            timerGPU.start();
         }
         cudaReal sumGPU = Reduce::sum(num, begin, end);

         // Check answer
         if (verbose() > 0) {
            timerGPU.stop();
            Log::file()
                 << "GPU wall time: " << Dbl(timerGPU.time()) << "\n"
                 << "Sum on CPU:    " << Dbl(sumCPU) << "\n"
                 << "Sum on GPU:    " << Dbl(sumGPU) << "\n"
                 << "Difference:    " << fabs(sumCPU - sumGPU) << "\n"
                 << std::endl;
         }

         // Check that error is at least 5 (10) orders of magnitude smaller
         // than the value of the sum for single (double) precision data
         TEST_ASSERT((fabs(sumCPU - sumGPU) / sumCPU) < tolerance_);
      }
   }

   /// Test sum of an array of cudaComplex elements.
   void testSumComplex()
   {
      printMethod(TEST_FUNC);

      IntVec<3> nVals;
      nVals[0] = 109857634; // large array
      nVals[1] = 1934896;   // medium array
      nVals[2] = 50022;     // small array

      for (int j = 0; j < 3; j++) {

         int n = nVals[j];

         if (verbose() > 0) {
            Log::file() << std::endl << "n = " << n << std::endl;
         }

         // Generate random test data,
         // normally distributed about 0.5 with stdev = 2
         DeviceArray<cudaReal> num_dr(2*n);
         HostDArray<cudaReal>  num_hr(2*n);
         rand_.normal(num_dr, (cudaReal)2.0, (cudaReal)0.5);
         num_hr = num_dr;

          // Copy data to cudaComplex arrays
         HostDArray<cudaComplex>  num_h(n);
         DeviceArray<cudaComplex> num_d(n);
         //cudaComplex sum0 = makeComplex(0.0, 0.0);
         for (int i = 0; i < n; ++i) {
            num_h[i].x = num_hr[2*i];
            num_h[i].y = num_hr[2*i + 1];
            //sum0.x += num_h[i].x;
            //sum0.y += num_h[i].y;
         }
         num_d = num_h;

         // Determine highest power of 2 less than n
         int nReduced = (int)(pow(2.0,floor(log2(n))) + 0.5);
         // note: 0.5 added to make sure it casts to the correct int value

         // Find sum on host using a binary tree
         // (numerical round-off error should match that from GPU sum)
         Timer timerCPU;
         if (verbose() > 0) {
            timerCPU.start();
         }

         for (int i = 0; i < nReduced; i++) {
            if (i + nReduced < n) {
               addEq(num_h[i], num_h[nReduced+i]);
            }
         }
         nReduced /= 2;
         for ( ; nReduced >= 1; nReduced /= 2) {
            for (int i = 0; i < nReduced; i++) {
               //num_h[i] += num_h[nReduced+i];
               addEq(num_h[i], num_h[nReduced+i]);
            }
         }
         std::complex<cudaReal> sumCPU;
	 assign(sumCPU, num_h[0]);
         //sumCPU = std::complex<cudaReal>(num_h[0].x, num_h[0].y);
         if (verbose() > 0) {
            timerCPU.stop();
            Log::file() << "CPU wall time: " << Dbl(timerCPU.time())
                        << std::endl;
            #if 0
            Log::file() << "sum0  : "
                        << Dbl(sum0.x) << "  " << Dbl(sum0.y) << "\n"
            Log::file() << "sumCPU: "
                        << Dbl(sumCPU.real()) << "  " << Dbl(sumCPU.imag())
                        << std::endl;
            #endif
         }

         // Call kernel wrapper to calculate sum on GPU
         Timer timerGPU;
         if (verbose() > 0) {
            timerGPU.start();
         }
	 std::complex<cudaReal> sumGPU = Reduce::sum(num_d);
	 std::complex<cudaReal> diff = sumGPU - sumCPU;

         // Check answer
         if (verbose() > 0) {
            timerGPU.stop();
            Log::file() 
               << "GPU wall time: " << Dbl(timerGPU.time()) << "\n"
               << "Sum on CPU:    "
               << Dbl(sumCPU.real()) << "  " << Dbl(sumCPU.imag()) << "\n"
               << "Sum on GPU:    "
               << Dbl(sumGPU.real()) << "  " << Dbl(sumGPU.imag()) << "\n"
               << "Difference:    " << std::abs(diff) << std::endl;
         }

         // Check that error is at least 5 (10) orders of magnitude smaller
         // than the value of the sum for single (double) precision data
         cudaReal relDiff = std::abs(diff) / std::abs(sumCPU);
         TEST_ASSERT( relDiff < tolerance_);
      }
   }

   /// Test sum of an array of cudaComplex elements.
   void testSumComplexSlice()
   {
      printMethod(TEST_FUNC);

      IntVec<3> nVals;
      nVals[0] = 574; // large array
      nVals[1] = 2378;   // medium array

      for (int j = 0; j < 2; j++) {

         int n = nVals[j];
	 int begin = 13;
	 int end = n - 1;

         if (verbose() > 0) {
            Log::file() << std::endl << "n = " << n << std::endl;
         }

         // Generate random test data,
         // normally distributed about 0.5 with stdev = 2
         DeviceArray<cudaReal> num_dr(2*n);
         HostDArray<cudaReal>  num_hr(2*n);
         rand_.normal(num_dr, (cudaReal)2.0, (cudaReal)0.5);
         num_hr = num_dr;

          // Copy data to cudaComplex arrays
         HostDArray<cudaComplex>  num_h(n);
         DeviceArray<cudaComplex> num_d(n);
         //cudaComplex sum0 = makeComplex(0.0, 0.0);
         for (int i = 0; i < n; ++i) {
            num_h[i].x = num_hr[2*i];
            num_h[i].y = num_hr[2*i + 1];
         }
         num_d = num_h;

         Timer timerCPU;
         if (verbose() > 0) {
            timerCPU.start();
         }

	 // Evaluate sum on CPU
	 cudaComplex sum;
	 assign(sum, 0.0);
	 for (int i = begin; i < end; ++i) {
            addEq(sum, num_h[i]);
         }
         std::complex<cudaReal> sumCPU;
	 assign(sumCPU, sum);

         if (verbose() > 0) {
            timerCPU.stop();
            Log::file() << "CPU wall time: " << Dbl(timerCPU.time())
                        << std::endl;
         }

         // Call kernel wrapper to calculate sum on GPU
         Timer timerGPU;
         if (verbose() > 0) {
            timerGPU.start();
         }

	 // Evaluate sum on GPU
	 std::complex<cudaReal> sumGPU = Reduce::sum(num_d, begin, end);
	 std::complex<cudaReal> diff = sumGPU - sumCPU;

         // Check answer
         if (verbose() > 0) {
            timerGPU.stop();
            Log::file() 
               << "GPU wall time: " << Dbl(timerGPU.time()) << "\n"
               << "Sum on CPU:    "
               << Dbl(sumCPU.real()) << "  " << Dbl(sumCPU.imag()) << "\n"
               << "Sum on GPU:    "
               << Dbl(sumGPU.real()) << "  " << Dbl(sumGPU.imag()) << "\n"
               << "Difference:    " << std::abs(diff) << std::endl;
         }

         // Check that error is at least 5 (10) orders of magnitude smaller
         // than the value of the sum for single (double) precision data
         cudaReal relDiff = std::abs(diff) / std::abs(sumCPU);
         TEST_ASSERT( relDiff < tolerance_);
      }
   }

   /// Test sum of an array of real elements.
   void testSumSq()
   {
      printMethod(TEST_FUNC);

      IntVec<3> nVals;
      nVals[0] = 50022;     // small array
      nVals[1] = 1934896;   // medium array
      nVals[2] = 109857634; // large array

      for (int j = 0; j < 3; j++) {

         int n = nVals[j];

         if (verbose() > 0) {
            Log::file() << std::endl << "n = " << n << std::endl;
         }

         // Generate random test data on host and device,
         // normally distributed about 0.001 with stdev = 1.0
         HostDArray<cudaReal> num_h(n);
         DeviceArray<cudaReal> num_d(n);
         rand_.normal(num_d, (cudaReal)1.0, (cudaReal)0.01);
         num_h = num_d;

         // Convert host array num_h to array of squared values
	 cudaReal val, valSq, sum0;
	 sum0 = 0.0;
         for (int i = 0; i < n; ++i) {
	    val = num_h[i];
	    valSq = val*val;
            sum0 += valSq;
            num_h[i] = valSq;
         }

         // Determine highest power of 2 less than n
         int nReduced = (int)(pow(2.0,floor(log2(n))) + 0.5);
         // note: 0.5 added to make sure it casts to the correct int value
	 UTIL_CHECK(2 * nReduced >= n);
	 UTIL_CHECK(nReduced <= n);

         // Find sum on host using a binary tree
         // (numerical round-off error should match that from GPU sum)
         Timer timerCPU;
         if (verbose() > 0) {
            timerCPU.start();
         }

         for (int i = 0; i < nReduced; i++) {
            if (i + nReduced < n) {
               addEq(num_h[i], num_h[nReduced+i]);
            }
         }
         nReduced /= 2;
         for ( ; nReduced >= 1; nReduced /= 2) {
            for (int i = 0; i < nReduced; i++) {
               addEq(num_h[i], num_h[nReduced+i]);
            }
         }
         cudaReal sumCPU;
	 assign(sumCPU, num_h[0]);
         if (verbose() > 0) {
            timerCPU.stop();
            #if 0
            Log::file() << "sum0  : " << Dbl(sum0) << "\n";
            Log::file() << "sumCPU: " << Dbl(sumCPU) << std::endl;
            #endif
         }

         // Call kernel wrapper to calculate sum on GPU
         Timer timerGPU;
         if (verbose() > 0) {
            timerGPU.start();
         }
	 cudaReal sumGPU = Reduce::sumSq(num_d);
	 cudaReal diff = sumGPU - sumCPU;

         // Check answer
         if (verbose() > 0) {
            timerGPU.stop();
            Log::file() 
               << "CPU wall time: " << Dbl(timerCPU.time()) << "\n"
               << "GPU wall time: " << Dbl(timerGPU.time()) << "\n"
               << "Sum on CPU:    " << Dbl(sumCPU) << "\n"
               << "Sum on GPU:    " << Dbl(sumGPU) << "\n"
               << "Difference:    " << std::abs(diff) << std::endl;
         }

         // Check that error is at least 5 (10) orders of magnitude smaller
         // than the value of the sum for single (double) precision data
         cudaReal relDiff = std::abs(diff) / std::abs(sumCPU);
         TEST_ASSERT( relDiff < tolerance_);
      }
   }

   /// Test sum of an array of cudaComplex elements.
   void testSumSqComplex()
   {
      printMethod(TEST_FUNC);

      IntVec<3> nVals;
      nVals[0] = 109857634; // large array
      nVals[1] = 1934896;   // medium array
      nVals[2] = 50022;     // small array

      for (int j = 0; j < 3; j++) {

         int n = nVals[j];

         if (verbose() > 0) {
            Log::file() << std::endl << "n = " << n << std::endl;
         }

         // Generate random test data,
         // normally distributed about 0.001 with stdev = 1.0
         DeviceArray<cudaReal> num_dr(2*n);
         HostDArray<cudaReal>  num_hr(2*n);
         rand_.normal(num_dr, (cudaReal)1.0, (cudaReal)0.001);
         num_hr = num_dr;

          // Copy data to cudaComplex arrays num_h and num_d
         HostDArray<cudaComplex>  num_h(n);
         DeviceArray<cudaComplex> num_d(n);
	 cudaReal valx, valy;
	 cudaComplex valSq;
         cudaComplex sum0 = makeComplex(0.0, 0.0);
         for (int i = 0; i < n; ++i) {
            num_h[i].x = num_hr[2*i];
            num_h[i].y = num_hr[2*i + 1];
         }
         num_d = num_h;

         // Convert num_h to array of complex square values
         for (int i = 0; i < n; ++i) {
	    valx = num_h[i].x;
	    valy = num_h[i].y;
            valSq.x = (valx * valx) - (valy * valy);
            valSq.y = 2.0 * valx * valy;
            sum0.x += valSq.x;
            sum0.y += valSq.y;
            num_h[i].x = valSq.x;
            num_h[i].y = valSq.y;
         }

         // Determine highest power of 2 less than n
         int nReduced = (int)(pow(2.0,floor(log2(n))) + 0.5);
         // note: 0.5 added to make sure it casts to the correct int value
	 UTIL_CHECK(2 * nReduced >= n);
	 UTIL_CHECK(nReduced <= n);

         // Find sum on host using a binary tree
         // (numerical round-off error should match that from GPU sum)
         Timer timerCPU;
         if (verbose() > 0) {
            timerCPU.start();
         }

         for (int i = 0; i < nReduced; i++) {
            if (i + nReduced < n) {
               addEq(num_h[i], num_h[nReduced+i]);
            }
         }
         nReduced /= 2;
         for ( ; nReduced >= 1; nReduced /= 2) {
            for (int i = 0; i < nReduced; i++) {
               addEq(num_h[i], num_h[nReduced+i]);
            }
         }
         std::complex<cudaReal> sumCPU;
	 assign(sumCPU, num_h[0]);
         if (verbose() > 0) {
            timerCPU.stop();
            Log::file() << "CPU wall time: " << Dbl(timerCPU.time())
                        << std::endl;
            #if 0
            Log::file() << "sum0  : "
                        << Dbl(sum0.x) << "  " << Dbl(sum0.y) << "\n";
            Log::file() << "sumCPU: "
                        << Dbl(sumCPU.real()) << "  " << Dbl(sumCPU.imag())
                        << std::endl;
            #endif
         }

         // Call kernel wrapper to calculate sum on GPU
         Timer timerGPU;
         if (verbose() > 0) {
            timerGPU.start();
         }
	 std::complex<cudaReal> sumGPU = Reduce::sumSq(num_d);
	 std::complex<cudaReal> diff = sumGPU - sumCPU;

         // Check answer
         if (verbose() > 0) {
            timerGPU.stop();
            Log::file() 
               << "GPU wall time: " << Dbl(timerGPU.time()) << "\n"
               << "Sum on CPU:    "
               << Dbl(sumCPU.real()) << "  " << Dbl(sumCPU.imag()) << "\n"
               << "Sum on GPU:    "
               << Dbl(sumGPU.real()) << "  " << Dbl(sumGPU.imag()) << "\n"
               << "Difference:    " << std::abs(diff) << std::endl;
         }

         // Check that error is at least 5 (10) orders of magnitude smaller
         // than the value of the sum for single (double) precision data
         cudaReal relDiff = std::abs(diff) / std::abs(sumCPU);
         TEST_ASSERT( relDiff < tolerance_);
      }
   }

   /// Test inner product of two real arrays.
   void testInnerProduct()
   {
      printMethod(TEST_FUNC);

      IntVec<3> nVals;
      nVals[0] = 109857634; // large array
      nVals[1] = 1934896;   // medium array
      nVals[2] = 50022;     // small array

      for (int j = 0; j < 3; j++) {

         int n = nVals[j];

         if (verbose() > 0) {
            Log::file() << std::endl << "n = " << n << std::endl;
         }

         // Generate random test data, normally distributed
         DeviceArray<cudaReal> a(n), b(n);
         rand_.normal(a, (cudaReal)2.0, (cudaReal)0.5);
         rand_.normal(b, (cudaReal)1.0, (cudaReal)2.0);

         // Copy test data to host
         HostDArray<cudaReal> a_h(n), b_h(n);
         a_h = a;
         b_h = b;

         // Determine highest power of 2 less than n
         int nReduced = (int)(pow(2.0,floor(log2(n))) + 0.5);
         // note: 0.5 added to make sure it casts to the correct int value

         // Find inner product on host using a binary tree
         // (numerical round-off error should match that from GPU sum)
         Timer timerCPU;
         if (verbose() > 0) {
            timerCPU.start();
         }

         for (int i = 0; i < nReduced; i++) {
            a_h[i] *= b_h[i];
            if (i + nReduced < n) {
               a_h[i] += a_h[nReduced+i] * b_h[nReduced+i];
            }
         }
         nReduced /= 2;
         for ( ; nReduced >= 1; nReduced /= 2) {
            for (int i = 0; i < nReduced; i++) {
               a_h[i] += a_h[nReduced+i];
            }
         }
         cudaReal ipCPU = a_h[0];
         if (verbose() > 0) {
            timerCPU.stop();
            Log::file() << "CPU wall time:     " << Dbl(timerCPU.time())
                        << std::endl;
         }

         // Call kernel wrapper to calculate inner product on GPU
         Timer timerGPU;
         if (verbose() > 0) {
            timerGPU.start();
         }
         cudaReal ipGPU = Reduce::innerProduct(a, b);

         // Check answer
         if (verbose() > 0) {
            timerGPU.stop();
            Log::file()
                 << "GPU wall time:     " << Dbl(timerGPU.time()) << "\n"
                 << "Inner prod on CPU: " << Dbl(ipCPU) << "\n"
                 << "Inner prod on GPU: " << Dbl(ipGPU) << "\n"
                 << "Difference:        " << fabs(ipCPU - ipGPU) << "\n"
                 << std::endl;
         }

         // Check that error is at least 5 (10) orders of magnitude smaller
         // than the value for single (double) precision data
         TEST_ASSERT((fabs(ipCPU - ipGPU) / ipCPU) < tolerance_);
      }
   }

   /// Test max of an array of real elements.
   void testMax()
   {
      printMethod(TEST_FUNC);

      IntVec<3> nVals;
      nVals[0] = 50022;     // small array
      nVals[1] = 1934896;   // medium array
      nVals[2] = 109857634; // large array

      for (int j = 0; j < 3; j++) {

         int n = nVals[j];

         if (verbose() > 0) {
            Log::file() << std::endl << "n = " << n << std::endl;
         }

         // Generate random test data,
         // normally distributed about 7.0 with stdev = 3
         DeviceArray<cudaReal> num(n);
         rand_.normal(num, (cudaReal)3.0, (cudaReal)7.0);

         // Copy test data to host
         HostDArray<cudaReal> num_h(n);
         num_h = num;

         // Find max on host
         Timer timerCPU;
         if (verbose() > 0) {
            timerCPU.start();
         }
         cudaReal maxCPU = num_h[0];
         for (int i = 1; i < n; i++) {
            if (num_h[i] > maxCPU) maxCPU = num_h[i];
         }
         if (verbose() > 0) {
            timerCPU.stop();
            Log::file() << "CPU wall time: " << Dbl(timerCPU.time())
                        << std::endl;
         }

         // Call kernel wrapper to calculate max on GPU
         Timer timerGPU;
         if (verbose() > 0) {
            timerGPU.start();
         }
         cudaReal maxGPU = Reduce::max(num);

         // Check answer
         if (verbose() > 0) {
            timerGPU.stop();
            Log::file()
                  << "GPU wall time: " << Dbl(timerGPU.time()) << "\n"
                  << "Max on CPU:    " << Dbl(maxCPU) << "\n"
                  << "Max on GPU:    " << Dbl(maxGPU) << "\n"
                  << "Difference:    " << fabs(maxCPU - maxGPU) << "\n"
                  << std::endl;
         }
         TEST_ASSERT((fabs(maxCPU - maxGPU)) < tolerance_);
      }
   }

   /// Test max of an array of real elements.
   void testMaxSlice()
   {
      printMethod(TEST_FUNC);

      IntVec<3> nVals;
      nVals[0] = 734;     // small array
      nVals[1] = 4378;    // medium array

      for (int j = 0; j < 2; j++) {

         int n = nVals[j];
	 int begin = 1;
	 int end = n;

         if (verbose() > 0) {
            Log::file() << std::endl << "n = " << n << std::endl;
         }

         // Generate random test data,
         // normally distributed about 7.0 with stdev = 3
         DeviceArray<cudaReal> num(n);
         rand_.normal(num, (cudaReal)3.0, (cudaReal)7.0);

         // Copy test data to host
         HostDArray<cudaReal> num_h(n);
         num_h = num;

         // Find max on host
         Timer timerCPU;
         if (verbose() > 0) {
            timerCPU.start();
         }
         cudaReal maxCPU = num_h[begin];
         for (int i = begin; i < end; i++) {
            if (num_h[i] > maxCPU) {
               maxCPU = num_h[i];
            }
         }
         if (verbose() > 0) {
            timerCPU.stop();
            Log::file() << "CPU wall time: " << Dbl(timerCPU.time())
                        << std::endl;
         }

         Timer timerGPU;
         if (verbose() > 0) {
            timerGPU.start();
         }

         // Evaluate on GPU
         cudaReal maxGPU = Reduce::max(num, begin, end);

         // Check answer
         if (verbose() > 0) {
            timerGPU.stop();
            Log::file()
                  << "GPU wall time: " << Dbl(timerGPU.time()) << "\n"
                  << "Max on CPU:    " << Dbl(maxCPU) << "\n"
                  << "Max on GPU:    " << Dbl(maxGPU) << "\n"
                  << "Difference:    " << fabs(maxCPU - maxGPU) << "\n"
                  << std::endl;
         }
         TEST_ASSERT((fabs(maxCPU - maxGPU)) < tolerance_);
      }
   }

   /// Test maxAbs for a real array.
   void testMaxAbs()
   {
      printMethod(TEST_FUNC);

      IntVec<3> nVals;
      nVals[0] = 50022;     // small array
      nVals[1] = 1934896;   // medium array
      nVals[2] = 109857634; // large array

      for (int j = 0; j < 3; j++) {

         int n = nVals[j];

         if (verbose() > 0) {
            Log::file() << std::endl << "n = " << n << std::endl;
         }

         // Generate random test data,
         // normally distributed about -1.0 with stdev = 3
         DeviceArray<cudaReal> num(n);
         rand_.normal(num, (cudaReal)3.0, (cudaReal)-1.0);

         // Copy test data to host
         HostDArray<cudaReal> num_h(n);
         num_h = num;

         // Find max on host
         Timer timerCPU;
         if (verbose() > 0) {
            timerCPU.start();
         }
         cudaReal maxCPU = fabs(num_h[0]);
         cudaReal val;
         for (int i = 1; i < n; i++) {
            val = fabs(num_h[i]);
            if (val > maxCPU) maxCPU = val;
         }
         if (verbose() > 0) {
            timerCPU.stop();
            Log::file() << "CPU wall time: " << Dbl(timerCPU.time())
                        << std::endl;
         }

         // Call kernel wrapper to calculate sum on GPU
         Timer timerGPU;
         if (verbose() > 0) {
            timerGPU.start();
         }
         cudaReal maxGPU = Reduce::maxAbs(num);

         // Check answer
         if (verbose() > 0) {
            timerGPU.stop();
            Log::file() << "GPU wall time: " << Dbl(timerGPU.time()) << "\n"
                        << "Max on CPU:    " << Dbl(maxCPU) << "\n"
                        << "Max on GPU:    " << Dbl(maxGPU) << "\n"
                        << "Difference:    " << fabs(maxCPU - maxGPU) << "\n"
                        << std::endl;
         }
         TEST_ASSERT((fabs(maxCPU - maxGPU)) < tolerance_);
      }
   }

   /// Test min for a real array.
   void testMin()
   {
      printMethod(TEST_FUNC);

      IntVec<3> nVals;
      nVals[0] = 50022;     // small array
      nVals[1] = 1934896;   // medium array
      nVals[2] = 109857634; // large array

      for (int j = 0; j < 3; j++) {

         int n = nVals[j];

         if (verbose() > 0) {
            Log::file() << std::endl << "n = " << n << std::endl;
         }

         // Generate random test data,
         // normally distributed about 7.0 with stdev = 3
         DeviceArray<cudaReal> num(n);
         rand_.normal(num, (cudaReal)3.0, (cudaReal)7.0);

         // Copy test data to host
         HostDArray<cudaReal> num_h(n);
         num_h = num;

         // Find min on host
         Timer timerCPU;
         if (verbose() > 0) {
            timerCPU.start();
         }
         cudaReal minCPU = num_h[0];
         for (int i = 1; i < n; i++) {
            if (num_h[i] < minCPU) minCPU = num_h[i];
         }
         if (verbose() > 0) {
            timerCPU.stop();
            Log::file() << "CPU wall time: " << Dbl(timerCPU.time())
                        << std::endl;
         }

         // Call kernel wrapper to calculate sum on GPU
         Timer timerGPU;
         if (verbose() > 0) {
            timerGPU.start();
         }
         cudaReal minGPU = Reduce::min(num);

         // Check answer
         if (verbose() > 0) {
            timerGPU.stop();
            Log::file()
                  << "GPU wall time: " << Dbl(timerGPU.time()) << "\n"
                  << "Min on CPU:    " << Dbl(minCPU) << "\n"
                  << "Min on GPU:    " << Dbl(minGPU) << "\n"
                  << "Difference:    " << fabs(minCPU - minGPU) << "\n"
                  << std::endl;
         }
         TEST_ASSERT((fabs(minCPU - minGPU)) < tolerance_);
      }
   }

   /// Test minAbs for a real array.
   void testMinAbs()
   {
      printMethod(TEST_FUNC);

      IntVec<3> nVals;
      nVals[0] = 50022;     // small array
      nVals[1] = 1934896;   // medium array
      nVals[2] = 109857634; // large array

      for (int j = 0; j < 3; j++) {

         int n = nVals[j];

         if (verbose() > 0) {
            Log::file() << std::endl << "n = " << n << std::endl;
         }

         // Random data -
         // normal distribution about -1.0 with stdev = 3
         DeviceArray<cudaReal> num(n);
         rand_.normal(num, (cudaReal)3.0, (cudaReal)-1.0);

         // Copy test data to host
         HostDArray<cudaReal> num_h(n);
         num_h = num;

         // Find min on host
         Timer timerCPU;
         if (verbose() > 0) {
            timerCPU.start();
         }
         cudaReal minCPU = fabs(num_h[0]);
         cudaReal val;
         for (int i = 1; i < n; i++) {
            val = fabs(num_h[i]);
            if (val < minCPU) minCPU = val;
         }
         if (verbose() > 0) {
            timerCPU.stop();
            Log::file() << "CPU wall time: " << Dbl(timerCPU.time())
                        << std::endl;
         }

         // Call kernel wrapper to calculate min on GPU
         Timer timerGPU;
         if (verbose() > 0) {
            timerGPU.start();
         }
         cudaReal minGPU = Reduce::minAbs(num);

         // Check answer
         if (verbose() > 0) {
            timerGPU.stop();
            Log::file()
                 << "GPU wall time: " << Dbl(timerGPU.time()) << "\n"
                 << "Min on CPU:    " << Dbl(minCPU) << "\n"
                 << "Min on GPU:    " << Dbl(minGPU) << "\n"
                 << "Difference:    " << fabs(minCPU - minGPU) << "\n"
                 << std::endl;
         }
         TEST_ASSERT((fabs(minCPU - minGPU)) < tolerance_);
      }
   }

};

TEST_BEGIN(CudaReduceTest)
TEST_ADD(CudaReduceTest, testSum)
TEST_ADD(CudaReduceTest, testSumSlice)
TEST_ADD(CudaReduceTest, testSumComplex)
TEST_ADD(CudaReduceTest, testSumComplexSlice)
TEST_ADD(CudaReduceTest, testSumSq)
TEST_ADD(CudaReduceTest, testSumSqComplex)
TEST_ADD(CudaReduceTest, testInnerProduct)
TEST_ADD(CudaReduceTest, testMax)
TEST_ADD(CudaReduceTest, testMaxSlice)
TEST_ADD(CudaReduceTest, testMaxAbs)
TEST_ADD(CudaReduceTest, testMin)
TEST_ADD(CudaReduceTest, testMinAbs)
TEST_END(CudaReduceTest)

#endif

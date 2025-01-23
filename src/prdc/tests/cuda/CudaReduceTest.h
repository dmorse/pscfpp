#ifndef PRDC_CUDA_REDUCE_TEST_H
#define PRDC_CUDA_REDUCE_TEST_H

#include <test/UnitTest.h>
#include <test/UnitTestRunner.h>

#include <prdc/cuda/types.h>
#include <prdc/cuda/Reduce.h>
#include <pscf/cuda/CudaRandom.h>
#include <pscf/cuda/DeviceArray.h>
#include <pscf/cuda/HostDArray.h>
#include <util/format/Dbl.h>
#include <util/misc/Timer.h>

#include <cstdlib>
#include <cmath>

using namespace Util;
using namespace Pscf;
using namespace Pscf::Prdc;
using namespace Pscf::Prdc::Cuda;

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
   CudaRandom rand_;

public:

   void setUp()
   {  
      setVerbose(0); 
      rand_.setSeed(0);
   }

   void tearDown()
   {}

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

         // Generate test data, normally distributed about 0.5 with stdev = 2
         DeviceArray<cudaReal> num(n);
         rand_.normal(num, (cudaReal)2.0, (cudaReal)0.5);

         // Copy test data to host
         HostDArray<cudaReal> num_h(n);
         num_h = num;

         // Determine highest power of 2 less than n
         int nReduced = (int)(pow(2.0,floor(log2(n))) + 0.5); 
         // note: 0.5 added to make sure it casts to the correct int value

         // Find sum on host using a binary tree 
         // (numerical round-off error should match that from the GPU summation)
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
            Log::file() << "GPU wall time: " << Dbl(timerGPU.time()) << "\n"
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

         // Generate test data, normally distributed about 7.0 with stdev = 3
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
         
         // Call kernel wrapper to calculate sum on GPU
         Timer timerGPU;
         if (verbose() > 0) {
            timerGPU.start();
         }
         cudaReal maxGPU = Reduce::max(num);
         
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

         // Generate test data, normally distributed about -1.0 with stdev = 3
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

         // Generate test data, normally distributed about 7.0 with stdev = 3
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
            Log::file() << "GPU wall time: " << Dbl(timerGPU.time()) << "\n"
                        << "Min on CPU:    " << Dbl(minCPU) << "\n"
                        << "Min on GPU:    " << Dbl(minGPU) << "\n"
                        << "Difference:    " << fabs(minCPU - minGPU) << "\n"
                        << std::endl;
         }
         TEST_ASSERT((fabs(minCPU - minGPU)) < tolerance_);
      }
   }

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

         // Generate test data, normally distributed about -1.0 with stdev = 3
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
         
         // Call kernel wrapper to calculate sum on GPU
         Timer timerGPU;
         if (verbose() > 0) {
            timerGPU.start();
         }
         cudaReal minGPU = Reduce::minAbs(num);
         
         // Check answer
         if (verbose() > 0) {
            timerGPU.stop();
            Log::file() << "GPU wall time: " << Dbl(timerGPU.time()) << "\n"
                        << "Min on CPU:    " << Dbl(minCPU) << "\n"
                        << "Min on GPU:    " << Dbl(minGPU) << "\n"
                        << "Difference:    " << fabs(minCPU - minGPU) << "\n"
                        << std::endl;
         }
         TEST_ASSERT((fabs(minCPU - minGPU)) < tolerance_);
      }
   }

   void testInnerProduct() 
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

         // Generate test data, normally distributed
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
         // (numerical round-off error should match that from the GPU summation)
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
            Log::file() << "GPU wall time:     " << Dbl(timerGPU.time()) << "\n"
                        << "Inner prod on CPU: " << Dbl(ipCPU) << "\n"
                        << "Inner prod on GPU: " << Dbl(ipGPU) << "\n"
                        << "Difference:        " << fabs(ipCPU - ipGPU) << "\n"
                        << std::endl;
         }

         // Check that error is at least 5 (10) orders of magnitude smaller 
         // than the value of the inner prod for single (double) precision data
         TEST_ASSERT((fabs(ipCPU - ipGPU) / ipCPU) < tolerance_);
      }
   }

};

TEST_BEGIN(CudaReduceTest)
TEST_ADD(CudaReduceTest, testSum)
TEST_ADD(CudaReduceTest, testMax)
TEST_ADD(CudaReduceTest, testMaxAbs)
TEST_ADD(CudaReduceTest, testMin)
TEST_ADD(CudaReduceTest, testMinAbs)
TEST_ADD(CudaReduceTest, testInnerProduct)
TEST_END(CudaReduceTest)

#endif

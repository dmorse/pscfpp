#ifndef PRDC_CUDA_REDUCE_TEST_H
#define PRDC_CUDA_REDUCE_TEST_H

#include <test/UnitTest.h>
#include <test/UnitTestRunner.h>

#include <pscf/cuda/Reduce.h>
#include <pscf/cuda/CudaRandom.h>
#include <pscf/cuda/KernelWrappers.h>
#include <pscf/cuda/GpuTypes.h>
#include <util/format/Dbl.h>
#include <util/misc/Timer.h>

#include <cstdlib>
#include <cmath>

using namespace Util;
using namespace Pscf;

class CudaReduceTest : public UnitTest
{

private:

   // Error tolerance for array equality
   #ifdef SINGLE_PRECISION
   constexpr static float tolerance_ = 1E-5;
   #else
   #ifdef DOUBLE_PRECISION
   constexpr static double tolerance_ = 1E-10;
   #endif
   #endif

   // Random number generator on the GPU
   CudaRandom rand_;

public:

   void setUp()
   {  
      setVerbose(1); 
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
         rand_.normal(num.cArray(), n, (cudaReal)2.0, (cudaReal)0.5);

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
         cudaReal sum = Reduce::sum(num);
         
         // Check answer
         if (verbose() > 0) {
            timerGPU.stop();
            Log::file() << "GPU wall time: " << Dbl(timerGPU.time()) << "\n"
                        << "Sum on CPU:    " << Dbl(num_h[0]) << "\n"
                        << "Sum on GPU:    " << Dbl(sum) << "\n"
                        << "Difference:    " << fabs(sum - num_h[0]) << "\n"
                        << std::endl;
         }
         TEST_ASSERT((fabs(sum - num_h[0]) / sum) < tolerance_);

         // Check answer against old version of the code
         Timer timerGPU2;
         if (verbose() > 0) {
            timerGPU2.start();
         }
         cudaReal sumOld = gpuSum(num.cArray(), n);

         // Check answer
         if (verbose() > 0) {
            timerGPU2.stop();
            Log::file() << "Old version:\n"
                        << "GPU wall time: " << Dbl(timerGPU2.time()) << "\n"
                        << "Sum on GPU:    " << Dbl(sumOld) << "\n"
                        << "Difference:    " << fabs(sum - sumOld) << "\n"
                        << std::endl;
         }
         TEST_ASSERT((fabs(sum - sumOld) / sum) < tolerance_);
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
         rand_.normal(num.cArray(), n, (cudaReal)3.0, (cudaReal)7.0);

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

         // Check answer against old version of the code
         Timer timerGPU2;
         if (verbose() > 0) {
            timerGPU2.start();
         }
         cudaReal maxOld = gpuMax(num.cArray(), n);

         // Check answer
         if (verbose() > 0) {
            timerGPU2.stop();
            Log::file() << "Old version:\n"
                        << "GPU wall time: " << Dbl(timerGPU2.time()) << "\n"
                        << "Max on GPU:    " << Dbl(maxOld) << "\n"
                        << "Difference:    " << fabs(maxGPU - maxOld) << "\n"
                        << std::endl;
         }
         TEST_ASSERT((fabs(maxGPU - maxOld)) < tolerance_);
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
         rand_.normal(num.cArray(), n, (cudaReal)3.0, (cudaReal)-1.0);

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

         // Check answer against old version of the code
         Timer timerGPU2;
         if (verbose() > 0) {
            timerGPU2.start();
         }
         cudaReal maxOld = gpuMaxAbs(num.cArray(), n);

         // Check answer
         if (verbose() > 0) {
            timerGPU2.stop();
            Log::file() << "Old version:\n"
                        << "GPU wall time: " << Dbl(timerGPU2.time()) << "\n"
                        << "Max on GPU:    " << Dbl(maxOld) << "\n"
                        << "Difference:    " << fabs(maxGPU - maxOld) << "\n"
                        << std::endl;
         }
         TEST_ASSERT((fabs(maxGPU - maxOld)) < tolerance_);
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
         rand_.normal(num.cArray(), n, (cudaReal)3.0, (cudaReal)7.0);

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

         // Check answer against old version of the code
         Timer timerGPU2;
         if (verbose() > 0) {
            timerGPU2.start();
         }
         cudaReal minOld = gpuMin(num.cArray(), n);

         // Check answer
         if (verbose() > 0) {
            timerGPU2.stop();
            Log::file() << "Old version:\n"
                        << "GPU wall time: " << Dbl(timerGPU2.time()) << "\n"
                        << "Min on GPU:    " << Dbl(minOld) << "\n"
                        << "Difference:    " << fabs(minGPU - minOld) << "\n"
                        << std::endl;
         }
         TEST_ASSERT((fabs(minGPU - minOld)) < tolerance_);
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
         rand_.normal(num.cArray(), n, (cudaReal)3.0, (cudaReal)-1.0);

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

         // Check answer against old version of the code
         Timer timerGPU2;
         if (verbose() > 0) {
            timerGPU2.start();
         }
         cudaReal minOld = gpuMinAbs(num.cArray(), n);

         // Check answer
         if (verbose() > 0) {
            timerGPU2.stop();
            Log::file() << "Old version:\n"
                        << "GPU wall time: " << Dbl(timerGPU2.time()) << "\n"
                        << "Min on GPU:    " << Dbl(minOld) << "\n"
                        << "Difference:    " << fabs(minGPU - minOld) << "\n"
                        << std::endl;
         }
         TEST_ASSERT((fabs(minGPU - minOld)) < tolerance_);
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
         rand_.normal(a.cArray(), n, (cudaReal)2.0, (cudaReal)0.5);
         rand_.normal(b.cArray(), n, (cudaReal)1.0, (cudaReal)2.0);

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
         cudaReal ip = Reduce::innerProduct(a, b);
         
         // Check answer
         if (verbose() > 0) {
            timerGPU.stop();
            Log::file() << "GPU wall time:     " << Dbl(timerGPU.time()) << "\n"
                        << "Inner prod on CPU: " << Dbl(a_h[0]) << "\n"
                        << "Inner prod on GPU: " << Dbl(ip) << "\n"
                        << "Difference:        " << fabs(ip - a_h[0]) << "\n"
                        << std::endl;
         }
         TEST_ASSERT((fabs(ip - a_h[0]) / ip) < tolerance_);

         // Check answer against old version of the code
         Timer timerGPU2;
         if (verbose() > 0) {
            timerGPU2.start();
         }
         cudaReal ipOld = gpuInnerProduct(a.cArray(), b.cArray(), n);

         // Check answer
         if (verbose() > 0) {
            timerGPU2.stop();
            Log::file() << "Old version:\n"
                        << "GPU wall time:     " << Dbl(timerGPU2.time()) << "\n"
                        << "Inner prod on GPU: " << Dbl(ipOld) << "\n"
                        << "Difference:        " << fabs(ip - ipOld) << "\n"
                        << std::endl;
         }
         TEST_ASSERT((fabs(ip - ipOld) / ip) < tolerance_);
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

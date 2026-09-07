#include "CudaVecRandom.h"
#include "ThreadArray.h"
#include <util/global.h>
#include <sys/time.h>
#include <string>

namespace Pscf {

   using namespace Util;

   // Anonymous namesapce, for functions that are only used in this file
   namespace {

      /*
      * Linear array transformation a[i] => c*a[i] + s (float).
      */
      __global__
      void _linearScale(float* a, float c, float s, int n)
      {
         int nThreads = blockDim.x * gridDim.x;
         int startID = blockIdx.x * blockDim.x + threadIdx.x;
         for (int i = startID; i < n; i += nThreads) {
            a[i] = a[i] * c + s;
         }
      }

      /*
      * Linear array transformation a[i] => c*a[i] + s (double).
      */
      __global__
      void _linearScale(double* a, double c, double s, int n)
      {
         int nThreads = blockDim.x * gridDim.x;
         int startID = blockIdx.x * blockDim.x + threadIdx.x;
         for (int i = startID; i < n; i += nThreads) {
            a[i] = a[i] * c + s;
         }
      }

      /*
      * Linear array transformation a[i] => c*a[i] + s (float).
      */
      void linearScale(DeviceArray<float>& a, float c, float s)
      {
         const int n = a.capacity();

         // GPU resources
         int nBlocks, nThreads;
         ThreadArray::setThreadsLogical(n, nBlocks, nThreads);

         // Launch kernel
         _linearScale<<<nBlocks, nThreads>>>(a.cArray(), c, s, n);
      }

      /*
      * Linear array transformation a[i] => c*a[i] + s (double).
      */
      void linearScale(DeviceArray<double>& a, double c, double s)
      {
         const int n = a.capacity();

         // GPU resources
         int nBlocks, nThreads;
         ThreadArray::setThreadsLogical(n, nBlocks, nThreads);

         // Launch kernel
         _linearScale<<<nBlocks, nThreads>>>(a.cArray(), c, s, n);
      }

   } // end anonymous namespace

   // Public class member functions

   /*
   * Constructor.
   */
   CudaVecRandom::CudaVecRandom()
    : gen_(),
      seed_(0),
      isInitialized_(false)
   {
      // Create pseudo-random number generator on gpu
      curandStatus_t status;
      status = curandCreateGenerator(&gen_, CURAND_RNG_PSEUDO_DEFAULT);
      errorCheck(status);
   }

   /*
   * Destructor.
   */
   CudaVecRandom::~CudaVecRandom()
   {}

   /*
   * Set random random number generator seed.
   *
   * \param seed value for random seed (private member variable seed)
   */
   void CudaVecRandom::setSeed(unsigned long long seed)
   {
      if (seed == 0) {
         timeval time;
         gettimeofday(&time, NULL);
         seed_ = time.tv_sec + 1123*time.tv_usec;
      } else {
         seed_ = seed;
      }
      curandStatus_t status;
      status = curandSetPseudoRandomGeneratorSeed(gen_, seed_);
      errorCheck(status);

      isInitialized_ = true;
   }

   /*
   * Populate array with uniform random floats in (0, 1].
   */
   void CudaVecRandom::uniform(DeviceArray<float>& data)
   {
      const int n = data.capacity();
      UTIL_CHECK(n > 0);
      if (!isInitialized_) {
         setSeed(0);
      }

      curandStatus_t status;
      status = curandGenerateUniform(gen_, data.cArray(), n);
      errorCheck(status);
   }

   /*
   * Populate array with uniform random doubles in (0, 1].
   */
   void CudaVecRandom::uniform(DeviceArray<double>& data)
   {
      const int n = data.capacity();
      UTIL_CHECK(n > 0);
      if (!isInitialized_) {
         setSeed(0);
      }

      curandStatus_t status;
      status = curandGenerateUniformDouble(gen_, data.cArray(), n);
      errorCheck(status);
   }

   /*
   * Populate array with uniform random floats in (min, max].
   */
   void CudaVecRandom::uniform(DeviceArray<float>& data, 
                            float min, float max)
   {
      UTIL_CHECK(max > min);
      uniform(data);
      linearScale(data, max - min, min);
   }

   /*
   * Populate array with uniform random doubles in (min, max].
   */
   void CudaVecRandom::uniform(DeviceArray<double>& data,
                            double min, double max)
   {
      UTIL_CHECK(max > min);
      uniform(data);
      linearScale(data, max - min, min);
   }

   /*
   * Populate array with normal-distributed random floats.
   */
   void CudaVecRandom::normal(DeviceArray<float>& data,
                           float stddev, float mean)
   {
      UTIL_CHECK(data.capacity() > 0);
      if (!isInitialized_) {
         setSeed(0);
      }

      int n = data.capacity();
      if (n % 2 == 1) {
         UTIL_THROW("normal() requires array size to be an even number.");
      }

      curandStatus_t status;
      status = curandGenerateNormal(gen_, data.cArray(), n, mean, stddev);
      errorCheck(status);
   }

   /*
   * Populate array with normal-distributed random doubles.
   */
   void CudaVecRandom::normal(DeviceArray<double>& data,
                           double stddev, double mean)
   {
      UTIL_CHECK(data.capacity() > 0);
      if (!isInitialized_) {
         setSeed(0);
      }

      int n = data.capacity();
      if (n % 2 == 1) {
         UTIL_THROW("normal() requires array size to be an even number.");
      }

      curandStatus_t status;
      status = curandGenerateNormalDouble(gen_, data.cArray(),
                                          n, mean, stddev);
      errorCheck(status);
   }

   /*
   * Check generator error status. 
   *
   * If not success, print info and throw Exception.
   */
   void CudaVecRandom::errorCheck(curandStatus_t const & error)
   {
      if (error == CURAND_STATUS_SUCCESS) {
         return;
      } else {
         std::string errString;
         switch (error)
         {
            default:
               errString = "UNKNOWN";
               break;
            case CURAND_STATUS_VERSION_MISMATCH:
               errString = "CURAND_STATUS_VERSION_MISMATCH";
               break;
            case CURAND_STATUS_NOT_INITIALIZED:
               errString = "CURAND_STATUS_NOT_INITIALIZED";
               break;
            case CURAND_STATUS_ALLOCATION_FAILED:
               errString = "CURAND_STATUS_ALLOCATION_FAILED";
               break;
            case CURAND_STATUS_TYPE_ERROR:
               errString = "CURAND_STATUS_TYPE_ERROR";
               break;
            case CURAND_STATUS_OUT_OF_RANGE:
               errString = "CURAND_STATUS_OUT_OF_RANGE";
               break;
            case CURAND_STATUS_LENGTH_NOT_MULTIPLE:
               errString = "CURAND_STATUS_LENGTH_NOT_MULTIPLE";
               break;
            case CURAND_STATUS_DOUBLE_PRECISION_REQUIRED:
               errString = "CURAND_STATUS_DOUBLE_PRECISION_REQUIRED";
               break;
            case CURAND_STATUS_LAUNCH_FAILURE:
               errString = "CURAND_STATUS_LAUNCH_FAILURE";
               break;
            case CURAND_STATUS_PREEXISTING_FAILURE:
               errString = "CURAND_STATUS_PREEXISTING_FAILURE";
               break;
            case CURAND_STATUS_INITIALIZATION_FAILED:
               errString = "CURAND_STATUS_INITIALIZATION_FAILED";
               break;
            case CURAND_STATUS_INTERNAL_ERROR:
               errString = "CURAND_STATUS_INTERNAL_ERROR";
               break;
            case CURAND_STATUS_ARCH_MISMATCH:
               errString = "CURAND_STATUS_ARCH_MISMATCH";
               break;
         }

         Log::file() << "CudaVecRandom error: " << errString << std::endl;
         UTIL_THROW("CudaVecRandom number generation failed.");
      }
   }

}

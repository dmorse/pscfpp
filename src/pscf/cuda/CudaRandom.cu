#include "CudaRandom.h"
#include <util/global.h>
#include <sys/time.h>
#include <string>

namespace Pscf {

   using namespace Util;

   /*
   * Constructor.
   */
   CudaRandom::CudaRandom()
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
   CudaRandom::~CudaRandom()
   {}

   /*
   * Sets of random seed, and initializes random number generator.
   *
   * \param seed value for random seed (private member variable seed)
   */
   void CudaRandom::setSeed(unsigned long long seed)
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
   * Populate array on device with random floats in (0, 1], uniform dist.
   */
   void CudaRandom::uniform(DeviceArray<float>& data)
   {
      UTIL_CHECK(data.capacity() > 0);
      if (!isInitialized_) {
         setSeed(0);
      }

      curandStatus_t status = curandGenerateUniform(gen_, data.cArray(), 
                                                    data.capacity());
      errorCheck(status);
   }

   /*
   * Populate array on device with random doubles in (0, 1], uniform dist.
   */
   void CudaRandom::uniform(DeviceArray<double>& data)
   {
      UTIL_CHECK(data.capacity() > 0);
      if (!isInitialized_) {
         setSeed(0);
      }
      
      curandStatus_t status = curandGenerateUniformDouble(gen_, data.cArray(), 
                                                          data.capacity());
      errorCheck(status);
   }

   /*
   * Populate array with normal-distributed random floats.
   */
   void CudaRandom::normal(DeviceArray<float>& data, float stddev, float mean)
   {
      UTIL_CHECK(data.capacity() > 0);
      if (!isInitialized_) {
         setSeed(0);
      }

      int n = data.capacity();
      if (n % 2 == 1) {
         UTIL_THROW("normal() requires array size to be an even number.");
      }
      
      curandStatus_t status = curandGenerateNormal(gen_, data.cArray(), 
                                                   n, mean, stddev);
      errorCheck(status);
   }

   /*
   * Populate array with normal-distributed random doubles.
   */
   void CudaRandom::normal(DeviceArray<double>& data, 
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

      curandStatus_t status = curandGenerateNormalDouble(gen_, data.cArray(), 
                                                         n, mean, stddev);
      errorCheck(status);
   }

   /*
   * Check generator error status. If not success, print info and throw error.
   */
   void CudaRandom::errorCheck(curandStatus_t const & error)
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

         Log::file() << "CudaRandom error: " << errString << std::endl;
         UTIL_THROW("CudaRandom number generation failed.");
      }
   }

}

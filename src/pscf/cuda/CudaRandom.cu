#include "CudaRandom.h"
#include "GpuTypes.h"

#include <util/global.h>

#include <curand.h>
#include <sys/time.h>

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
      UTIL_CHECK(status == CURAND_STATUS_SUCCESS);
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
      UTIL_CHECK(status == CURAND_STATUS_SUCCESS);
      isInitialized_ = true;
   }

   /*
   * Return uniformly distributed random number in [0,1]
   */
   double CudaRandom::uniform(cudaReal* data, int n)
   {
      if (!isInitialized_) {
         setSeed(0);
      }
      #ifdef SINGLE_PRECISION
      curandStatus_t gen_error 
           = curandGenerateUniform(gen_, data, int n);
      #else
      curandStatus_t gen_error 
            = curandGenerateUniformDouble(gen_, data, n);
      #endif
   }

}

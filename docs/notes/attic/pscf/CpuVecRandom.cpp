#include "CpuVecRandom.h"
#include <util/random/Random.h>
#include <util/containers/Array.h>
#include <util/global.h>

namespace Pscf {

   using namespace Util;

   /*
   * Default constructor.
   */
   CpuVecRandom::CpuVecRandom()
    : randomPtr_(nullptr)
   {}

   /*
   * Constructor - create association with a scalar RNG.
   */
   CpuVecRandom::CpuVecRandom(Random& random)
    : randomPtr_(&random)
   {}

   /*
   * Destructor.
   */
   CpuVecRandom::~CpuVecRandom()
   {}

   /*
   * Create an association after construction
   */
   void CpuVecRandom::associate(Random& random)
   {  randomPtr_ = &random; }

   /*
   * Populate array on device with random doubles in (0, 1], uniform dist.
   */
   void CpuVecRandom::uniform(Array<double>& data)
   {
      UTIL_CHECK(randomPtr_);
      UTIL_CHECK(data.capacity() > 0);
      const int n = data.capacity();
      for (int i = 0; i < n; ++i) {
         data[i] = randomPtr_->uniform();
      } 
   }

   /*
   * Populate array with random doubles uniform dist in (min, max].
   */
   void CpuVecRandom::uniform(Array<double>& data, double min, double max)
   {
      UTIL_CHECK(randomPtr_);
      UTIL_CHECK(data.capacity() > 0);
      UTIL_CHECK(max > min);
      const int n = data.capacity();
      for (int i = 0; i < n; ++i) {
         data[i] = randomPtr_->uniform(min, max);
      } 
   }

   /*
   * Populate array with normal-distributed random doubles.
   */
   void CpuVecRandom::normal(Array<double>& data, 
                             double stddev, double mean)
   {
      UTIL_CHECK(randomPtr_);
      UTIL_CHECK(data.capacity() > 0);
      const int n = data.capacity();
      for (int i = 0; i < n; ++i) {
         data[i] = mean + stddev * randomPtr_->gaussian();
      }
   }

} // namespace Pscf

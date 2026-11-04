#include "VecRandom.h"
#include <util/random/Random.h>
#include <util/containers/Array.h>
#include <util/global.h>

namespace Pscf {

   using namespace Util;

   /*
   * Default constructor.
   */
   VecRandom<CPT>::VecRandom()
    : randomPtr_(nullptr)
   {}

   /*
   * Constructor - create association with a scalar RNG.
   */
   VecRandom<CPT>::VecRandom(Util::Random& random)
    : randomPtr_(&random)
   {}

   /*
   * Destructor.
   */
   VecRandom<CPT>::~VecRandom()
   {}

   /*
   * Create an association with a scalar RNG after construction
   */
   void VecRandom<CPT>::associate(Util::Random& random)
   {  randomPtr_ = &random; }

   /*
   * Populate array on device with random doubles in (0, 1], uniform dist.
   */
   void VecRandom<CPT>::uniform(Util::Array<double>& data)
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
   void VecRandom<CPT>::uniform(Util::Array<double>& data, double min, double max)
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
   void VecRandom<CPT>::normal(Util::Array<double>& data, 
                               double stddev, double mean)
   {
      UTIL_CHECK(randomPtr_);
      UTIL_CHECK(data.capacity() > 0);
      const int n = data.capacity();
      for (int i = 0; i < n; ++i) {
         data[i] = mean + stddev * randomPtr_->gaussian();
      }
   }

   // Free functions

   /*
   * Associate vector and scalar random number generators, if needed.
   *
   * GPU-enabled code would instead do nothing.
   */
   template <>
   void linkVecRandom<CPT>(VecRandom<CPT>& vr, Random & sr)
   {  vr.associate(sr); }

   /*
   * Set vector RNG seed, if needed (do nothing for CPU code).
   *
   * GPU-enabled code would instead set the seed.
   */
   template <>
   void seedVecRandom<CPT>(VecRandom<CPT>& vr, long seed)
   {}


} // namespace Pscf

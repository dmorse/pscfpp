#ifndef PSCF_CUDA_RANDOM_H
#define PSCF_CUDA_RANDOM_H

#include "GpuTypes.h"

#include <curand.h>

/*
* PSCF Package - Polymer Self-Consistent Field 
*
* Copyright 2016 - 2023, The Regents of the University of Minnesota
* Distributed under the terms of the GNU General Public License.
*/

namespace Pscf {

   /**
   * CudaRandom GPU number generator.
   *
   * The generator may be seeded either by reading a seed from 
   * file, using the readParam() method, or by using setSeed()
   * to set or reset it explicitly. In either case, inputting 
   * a positive integer causes that value to be used as a seed,
   * but inputting a value of 0 causes the use of a seed that 
   * is generated from the system clock. 
   *
   * \ingroup Cuda_Module
   */
   class CudaRandom 
   {
 
   public:

      typedef unsigned long SeedType;

      /**
      * Constructor.
      */   
      CudaRandom();

      /**
      * Destructor.
      */   
      virtual ~CudaRandom();
   
      void setSeed(unsigned long long seed);

      /**
      * Return an array of uniform random floating point numbers.
      *   
      * range 0 <= x < 1.
      */
      double uniform(cudaReal* data, int n);
   
      /**
      * Returns value of random seed (private member variable seed_).
      *
      * \return value of random number generator seed.
      */
      long seed();
   
   private:

      /// GPU random number generator
      curandGenerator_t gen_;

      /// Seed value.
      unsigned long long seed_;

      /// Has a seed been set by readParam() or setSeed()?
      bool isInitialized_;
   
   };
 
   /* 
   * Returns value of random seed (private member variable idum)
   */
   inline long CudaRandom::seed() 
   {  return seed_; }

} 
#endif 

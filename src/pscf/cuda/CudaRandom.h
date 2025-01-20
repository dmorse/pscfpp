#ifndef PSCF_CUDA_RANDOM_H
#define PSCF_CUDA_RANDOM_H

#include "GpuTypes.h"
#include "DeviceArray.h"

#include <curand.h>

/*
* PSCF Package - Polymer Self-Consistent Field 
*
* Copyright 2016 - 2023, The Regents of the University of Minnesota
* Distributed under the terms of the GNU General Public License.
*/

namespace Pscf {

   /**
   * Random number generator on GPU.
   *
   * The generator may be seeded either by reading a seed from 
   * file, using the readParam() method, or by using setSeed()
   * to set or reset it explicitly. In either case, inputting 
   * a positive integer causes that value to be used as a seed,
   * but inputting a value of 0 causes the use of a seed that 
   * is generated from the system clock. 
   *
   * \ingroup Pscf_Cuda_Module
   */
   class CudaRandom 
   {
 
   public:

      /**
      * Constructor.
      */   
      CudaRandom();

      /**
      * Destructor.
      */   
      virtual ~CudaRandom();
   
      /**
      * Set value of random seed (private member variable seed_).
      *
      * \param seed  value of random number generator seed.
      */
      void setSeed(unsigned long long seed);

      /**
      * Populate array with uniform random floating point numbers.
      *   
      * The array, which is stored on the device, is filled with random 
      * numbers in the range 0 < x <= 1 (zero excluded and 1 is included)
      *  
      * \param data  array to populate
      */
      void uniform(DeviceArray<cudaReal>& data);
   
      /**
      * Populate array with normal-distributed random floating point numbers.
      *   
      * The array, which is stored on the device, is filled with random 
      * numbers chosen from a normal distribution with the specified mean 
      * value and standard deviation.
      * 
      * Note: the input array must have an even number of elements. This is a 
      * requirement imposed by cuRAND, the random number generator software 
      * used by CudaRandom.
      *
      * \param data  array to populate
      * \param stddev  standard deviation (input)
      * \param mean  mean value (input, default = 0.0)
      */
      void normal(DeviceArray<cudaReal>& data, 
                  cudaReal stddev, cudaReal mean = 0.0);
   
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

      void errorCheck(curandStatus_t const & error);
   
   };
 
   /* 
   * Returns value of random seed (private member variable idum)
   */
   inline long CudaRandom::seed() 
   {  return seed_; }

} 
#endif 

#ifndef PSCF_CUDA_RANDOM_H
#define PSCF_CUDA_RANDOM_H

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
      * Populate array on device with random floats in (0, 1], uniform dist.
      *  
      * \param data  array to populate
      */
      void uniform(DeviceArray<float>& data);

      /**
      * Populate array on device with random doubles in (0, 1], uniform dist.
      *  
      * \param data  array to populate
      */
      void uniform(DeviceArray<double>& data);
   
      /**
      * Populate array on device with normal-distributed random floats.
      * 
      * Note: the input array must have an even number of elements. This is a 
      * requirement imposed by cuRAND, the random number generator software 
      * used by CudaRandom.
      *
      * \param data  array to populate
      * \param stddev  standard deviation (input)
      * \param mean  mean value (input, default = 0.0)
      */
      void normal(DeviceArray<float>& data, float stddev, float mean = 0.0);

      /**
      * Populate array on device with normal-distributed random doubles.
      * 
      * Note: the input array must have an even number of elements. This is a 
      * requirement imposed by cuRAND, the random number generator software 
      * used by CudaRandom.
      *
      * \param data  array to populate
      * \param stddev  standard deviation (input)
      * \param mean  mean value (input, default = 0.0)
      */
      void normal(DeviceArray<double>& data, double stddev, double mean = 0.0);
   
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

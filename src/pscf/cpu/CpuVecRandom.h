#ifndef PSCF_CPU_VEC_RANDOM_H
#define PSCF_CPU_VEC_RANDOM_H

/*
* PSCF - Polymer Self-Consistent Field
*
* Copyright 2015 - 2026, The Regents of the University of Minnesota
* Distributed under the terms of the GNU General Public License.
*/

// Forward declarations
namespace Util {
   class Random;
   template <typename T> class Array;
}

namespace Pscf {

   using namespace Util;

   /**
   * Random number generator for arrays of random numbers on CPU.
   *
   * A CpuVecRandom generates arrays of random real numbers on a CPU. 
   * It uses an associated Util::Random scalar random number generator
   * to generate these random numbers. 
   *
   * CpuVecRandom has an interface that is analogous to that of the
   * CudaVecRandom class, which uses a GPU to generates arrays of random 
   * numbers in global GPU memory. This similarity is intended to allow 
   * creation of template code that can use either type of vector 
   * random number object interchangably for use in analogous programs 
   * that run on CPU or GPU-accelerated hardware.
   *
   * \ingroup Pscf_Cpu_Module
   */
   class CpuVecRandom
   {

   public:

      /**
      * Default constructor.
      */
      CpuVecRandom();

      /**
      * Constructor - creates association with a scalar RNG.
      *
      * \param random  associated scalar random number generator (RNG)
      */
      CpuVecRandom(Util::Random& random);

      /**
      * Destructor.
      */
      virtual ~CpuVecRandom();

      /**
      * Create an association with a Util::Random scalar RNG.
      *
      * \param random  associated scalar random number generator
      */
      void associate(Util::Random& random);

      /**
      * Generate uniform random doubles in (0, 1].
      *
      * \param data  array to populate
      */
      void uniform(Array<double>& data);

      /**
      * Generate uniform random distribution in range (min, max].
      *
      * \param data  array to populate with random numbers
      * \param min  minimum of range
      * \param max  maximum of range
      */
      void uniform(Array<double>& data, double min, double max);

      /**
      * Generate normal-distributed random doubles.
      *
      * \param data  array to populate
      * \param stddev  standard deviation (input)
      * \param mean  mean value (input, default = 0.0)
      */
      void normal(Array<double>& data, double stddev, double mean = 0.0);

   private:

      /// Pointer to associated scalar random number generator.
      Util::Random* randomPtr_;

   };

   #if 0
   /*
   * Returns value of random seed (private member variable idum)
   */
   inline long CpuVecRandom::seed()
   {  return randomPtr_->seed(); }
   #endif

}
#endif

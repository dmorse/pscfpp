#ifndef PSCF_VEC_RANDOM_CP_H
#define PSCF_VEC_RANDOM_CP_H

/*
* PSCF - Polymer Self-Consistent Field
*
* Copyright 2015 - 2026, The Regents of the University of Minnesota
* Distributed under the terms of the GNU General Public License.
*/

#include <pscf/backend/cpp/CPT.h>

// Forward declarations
namespace Util {
   class Random;
   template <typename T> class Array;
}

namespace Pscf {

   using namespace Util;

   // Primary template declaration
   template <typename T> class VecRandom;

   /**
   * Random number generator for arrays of random numbers on CPU.
   *
   * A VecRandom<CPT> generates arrays of random real numbers on a CPU. 
   * It uses an associated Util::Random scalar random number generator
   * to generate these random numbers. 
   *
   * \ingroup Pscf_Backend_Cpp_Module
   */
   template <>
   class VecRandom<CPT>
   {

   public:

      /**
      * Default constructor.
      */
      VecRandom();

      /**
      * Constructor - creates association with a scalar RNG.
      *
      * \param random  associated scalar random number generator (RNG)
      */
      VecRandom(Util::Random& random);

      /**
      * Destructor.
      */
      ~VecRandom();

      // Prohibit copying and assignment
      VecRandom(VecRandom<CPT> const &) = delete;
      VecRandom<CPT> operator = (VecRandom<CPT> const &) = delete;

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
      void uniform(Util::Array<double>& data);

      /**
      * Generate uniform random distribution in range (min, max].
      *
      * \param data  array to populate with random numbers
      * \param min  minimum of range
      * \param max  maximum of range
      */
      void uniform(Util::Array<double>& data, double min, double max);

      /**
      * Generate normal-distributed random doubles.
      *
      * \param data  array to populate
      * \param stddev  standard deviation (input)
      * \param mean  mean value (input, default = 0.0)
      */
      void normal(Util::Array<double>& data, double stddev, double mean = 0.0);

   private:

      /// Pointer to associated scalar random number generator (non-owning).
      Util::Random* randomPtr_;

   };

}
#endif

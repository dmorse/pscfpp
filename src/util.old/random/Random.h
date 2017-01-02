#ifndef UTIL_RANDOM_H
#define UTIL_RANDOM_H

/*
* Simpatico - Simulation Package for Polymeric and Molecular Liquids
*
* Copyright 2010 - 2014, The Regents of the University of Minnesota
* Distributed under the terms of the GNU General Public License.
*/

#include <util/param/ParamComposite.h>
#include <util/space/Vector.h>

#include <cmath>

#include <util/random/mersenne/mtrand.h>
#define UTIL_ENGINE MTRand_int32

namespace Util
{

   class Vector;

   /**
   * Random number generator.
   *
   * This class provides functions that return several forms of
   * random numbers, using an internal Mersenne-Twister random 
   * number generator.  
   *
   * The generator may be seeded either by reading a seed from 
   * file, using the readParam() method, or by using setSeed()
   * to set or reset it explicitly. In either case, inputting 
   * a positive integer causes that value to be used as a seed,
   * but inputting a value of 0 causes the use of a seed that 
   * is generated from the system clock. 
   *
   * If the program is compiled with MPI, and MPI is initialized, 
   * then any automatically generated seed is also offset by a 
   * value that depends on the rank of the processor within the 
   * MPI world communicator, so that different processor use
   * different seeds. 
   *
   * \ingroup Random_Module
   */
   class Random : public ParamComposite
   {
 
   public:

      typedef unsigned long SeedType;

      /**
      * Constructor.
      */   
      Random();

      /**
      * Destructor.
      */   
      virtual ~Random();

      /**
      * Read seed from file, initialize RNG.
      *
      * \param in input stream.
      */
      virtual void readParameters(std::istream &in);
   
      /**
      * Load internal state from file.
      *
      * \param ar input/loading archive
      */
      virtual void loadParameters(Serializable::IArchive &ar);
   
      /**
      * Save internal state to file. 
      *
      * \param ar output/saving archive
      */
      virtual void save(Serializable::OArchive &ar);
   
      /**
      * Sets of random seed, and initializes random number generator.
      *
      * \param seed value for random seed (private member variable idum)
      */
      void setSeed(SeedType seed);
   
      /**
      * Return a random floating point number x, uniformly distributed in the
      * range 0 <= x < 1.
      *
      * \return random double precision number
      */
      double uniform();
   
      /**
      * Return a random floating point number x, uniformly distributed in the
      * range range1 <= x < range2.
      *
      * \return random double precision number
      */
      double uniform(double range1, double range2);
   
      /**
      * Return random long int x uniformly distributed in range1 <= x < range2.
      *
      * Parameters range1 and range2 must be within the range of long integers.
      *
      * \return random integer
      */
      long uniformInt(long range1, long range2);
    
      /**
      * Generate a random point in a box.
      *
      * \param  minR[] array of minimum coordinate values along three axes
      * \param  maxR[] array of maximum coordinate values along three axes
      * \param  r[]    random position such that minR[axis] < r[axis] < maxR[axis]
      */
      void getPoint(double minR[], double maxR[], double r[]);
   
      /**
      * Return a Gaussian random number with zero average and unit variance.
      *
      * \return Gaussian distributed random number.
      */
      double gaussian(void);
   
      /**
      * Generate unit vector with uniform probability over the unit sphere
      *
      * \param v random unit vector (upon return)
      */
      void unitVector(Vector& v);
   
      /**
      * Metropolis algorithm for whether to accept a MC move.
      *
      * If ratio > 1, this function return true. If 0 < ratio < 1,
      * this function returns true with probability ratio, and false
      * with probability 1 - ratio.
      *
      * \param  ratio ratio of old to new equilibrium weights
      * \return true if accepted, false if rejected
      */
      bool metropolis(double ratio);

      /** 
      * Choose one of several outcomes with a specified set of probabilities.
      *
      * Precondition: Elements of probability array must add to 1.0
      *
      * \param probability[] array of probabilities, for indices 0,...,size-1
      * \param size          number of options
      *
      * \return random integer index of element of probability[] array
      */
      long drawFrom(double probability[], long size);
    
      /**
      * Serialize to/from an archive.
      */
      template <class Archive>
      void serialize(Archive& ar, const unsigned int version);

      /**
      * Returns value of random seed (private member variable seed_).
      *
      * \return value of random number generator seed.
      */
      long seed();
   
   private:

      /**
      * Uniform random number generator engine.
      */
      UTIL_ENGINE engine_;

      /**
      * Seed value.
      */
      SeedType seed_;

      /// Has a seed been set by readParam() or setSeed()?
      bool isInitialized_;
   
      /**
      * Initialize random number generator.
      */
      void setSeed();
   
   };
  
   // Inline methods 

   /*
   * Get a uniform double precision number.
   */
   inline double Random::uniform()
   {
      // Divide 32 bit integer engine output by 2^32
      return static_cast<double>(engine_()) * (1.0/ 4294967296.0); 
   }
 
   /*
   * Get a uniform double precision number.
   */
   inline double Random::uniform(double range1, double range2)
   {
      // Divide 32 bit integer engine output by 2^32
      double frac = static_cast<double>(engine_()) * (1.0/ 4294967296.0); 
      return range1 + (range2 - range1)*frac; 
   }
 
   /* 
   * Return a random long integer x uniformly distributed in range1 <= x < range2.
   *
   * Parameters range1 and range2 must be within the range of long integers.
   */
   inline long Random::uniformInt(long range1, long range2) {
      double x = uniform(range1, range2);
      return long(x);
   }
   
   /* 
   * Choose one of several outcomes with a specified set of probabilities.
   *
   * Precondition: Elements of probability array must add to 1.0
   *
   * \param probability[] array of probabilities, for indices 0,...,size-1
   * \param size          number of options
   */
   inline long Random::drawFrom(double probability[], long size) {
      double roulette = uniform();
      long   n=0;
      double cumProb = probability[0];
      while ( cumProb < roulette && n < size ) {
         ++n;
         cumProb += probability[n];
      }
      return(n);
   }
   
   /*
   * Get random point within a rectangular prism.
   */
   inline void Random::getPoint(double minR[], double maxR[], double r[]) {
      r[0] = uniform(minR[0], maxR[0]);
      r[1] = uniform(minR[1], maxR[1]);
      r[2] = uniform(minR[2], maxR[2]);
   }
   
   /* 
   * Metropolis criteria for whether to accept a MC move.
   */
   inline bool Random::metropolis(double ratio)
   {
      if ( ratio > 1.0 ) {
         return true;
      } else {
         double ran = uniform();
         if ( ran < ratio ) {
            return true;
         } else {
            return false;
         }
      }
   }

   /* 
   * Returns value of random seed (private member variable idum)
   */
   inline long Random::seed() 
   {  return seed_; }

   /*
   * Serialize to/from an archive.
   */
   template <class Archive>
   void Random::serialize(Archive& ar, const unsigned int version)
   {
      ar & engine_;
      ar & seed_;
   }

} 
#undef UTIL_ENGINE 
#endif // ifndef UTIL_RANDOM_H

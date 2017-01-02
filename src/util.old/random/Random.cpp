#include "Random.h"

#include <cmath>
#include <cstdio>
#include <sys/time.h>

namespace Util
{

   /*
   * Constructor.
   */
   Random::Random()
    : engine_(),
      seed_(0),
      isInitialized_(false)
   {  setClassName("Random"); }

   /*
   * Destructor.
   */
   Random::~Random()
   {}

   /*
   * Read random seed and initialize.
   */
   void Random::readParameters(std::istream &in)
   {
      read<SeedType>(in, "seed", seed_);
      setSeed();
   }

   /*
   * Load random seed and internal state.
   */
   void Random::loadParameters(Serializable::IArchive& ar)
   {
      loadParameter<SeedType>(ar, "seed", seed_);
      ar >> engine_;
   }

   /*
   * Save internal state.
   */
   void Random::save(Serializable::OArchive& ar)
   {
      ar << seed_;
      ar << engine_;
   }

   /*
   * Sets of random seed, and initializes random number generator.
   *
   * \param seed value for random seed (private member variable seed)
   */
   void Random::setSeed(Random::SeedType seed)
   {
      seed_ = seed;
      setSeed();
   }

   /*
   * If seed_ == 0, use a seed taken from the clock. Otherwise, use
   * seed_ as an integer seed for the random number generator.
   */
   void Random::setSeed()
   {
      // If seed  == 0, replace with random seed from clock.
      if (seed_ == 0) {
         SeedType temp;
         //temp = static_cast<SeedType>(std::time(0)); 
         timeval time;
         gettimeofday(&time, NULL);
         temp = time.tv_sec + 1123*time.tv_usec;
         #ifdef UTIL_MPI
         if (MPI::Is_initialized()) {
            SeedType rank = MPI::COMM_WORLD.Get_rank();
            temp += rank*(31 + time.tv_usec);
         }
         #endif
         engine_.seed(temp);
      } else {
         engine_.seed(seed_);
      }
      isInitialized_ = true;
   }

   /*
   * Return Gaussian distributed random number.
   */
   double Random::gaussian(void)
   {
      static bool   iset = false;
      static double zGauss1;
      double zGauss2;

      if (iset == false) {
         double v1;
         double v2;
         double rsq=2.0;

         while( rsq >= 1.0) {
            v1 = 2.0*uniform(-1.0, 1.0);
            v2 = 2.0*uniform(-1.0, 1.0);
            rsq = v1*v1 + v2*v2;
         }
         double fac = sqrt(-2.0*log(rsq)/rsq);
         zGauss1 = v1*fac;
         zGauss2 = v2*fac;
         iset = true;
      } else {
         zGauss2 = zGauss1;
         iset = false;
      }
      return zGauss2;
   }

   /*
   * unit vector with uniform probability over the unit sphere
   */
   void Random::unitVector(Vector& v)
   {
      double ran1=0, ran2=0, ranh;
      double ransq;

      ransq=2.0;
      while (ransq >= 1.0) {
         ran1  = 1.0-2.0*uniform(0.0, 1.0);
         ran2  = 1.0-2.0*uniform(0.0, 1.0);
         ransq = ran1*ran1 + ran2*ran2;
      }
      ranh= 2.0*sqrt(1.0-ransq);
      v[0] = ran1*ranh;
      v[1] = ran2*ranh;
      v[2] = 1.0 - 2.0*ransq;
   }


}

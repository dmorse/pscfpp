#ifndef UTIL_AR1_PROCESS_H
#define UTIL_AR1_PROCESS_H

/*
* Simpatico - Simulation Package for Polymeric and Molecular Liquids
*
* Copyright 2010 - 2014, The Regents of the University of Minnesota
* Distributed under the terms of the GNU General Public License.
*/

#include <util/random/Random.h>
#include <cmath>

namespace Util
{

   /**
   * Generator for a discrete AR(1) Markov process.
   *
   * An auto-regressive AR(1) process is a discrete stationary Markov process 
   * x(n) with an autocorrelation function <x(n)*x(n+m)> = exp(-m/tau), where
   * tau is a decay time. It is a discrete version of the Ornstein-Uhlenbeck
   * continuous Markov process.
   *
   * \ingroup Random_Module
   */
   class Ar1Process 
   {
   
   public:

      /**
      * Constructor.
      */   
      Ar1Process();

      /**
      * Constructor.
      *
      * \param random associated random number generator.
      */   
      Ar1Process(Random& random);

      /**
      * Associate a random number generator.
      *
      * \param random associated random number generator.
      */
      void setRNG(Random& random);

      /**
      * Initialize process
      *
      * \param tau decay time (in discrete steps)
      */
      void init(double tau);

      /**
      * Generate and return a new value.
      */
      double operator () ();

   private:

      double x_;

      double B_;

      double C_;

      Random* randomPtr_;

      bool isInitialized_;

   };

   inline double Ar1Process::operator() ()
   {
      assert(isInitialized_);
      x_ *= C_;
      x_ += B_*randomPtr_->gaussian();
      return x_;
   }

} 
#endif

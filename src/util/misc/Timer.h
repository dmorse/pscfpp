#ifndef UTIL_TIMER_H
#define UTIL_TIMER_H

/*
* Simpatico - Simulation Package for Polymeric and Molecular Liquids
*
* Copyright 2010 - 2014, The Regents of the University of Minnesota
* Distributed under the terms of the GNU General Public License.
*/

#include <util/global.h>

#include <ctime>

namespace Util
{

   /**
   * Wall clock timer.
   *
   * A Timer keeps track of the time elapsed during one or more
   * interval. Each interval begins when start() is called and
   * ends when stop() is called. If start() and stop() are invoked 
   * repeatedly, the timer accumulates the time elapses in multiple 
   * intervals.  The accumulated time is returned by the time()
   * method, and can be reset to zero by the clear() method.
   *
   * \ingroup Misc_Module
   */
   class Timer
   {
   
   public:
	
      /// Constructor.
      Timer();

      /// Start the clock.
      void start();

      /// Stop the clock, increment time.
      void stop();
      
      /// Set accumulated time to zero.
      void clear();

      /// Get the accumulated time, in seconds
      double time();

      /// Is the timer running?
      bool isRunning();

   private:

      /// Accumulated time.
      double time_;
 
      /// Beginning of interval when Timer is running
      double begin_;
 
      /// Is the timer running now? 
      bool   isRunning_; 

   };

   /**
   * Constructor
   */
   inline Timer::Timer()
    : time_(0.0),
      begin_(0.0),
      isRunning_(false)
   {}
  
   /*
   * Start the timer.
   */ 
   inline void Timer::start()
   {
      if (isRunning_) 
         UTIL_THROW("Attempt to restart an active Timer");
      isRunning_ = true;
      begin_ = clock(); 
   }
   
   /*
   * Stop the timer.
   */ 
   inline void Timer::stop()
   {
      double end = clock();
      if (!isRunning_) 
         UTIL_THROW("Attempt to stop an inactive Timer");
      isRunning_ = false;
      time_ += double(end - begin_)/double(CLOCKS_PER_SEC); 
   }
   
   /*
   * Clear the timer.
   */ 
   inline void Timer::clear()
   {
      time_      = 0.0;
      isRunning_ = false;
   }

   /*
   * Get the current time.
   */ 
   inline double Timer::time()
   {  return time_; }

   /*
   * Is this timer running?
   */ 
   inline bool Timer::isRunning()
   {  return isRunning_; }

}
#endif

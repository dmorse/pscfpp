/*
* Simpatico - Simulation Package for Polymeric and Molecular Liquids
*
* Copyright 2010 - 2014, The Regents of the University of Minnesota
* Distributed under the terms of the GNU General Public License.
*/

#include "Signal.h"

namespace Util
{

   /*
   * Destructor.
   */
   Signal<>::~Signal()
   {  clear(); }

   /* 
   * Notify observers (call associated methods).
   */
   void Signal<>::notify()
   {
      std::list< IFunctor<>* >::iterator pos;
      pos = functorPtrs_.begin();
      while (pos != functorPtrs_.end())
      {
         (**pos)();
         ++pos;
      }
   }

   /* 
   * Notify observers (call associated methods).
   */
   void Signal<>::clear()
   {
      std::list< IFunctor<>* >::iterator pos;
      pos = functorPtrs_.begin();
      while (pos != functorPtrs_.end())
      {
         delete *pos;
         ++pos;
      }
      functorPtrs_.clear();
   }

   /* 
   * Get number of registered observers.
   */
   int Signal<>::nObserver() const
   { return functorPtrs_.size(); }

}

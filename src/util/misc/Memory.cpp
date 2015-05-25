/*
* Simpatico - Simulation Package for Polymeric and Molecular Liquids
*
* Copyright 2010 - 2014, The Regents of the University of Minnesota
* Distributed under the terms of the GNU General Public License.
*/

#include "Memory.h"

namespace Util
{

   /// Number of calls to allocate.
   int Memory::nAllocate_ = 0;

   /// Number of calls to de-allocate.
   int Memory::nDeallocate_ = 0;

   /// Total amount of memory allocated, in bytes.
   int Memory::total_ = 0;

   /// Maximum of total over course of simulation.
   int Memory::max_ = 0;

   /*
   * Call this to ensure compilation of this file. 
   */
   void Memory::initStatic()
   {  max_ = 0; }  

   /*
   * Return number of calls to allocate.
   */
   int Memory::nAllocate()
   {  return nAllocate_; }

   /*
   * Return number of calls to deallocate.
   */
   int Memory::nDeallocate()
   {  return nDeallocate_; }

   /*
   * Return total amount of memory allocated thus far.
   */
   int Memory::total()
   {  return total_; }

   /*
   * Return maximum amount of allocated memory thus far.
   */
   int Memory::max()
   {  return max_; }

   #ifdef UTIL_MPI
   int Memory::max(MPI::Intracomm& communicator)
   { 
      int maxGlobal;
      int maxLocal = max_;
      communicator.Allreduce(&maxLocal, &maxGlobal, 1, MPI::INT, MPI::MAX);
      return maxGlobal;
   }
   #endif

} 

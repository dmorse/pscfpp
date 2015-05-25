#ifndef UTIL_MEMORY_H
#define UTIL_MEMORY_H

/*
* Simpatico - Simulation Package for Polymeric and Molecular Liquids
*
* Copyright 2010 - 2014, The Regents of the University of Minnesota
* Distributed under the terms of the GNU General Public License.
*/

#include <util/global.h>
#include <stddef.h>
#include <iostream>
#include <new>

namespace Util
{

   /**
   * Provides method to allocate array.
   *
   * The Memory::allocate() method invokes the new operator within
   * a try catch block, and keeps track of the total memory 
   * allocated.
   *
   * \ingroup Misc_Module
   */
   class Memory
   { 
   public:

      /**
      * Allocate a C array.
      *
      * Allocates a Data array of size elements, assigns ptr the
      * address of the first element. 
      */
      template <typename Data>
      static void allocate(Data*& ptr, size_t size);

      /**
      * Allocate a C array.
      *
      * Allocates a Data array of size elements, assigns ptr the
      * address of the first element. 
      */
      template <typename Data>
      static void deallocate(Data*& ptr, size_t size);

      /**
      * Return number of times allocated was called.
      */
      static int nAllocate();

      /**
      * Return number of times deallocate was called.
      */
      static int nDeallocate();

      /**
      * Return total amount of memory currently allocated.
      */
      static int total();

      /**
      * Return the maximum amount of allocated heap memory thus far.
      *
      * This function returns the temporal maximum of total().
      */
      static int max();

      #ifdef UTIL_MPI
      /**
      * Return max for any processor in communicator.
      */
      static int max(MPI::Intracomm& communicator);
      #endif

      /**
      * Call this just to guarantee initialization of static memory.
      */
      static void initStatic();
   
   private: 

      /// Total amount of memory allocated, in bytes. 
      static int total_;
   
      /// Maximum amount of memory allocated, in bytes. 
      static int max_;
   
      /// Number of calls to allocate.
      static int nAllocate_;
   
      /// Number of calls to deallocate.
      static int nDeallocate_;
   
   };
   
   /*
   * Allocate a C array.
   */
   template <typename Data>
   void Memory::allocate(Data*& ptr, size_t size)
   {
     if (ptr) {
         UTIL_THROW("Attempt to allocate to non-null pointer");
      }
      try {
         ptr = new Data[size];
         total_ += (size*sizeof(Data));
         ++nAllocate_;
         if (total_ > max_) max_ = total_;
      } catch (std::bad_alloc&) {
         std::cout << "Allocation error" << std::endl;
         throw;
      }
   }

   /*
   * De-allocate a C array.
   */
   template <typename Data>
   void Memory::deallocate(Data*& ptr, size_t size)
   {
      if (ptr) {
         delete [] ptr;
         ptr = 0;
         total_ -= size*sizeof(Data);
         ++nDeallocate_;
      } else {
         UTIL_THROW("Attempt to de-allocate null pointer");
      }
   }

} 
#endif

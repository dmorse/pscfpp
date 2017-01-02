#ifndef UTIL_SCOPED_PTR_H
#define UTIL_SCOPED_PTR_H

/*
* Simpatico - Simulation Package for Polymeric and Molecular Liquids
*
* Copyright 2010 - 2014, The Regents of the University of Minnesota
* Distributed under the terms of the GNU General Public License.
*/

#include <util/global.h>

namespace Util
{

   /**
   * A very simple RAII pointer.
   *
   * A ScopedPtr mimics a built-in pointer, except that guarantees destruction 
   * of the object to which it points when the ScopedPtr goes out of scope. 
   * It accepts ownership of a built-in pointer either upon construction or by 
   * the reset() method, and deletes the associated object in its destructor.
   * A ScopedPtr cannot be copy constructed or assigned.
   *
   * Similar to boost::scoped_ptr, with minor differences. It takes the same
   * amount of memory as a built-in pointer, and should be equally fast.
   */
   template <typename T>
   class ScopedPtr
   {

   public:

      /// Type of object pointed to.
      typedef T element_type;

      /// Constructor.
      explicit ScopedPtr(T* p = 0)
      {  ptr_ = p; }

      /// Destructor, destroys object pointed to, if any.
      ~ScopedPtr()
      {
         if (ptr_ != 0) {  
            delete ptr_; 
         }
      }

      /**
      * Acquire ownership of a built-in pointer.
      * 
      * \param p built-in pointer to acquire.
      */
      void reset(T* p = 0)
      {
         if (ptr_ != 0) {
            delete ptr_;
         }  
         ptr_ = p;
      }

      /// Dereference.
      T& operator*() const
      {  return *ptr_; }

      /// Member access.
      T* operator->() const
      {  return ptr_; }

      /// Return enclosed built-in pointer.
      T* get() const
      {  return ptr_; }

   private:

      T* ptr_;

      /// Copy constructor - private and not implemented.
      ScopedPtr(const ScopedPtr&);

      /// Assignment - private and not implemented.
      ScopedPtr& operator = (const ScopedPtr& );

   };

   /** 
   * Return true iff the enclosed built-in pointer is null.
   */
   template <typename T>
   inline bool isNull(ScopedPtr<T> p)
   {  return (p.get() == 0); }

}
#endif 

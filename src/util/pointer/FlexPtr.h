#ifndef UTIL_FLEX_PTR_H
#define UTIL_FLEX_PTR_H

/*
* Simpatico - Simulation Package for Polymeric and Molecular Liquids
*
* Copyright 2010 - 2014, The Regents of the University of Minnesota
* Distributed under the terms of the GNU General Public License.
*/

#include "isNull.h"
#include <util/global.h>

namespace Util
{

   /**
   * A pointer that may or may not own the object to which it points.
   *
   * FlexPtr overloads * and ->, and thus mimics a built-in pointer in 
   * most respects.
   *
   * The acquire(T*) method copies a built-in pointer and accept ownership
   * of the object to which it points, i.e., accepts responsibility for 
   * deleting the object, normally when the FlexPtr goes out of scope. 
   *
   * The copy(T*) method copies a built-in pointer without accepting 
   * ownership, i.e., without accepting responsibility for deleting the
   * object to which it points.
   * 
   * Both acquire() and copy() destroy any object that is already owned
   * by this FlexPtr before copying of a new pointer.
   */
   template <typename T>
   class FlexPtr
   {

   public:

      /// Type of object pointed to.
      typedef T element_type;

      /**
      * Constructor.
      */
      FlexPtr()
       : ptr_(0),
         isOwner_(0)
      {}

      /**
      * Destructor.
      * 
      * Deletes any object that is owned by this FlexPtr.
      */
      ~FlexPtr()
      {
         if (ptr_ != 0 && isOwner_) {  
            delete ptr_; 
         }
      }

      /**
      * Copy a built-in pointer, and accept ownership.
      *
      * If this FlexPtr already owns an object, it will be deleted before
      * acquiring a new pointer.
      *
      * Throws an Exception if p is null.
      *
      * \param p Built-in pointer to be acquired.
      */
      void acquire(T* p)
      {
         if (p == 0) {
            UTIL_THROW("Cannot acquire a null pointer");
         }  
         if (ptr_ != 0 && isOwner_) {
            delete ptr_; 
         }  
         ptr_     = p;
         isOwner_ = 1;
      }

      /**
      * Copy a built-in pointer, without accepting ownership.
      * 
      * If this FlexPtr already owns an object, it will be deleted
      * before copying a new pointer.
      *
      * \param p Built-in pointer to be copied.
      */
      void copy(T* p)
      {
         if (p == 0) {
            UTIL_THROW("Cannot copy a null pointer");
         }  
         if (ptr_ != 0 && isOwner_) {
            delete ptr_; 
         }  
         ptr_     = p;
         isOwner_ = 0;
      }

      /**
      * Dereference.
      */
      T& operator*() const
      {  return *ptr_; }

      /**
      * Member access.
      */
      T* operator->() const
      {  return ptr_; }

      /**
      * Return the built-in pointer.
      */
      T* get() const
      {  return ptr_; }

   private:

      /// Built-in pointer.
      T*  ptr_;

      /// True iff ptr_ should be deleted in the destructor.
      int isOwner_;

      /// Copy constructor - private and not implemented.
      FlexPtr(const FlexPtr&);

      /// Assignment - private and not implemented.
      FlexPtr& operator = (const FlexPtr& );

   };

   /** 
   * Return true iff the enclosed built-in pointer is null.
   */
   template <typename T>
   inline bool isNull(FlexPtr<T> p)
   {  return (p.get() == 0); }

}

#endif

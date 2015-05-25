#ifndef UTIL_ARRAY_STACK_H
#define UTIL_ARRAY_STACK_H

/*
* Simpatico - Simulation Package for Polymeric and Molecular Liquids
*
* Copyright 2010 - 2014, The Regents of the University of Minnesota
* Distributed under the terms of the GNU General Public License.
*/

#include <util/misc/Memory.h>
#include <util/global.h>

namespace Util
{

   /**
   * A stack of fixed capacity.
   *
   * Pointers to elements are stored in an allocatable, non-resizable array.
   *
   * \ingroup Pointer_Array_Module
   */
   template <typename Data>
   class ArrayStack
   { 

   public: 
   
      /**
      * Default constructor.
      */
      ArrayStack();
   
      /**
      * Destructor.
      */
      virtual ~ArrayStack();
    
      /**
      * Initialize and allocate required memory.
      *
      * \param capacity maximum size of stack.
      */
      void allocate(int capacity);
 
      /// \name Mutators
      //@{

      /**
      * Push an element onto the Stack.
      * 
      * \param data element to be added to stack.
      */
      void push(Data& data);

      /**
      * Pop an element off the stack.
      *
      * \return the top element (which is popped off stack).
      */
      Data& pop();

      //@}
      /// \name Accessors
      //@{
 
      /**
      * Return capacity of the underlying array.
      *
      * \return Number of elements allocated in array.
      */
      int capacity() const;

      /**
      * Get the number of elements in the stack.
      */
      int size() const; 

      /**
      * Return a reference to the top element (don't pop).
      */
      Data& peek();

      /**
      * Return a const ref to the top element (don't pop).
      */
      const Data& peek() const;

      /**
      * Return true if the ArrayStack is valid, or throw an exception.
      */
      bool isValid() const;
   
      /**
      * Return true only if the ArrayStack has been allocated.
      */
      bool isAllocated() const;
   
      //@}

   private:

      // C array of pointers to elements of data_ array
      Data **ptrs_;
  
      // Allocated size of ptrs_ array (maximum size of stack).
      int  capacity_;
  
      // Size of stack.
      int  size_; 

      /// Copy constructor, declared private to prohibit copying.
      ArrayStack(const ArrayStack&);

      /// Assignment operator, declared private to prohibit copying.
      ArrayStack& operator = (const ArrayStack&);

   }; 


   /* 
   * Default constructor.
   */
   template <typename Data>
   ArrayStack<Data>::ArrayStack()
    : ptrs_(0),    
      capacity_(0),
      size_(0)
   {}

   /* 
   * Destructor.
   */
   template <typename Data>
   ArrayStack<Data>::~ArrayStack()
   {
      if (ptrs_) {
         Memory::deallocate<Data*>(ptrs_, capacity_);
      }
   }
 
   /* 
   * Allocate and initialize required memory.
   */
   template <typename Data>
   void ArrayStack<Data>::allocate(int capacity) 
   {
      if (capacity <=0) {
         UTIL_THROW("Cannot allocate ArrayStack with capacity <= 0");
      }
      if (ptrs_) {
         UTIL_THROW("Attempt to re-allocate an ArrayStack");
      }
      Memory::allocate<Data*>(ptrs_, capacity);
      capacity_ = capacity;
      size_ = 0;

      // Nullify all Data* pointer elements
      for (int i = 0; i < capacity_; ++i) {
         ptrs_[i] =  0;
      }
   }

   /*
   * Return capacity of the associated Array.
   */
   template <typename Data>
   int ArrayStack<Data>::capacity() const
   {  return capacity_; }

   /*
   * Get the current size.
   */
   template <typename Data>
   inline int ArrayStack<Data>::size() const
   {  return size_; }

   /*
   * Push a Data object onto the stack.
   */
   template <typename Data>
   void ArrayStack<Data>::push(Data& data)
   {
      if (size_ == capacity_) {
         UTIL_THROW("Attempt to push onto full stack");
      }
      ptrs_[size_] = &data;
      ++size_;
   }

   /*
   * Pop a Data object off the stack.
   */
   template <typename Data>
   Data& ArrayStack<Data>::pop()
   {
      if (size_ == 0) {
         UTIL_THROW("Attempt to pop from empty stack");
      }
      Data *ptr = ptrs_[size_-1];
      ptrs_[size_-1] = 0;
      --size_;
      return *ptr;
   }

   /* 
   * Return a reference to the top element, without popping.
   */
   template <typename Data>
   inline Data& ArrayStack<Data>::peek()
   {  return *ptrs_[size_-1]; }

   /* 
   * Return a const reference to the top element, without popping.
   */
   template <typename Data>
   inline const Data& ArrayStack<Data>::peek() const
   {  return *ptrs_[size_-1]; }

   /* 
   * Return true if valid, or throw an exception.
   */
   template <typename Data>
   bool ArrayStack<Data>::isValid() const
   {
      if (capacity_ < 0) {
         UTIL_THROW("Negative capacity_");
      }
      if (size_ < 0) {
         UTIL_THROW("Negative size_");
      }
      if (size_ > capacity_) {
         UTIL_THROW("size_ > capacity_");
      }

      if (ptrs_ != 0) {
         int i;
         for (i = 0; i < size_ ; ++i) {
            if (ptrs_[i] == 0) {
               UTIL_THROW("Null ptrs_[i] for i < size_");
            }
         }
         for (i = size_; i < capacity_ ; ++i) {
            if (ptrs_[i] != 0) {
               UTIL_THROW("Non-null ptrs_[i] for i >= size_");
            }
         }
      } else {
         if (capacity_ != 0) {
            UTIL_THROW("Unallocated stack but capacity != 0");
         }
      }
      return true;
   }

   /* 
   * Return true if the ArrayStack has been allocated, false otherwise.
   */
   template <typename Data>
   inline bool ArrayStack<Data>::isAllocated() const
   {  return (ptrs_ != 0); }

} 
#endif

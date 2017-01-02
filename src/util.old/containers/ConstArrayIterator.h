#ifndef UTIL_CONST_ARRAY_ITERATOR_H
#define UTIL_CONST_ARRAY_ITERATOR_H

/*
* Simpatico - Simulation Package for Polymeric and Molecular Liquids
*
* Copyright 2010 - 2014, The Regents of the University of Minnesota
* Distributed under the terms of the GNU General Public License.
*/

namespace Util
{

   /**
   * Forward const iterator for an Array or a C array.
   *
   * A ConstArrayIterator is a forward iterator for an array of const data.
   * It is an abstraction of a pointer to const, similar to an STL const
   * forward iterator. The * operator returns a const reference to an associated 
   * Data object the -> operator returns a const pointer to that object. The 
   * ++ operator increments the current pointer by one array element.
   *
   * Unlike an STL const forward iterator, an ConstArrayIterator contains the 
   * address of the end of the array. The isEnd() method can be used to test 
   * for termination of a for or while loop. When isEnd() is true, the current
   * pointer is one past the end of the array, and thus the iterator has no
   * current value.
   *
   * A ConstArrayIterator behave like a pointer to constant data, and provides
   * read-only access to the object to which it points. Use an ArrayIterator
   * if you need read-write access. 
   *
   * \ingroup Array_Module
   * \ingroup Iterator_Module
   */
   template <typename Data>
   class ConstArrayIterator
   {

   public:

      /**
      * Default constructor.  
      *  
      * Constructs an uninitialized iterator.  
      */
      ConstArrayIterator()
       : current_(0),
         end_(0)
      {}

      /**
      * Set the current pointer value.
      *
      * \param ptr Pointer to current element of the array.
      */
      void setCurrent(Data *ptr)
      {  current_ = ptr; }

      /**
      * Set the value of the end pointer.
      *
      * \param ptr Pointer to one element past end of array.
      */
      void setEnd(Data* ptr)
      {  end_ = ptr; }

      /**
      * Has the end of the array been reached?
      *
      * \return true if at end, false otherwise.
      */
      bool isEnd() const
      {  return (current_ == end_); }

      /**
      * Is this not the end of the array?
      *
      * \return true if not at end, false otherwise.
      */
      bool notEnd() const
      {  return (current_ != end_); }

      /**
      * Return a pointer to the current data.
      *
      * \return true if at end, false otherwise.
      */
      const Data* get() const
      {  return current_; }

      /// \name Operators
      //@{
      
      /**
      * Get a reference to the current Data.
      *
      * \return reference to associated Data object
      */
      const Data& operator* () const
      {  return *current_; }

      /**
      * Provide a pointer to the current Data object.
      *
      * \return const pointer to the Data object
      */
      const Data* operator -> () const
      {  return current_; }

      /**
      * Increment the current pointer.
      *
      * \return this ConstArrayIterator, after modification.
      */
      ConstArrayIterator<Data>& operator++ ()
      {
         ++current_;
         return *this;
      }

      //@}
      
   private:

      /// Pointer to the current element.
      Data* current_;

      /// Pointer to element one past last in the array.
      Data* end_;

   };

} 
#endif

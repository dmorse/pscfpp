#ifndef UTIL_CONST_P_ARRAY_ITERATOR_H
#define UTIL_CONST_P_ARRAY_ITERATOR_H

/*
* Simpatico - Simulation Package for Polymeric and Molecular Liquids
*
* Copyright 2010 - 2014, The Regents of the University of Minnesota
* Distributed under the terms of the GNU General Public License.
*/

namespace Util
{

   /**
   * Forward iterator for a PArray.
   *
   * An ConstPArrayIterator is an abstraction of a pointer, similar to an STL
   * forward iterator. The * operator returns a reference to an associated 
   * Data object, the -> operator returns a pointer to that object. The ++ 
   * operator increments the current pointer by one array element.
   *
   * Unlike an STL forward iterator, an ConstPArrayIterator contains the address 
   * of the end of the array. The isEnd() method can be used to test for
   * termination of a for or while loop. When isEnd() is true, the iterator
   * has no current value, and cannot be incremented further. The isEnd()
   * method returns true either if the iterator: i) has already been 
   * incremented one past the end of an associated PArray, or ii) is in a
   * null state that is produced by the constructor and the clear() method. 
   *
   * \ingroup Pointer_Array_Module
   * \ingroup Iterator_Module
   */
   template <typename Data>
   class ConstPArrayIterator
   {

   public:

      /**
      * Default constructor.  
      *  
      * Constructs a null iterator.  
      */
      ConstPArrayIterator()
       : current_(0),
         end_(0),
         data_(0)
      {}

      /**
      * Set the current pointer value.
      *
      * \param ptr Pointer to current element of array of Data* pointers.
      */
      void setCurrent(Data** ptr)
      { 
         current_ =  ptr;
         data_    = *ptr; 
      }

      /**
      * Set the value of the end pointer.
      *
      * \param ptr Pointer to one element past end of array of Data* pointers.
      */
      void setEnd(Data** ptr)
      {  end_ = ptr; }

      /**
      * Nullify the iterator.
      */
      void setNull()
      {
         current_ = 0; 
         end_ = 0; 
         data_ = 0; 
      }

      /**
      * Is the current pointer at the end of the array.
      *
      * \return true if at end, false otherwise.
      */
      bool isEnd() const
      { return (current_ == end_); }

      /**
      * Is the current pointer not at the end of the array?
      *
      * \return true if not at end, false otherwise.
      */
      bool notEnd() const
      { return (current_ != end_); }

      /**
      * Return a pointer to const current data.
      *
      * \return true if at end, false otherwise.
      */
      const Data* get() const
      { return data_; }

      /// \name Operators
      //@{
      
      /**
      * Return a const refererence to the current Data.
      *
      * \return const reference to the Data object
      */
      const Data& operator* () const
      { return *data_; }

      /**
      * Provide a pointer to the current Data object.
      *
      * \return pointer to the Data object
      */
      const Data* operator -> () const
      {  return data_; }

      /**
      * Increment the current pointer.
      *
      * \return this ConstPArrayIterator, after modification.
      */
      ConstPArrayIterator<Data>& operator++ ()
      {
         ++current_;
         if (current_ != end_) {
           data_ = *current_;
         } else {
           data_ = 0;
         }
         return *this;
      }

      //@}
      
   private:

      // Pointer to the current Data* pointer.
      Data** current_;

      // Pointer to one element one past last Data* pointer in the array.
      Data** end_;

      // Pointer to current Data object.
      Data*  data_;

   };

} 
#endif

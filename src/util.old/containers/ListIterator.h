#ifndef UTIL_LIST_ITERATOR_H
#define UTIL_LIST_ITERATOR_H

#include <util/containers/Node.h>

namespace Util
{

   template <typename Data> class List;

   /**
   * Bidirectional iterator for a List.
   *
   * A ListIterator provides bidirectional input/output access to a linked
   * list, similar to an STL bidirectional iterator. An * operator returns 
   * a reference to an associated Data object.  The ++ and -- operators 
   * change the current pointer to the next or prev element in a list. 
   *
   * The isEnd() method returns true if either end of the list has already
   * been passed by a previous ++ or -- operation. When isEnd() is true, the 
   * iterator is no longer usable, since it no longer points to a Node and
   * cannot be incremented or decremented.
   *
   * \ingroup Iterator_Module
   */
   template <typename Data>
   class ListIterator
   {

   public:

      /**
      * Default constructor.
      *
      * Creates a "dead" iterator, for which isEnd()==true. Before it
      * can be used, such an iterator must be initialized by either the
      * ListIterator<Data>::setCurrent() method or the List<Data>::begin() 
      * method of an associated List.
      */
      ListIterator()
       : current_(0)
      {}

      /**
      * Constructor for initialized iterator.
      *
      * Creates an iterator that points to the front of a List.
      * Calls List<Data>::begin(*this) internally.
      *
      * \param list parent List
      */
      explicit ListIterator(const List<Data>& list)
       : current_(0)
      { list.begin(*this); }

      /**
      * Point the iterator at a Node.
      *
      * \param nodePtr pointer to current Node in a List, or null.
      */
      void setCurrent(Node<Data> *nodePtr)
      { current_ = nodePtr; }

      /**
      * Has the end of the list been reached?
      *
      * Return true if the current pointer is null, indicating that
      * the previous increment or decrement passed an end of the list.
      *
      * \return true if current node is null, false otherwise.
      */
      bool isEnd() const
      { return (current_ == 0); }

      /**
      * Is this the back of the List?
      *
      * \return true if current node is the back Node, false otherwise.
      */
      bool isBack() const
      { return (current_->next() == 0); }

      /**
      * Is this the front of the List?
      *
      * \return true if current node is the front Node, false otherwise.
      */
      bool isFront() const
      { return (current_->previous() == 0); }

      /// \name Operators
      //@{

      /**
      * Get a const reference to the associated Data object.
      *
      * \return const reference to the associated Data object
      */
      const Data& operator* () const
      { return current_->data(); }

      /**
      * Get the associated Data object.
      *
      * \return reference to associated Data object
      */
      Data& operator* ()
      { return current_->data(); }

      /**
      * Get a pointer to const to the associated Data object.
      *
      * \return pointer to associated Data object
      */
      const Data * operator -> () const
      { return &(current_->data()); }

      /**
      * Get a pointer to the associated Data object.
      *
      * \return pointer to associated Data object
      */
      Data* operator -> ()
      { return &(current_->data()); }

      /**
      * Go to the next element in a linked list.
      *
      * This method assigns the current pointer to the address of the
      * next Node in the list, and then returns *this. If there is no
      * next Node, the current pointer is set null, and any subsequent
      * call to isEnd() will return true. 
      *
      * \return this ListIterator, after modification.
      */
      ListIterator<Data>& operator++ ()
      {
         current_ = current_->next();
         return *this;
      }

      /**
      * Go to the previous element in a linked list.
      *
      * This method assigns the current Node to the previous in the List,
      * and returns a reference to *this.
      *
      * \return this ListIterator
      */
      ListIterator<Data>& operator-- ()
      {
         current_ = current_->prev();
         return *this;
      }

      //@}
      
   private:

      /// Pointer to current Node.
      Node<Data>* current_;

   };

} 
#endif

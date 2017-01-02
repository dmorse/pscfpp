#ifndef UTIL_NODE_H
#define UTIL_NODE_H

/*
* Simpatico - Simulation Package for Polymeric and Molecular Liquids
*
* Copyright 2010 - 2014, The Regents of the University of Minnesota
* Distributed under the terms of the GNU General Public License.
*/

namespace Util
{

   template <typename Data> class List;

   /**
   * Linked List Node, class template.
   *
   * \ingroup List_Module
   */
   template <typename Data>
   class Node
   {

   public:

      /**
      * Default constructor.
      */
      Node()
       : prev_(0),
         next_(0),
         list_(0),
         datum_()
      {}

      /**
      * Copy constructor.
      */
      Node(const Node<Data>& other)
       : prev_(other.prev_),
         next_(other.next_),
         list_(other.list_),
         datum_(other.datum_)
      {}

      /**
      * Get the next pointer.
      *
      * \return pointer to next Node in List.
      */
      Node<Data>* next() const
      { return next_; }

      /**
      * Get the previous pointer.
      *
      * \return pointer to previous Node in List.
      */
      Node<Data>* prev() const
      { return prev_; }

      /**
      * Get a const reference to the associated Data.
      *
      * \return Data object associated with this Node.
      */
      const Data& data() const
      { return datum_; }

      /**
      * Get a reference to the associated Data object.
      *
      * \return Data object associated with this Node.
      */
      Data& data()
      { return datum_; }

      /**
      * Get a reference to the List.
      *
      * \return Reference to the list to which this Node belongs.
      */
      List<Data>& list() const
      { return *list_; }

      /**
      * Set pointer to the next Node.
      *
      * \param next pointer to next Node
      */
      void setNext(Node<Data>* next)
      { next_ = next; }

      /**
      * Set pointer to the previous Node.
      *
      * \param prev pointer to previous Node
      */
      void setPrev(Node<Data>* prev)
      { prev_ = prev; }

      /**
      * Set the list.
      *
      * \param list associated List object
      */
      void setList(List<Data>& list)
      { list_ = &list; }

      /**
      * Set the list.
      *
      * \param list pointer to an associated List object
      */
      void setList(List<Data>* list)
      { list_ = list; }

      /**
      * Set pointers connecting the other node after this node.
      *
      * This method sets the next pointer of this node to other and
      * the previous node of other to this. It also also sets the list
      * of the other node to this list.
      *
      * It does not reset the next pointer of the other node, and so
      * does not finish splicing the other node into the middle of this
      * list. The next pointer of the other node is left unchanged.
      *
      * \param other a Node to connect to this one.
      */
      void attachNext(Node<Data>& other)
      {
         next_       = &other;
         other.prev_ = this;
         other.list_ = list_;
      }

      /**
      * Set pointers connecting the other node before this node.
      *
      * This method sets the previous pointer of this node to other and
      * the next node of other to this. It also also sets the list of
      * the other node to this list.
      *
      * It does not reset the previous pointer of the other node, and so
      * does not finish splicing the other node into the middle of this
      * list. The previous pointer of the other node is left unchanged.
      *
      * \param other a Node to connect to this one.
      */
      void attachPrev(Node<Data>& other)
      {
         prev_          = &other;
         other.next_   = this;
         other.list_   = list_;
      }

      /**
      * Nullify previous and next pointers, and nullify the list pointer.
      *
      * This method disconnects the Node from any List, but does not modify
      * the datum_.
      */
      void clear()
      {
         prev_ = 0;
         next_ = 0;
         list_ = 0;
      }

   private:

      /**
      * Pointer to previous Node in List.
      *
      * This should be null for the front node.
      */
      Node<Data>*   prev_;

      /**
      * Pointer to next Node in List.
      *
      * This should be null for the back node.
      */
      Node<Data>*   next_;

      /**
      * Pointer to an associated List object.
      *
      * This should be null when the Node is not part of a list.
      */
      List<Data>*   list_;

      /**
      * Pointer to an associated Data object.
      *
      * Arrays of Data and Node objects are created by an ArrayList, with
      * the same array index for associated Node and Data elements.
      */
      Data          datum_;

   };

} 
#endif

#ifndef UTIL_LIST_H
#define UTIL_LIST_H

/*
* Simpatico - Simulation Package for Polymeric and Molecular Liquids
*
* Copyright 2010 - 2014, The Regents of the University of Minnesota
* Distributed under the terms of the GNU General Public License.
*/

#include <util/containers/Node.h>
#include <util/containers/ListIterator.h>
#include <util/global.h>

class ListTest;

namespace Util
{

   template <typename Data> class ListIterator;

   /**
   * Linked list class template.
   *
   * This list implementation is based on an underlying C array of Node<Data>
   * objects.  This array may be used by several List objects, and so must be
   * allocated outside the Link class and provided via the initialize method.
   *
   * \ingroup List_Module
   */
   template <typename Data>
   class List
   {

   public:

      /**
      * Default constructor.
      */
      List();

      /**
      * Destructor (does nothing).
      */
      virtual ~List()
      {}

      /**
      * Provide an array of Node<Data> objects for this List.
      */
      void initialize(Node<Data>* nodes, int capacity);

      /**
      * Get the number of elements.
      *
      * \return Number of elements in this list.
      */
      int size() const;

      /**
      * Get capacity of the array.
      *
      * \return Number of elements allocated in the associated arrays.
      */
      int capacity() const;

      /**
      * Push a node onto the the back of the List.
      *
      * \param node Node object from associated node array.
      */
      void pushBack(Node<Data>& node);

      /**
      * Push a node onto the the front of the List.
      *
      * \param node Node object from associated node array.
      */
      void pushFront(Node<Data>& node);

      /**
      * Remove a node from the back of the list.
      *
      * \return Node that was removed from this list.
      */
      Node<Data>& popBack();

      /**
      * Remove a node from the front of the list.
      *
      * \return Node that was removed from this list.
      */
      Node<Data>& popFront();

      /**
      * Insert newNode into list after node.
      *
      * \param node    Node in the existing list.
      * \param newNode new node, to be inserted as the next after node.
      */
      void insertNext(Node<Data>& node, Node<Data>& newNode);

      /**
      * Insert newNode into list before node.
      *
      * \param node    Node in the existing list.
      * \param newNode new Node, to be inserted previous to node.
      */
      void insertPrev(Node<Data>& node, Node<Data>& newNode);

      /**
      * Insert node into list in sequential order.
      *
      * \param node Node to be inserted into the list.
      */
      void insert(Node<Data> &node);

      /**
      * Remove node from list.
      *
      * \param node  Node to be removed from the list.
      */
      void remove(Node<Data> &node);

      /**
      * Set an iterator to the front of this List.
      *
      * \param iterator ListIterator, initialized on output.
      */
      void begin(ListIterator<Data> &iterator) const;

      /**
      * Check validity of linked list.
      *
      * \return true if the list is valid, false otherwise.
      */
      bool isValid() const;

   private:

      /// Array of nodes.
      Node<Data>*  nodes_;

      /// Pointer to front node of list.
      Node<Data>*  front_;

      /// Pointer to back node of list.
      Node<Data>*  back_;

      /// Pointer to lowest address node of list (or last element of array).
      Node<Data>*  lower_;

      /// Pointer to highest address node of list (or first element of array).
      Node<Data>*  upper_;

      /// Number of elements allocated for array nodes_ .
      int capacity_;

      /// Number of elements currently in this linked list.
      int size_;

      /**
      * Set List to empty state.
      */
      void setEmpty();

      /**
      * Update upper_ and lower_ bounds for a new node.
      *
      * \param node new Node being added.
      */
      void expandBounds(Node<Data>& node);

      /**
      * Update upper_ and lower_ bounds to reflect removal of a node.
      *
      * \param node Node being removed.
      */
      void contractBounds(Node<Data>& node);

   //friends

      friend class ::ListTest;

   }; // end class List


   /*
   * Default constructor.
   */
   template <typename Data>
   List<Data>::List()
    : nodes_(0),
      front_(0),
      back_(0),
      lower_(0),
      upper_(0),
      capacity_(0),
      size_(0)
   {}

   /*
   * Provide a C array of Node<Data> objects for this List.
   */
   template <typename Data>
   void List<Data>::initialize(Node<Data>* nodes, int capacity)
   {
      nodes_    = nodes;
      capacity_ = capacity;
      setEmpty();
   }

   /*
   * Get number of elements in list.
   */
   template <typename Data>
   inline int List<Data>::size() const
   {  return size_; }

   /*
   * Get capacity of underlying array.
   */
   template <typename Data>
   inline int List<Data>::capacity() const
   {  return capacity_; }

   /*
   * Push a node onto the the back of the List.
   */
   template <typename Data>
   void List<Data>::pushBack(Node<Data>& node)
   {
      if (back_ == 0) {
         node.setList(*this);
         front_ = &node;
         lower_ = &node;
         upper_ = &node;
      } else {
         back_->attachNext(node);
         expandBounds(node);
      }
      back_ = &node;
      node.setNext(0);
      ++size_;
   }

   /*
   * Push a node onto the the front of the List.
   */
   template <typename Data>
   void List<Data>::pushFront(Node<Data>& node)
   {
      if (front_ == 0) {
         node.setList(*this);
         front_ = &node;
         back_ = &node;
         lower_ = &node;
         upper_ = &node;
         size_  = 1;
      } else {
         front_->attachPrev(node);
         front_ = &node;
         expandBounds(node);
         ++size_;
      }
      node.setPrev(0);
   }

   /*
   * Remove a node from the back of the list.
   */
   template <typename Data>
   Node<Data>& List<Data>::popBack()
   {
      if (size_ == 0) {
         UTIL_THROW("Attempted popBack from empty List");
      }
      assert(back_  != 0);
      assert(front_ != 0);
      Node<Data>& oldBack = *back_;
      if (back_ == front_) {
         back_->clear();
         setEmpty();
      } else {
         back_ = oldBack.prev();
         back_->setNext(0);
         contractBounds(oldBack);
         oldBack.clear();
         --size_;
      }
      return oldBack;
   }

   /*
   * Remove a node from the back of the list.
   */
   template <typename Data>
   Node<Data>& List<Data>::popFront()
   {
      if (size_ == 0) {
         UTIL_THROW("Attempted popFront from empty List");
      }
      assert(front_ != 0);
      assert(back_  != 0);
      Node<Data>& oldFront = *front_;
      if (front_ == back_) {
         front_->clear();
         setEmpty();
      } else {
         front_ = oldFront.next();
         front_->setPrev(0);
         contractBounds(oldFront);
         oldFront.clear();
         --size_;
      }
      return oldFront;
   }

   /*
   * Insert newNode into the list after previous.
   */
   template <typename Data>
   void List<Data>::insertNext(Node<Data>& previous, Node<Data>& newNode)
   {
      if (&previous == back_) {
         back_ = &newNode;
         newNode.setNext(0);
      } else {
         previous.next()->attachPrev(newNode);
      }
      previous.attachNext(newNode);
      expandBounds(newNode);
      ++size_;
   }

   /*
   * Insert newNode into the list before next.
   */
   template <typename Data>
   void List<Data>::insertPrev(Node<Data>& next, Node<Data>& newNode)
   {
      if (&next == front_ ) {
         front_ = &newNode;
         newNode.setPrev(0);
      } else {
         next.prev()->attachNext(newNode);
      }
      next.attachPrev(newNode);
      expandBounds(newNode);
      ++size_;
   }

   /*
   * Remove a specific node from the list.
   */
   template <typename Data>
   void List<Data>::remove(Node<Data>& node)
   {
      assert( &node.list() == this );
      if (&node == back_) {
         if (back_ == front_) {
            back_->clear();
            setEmpty();
            return;
         } else {
            back_ = node.prev();
            back_->setNext(0);
         }
      } else
      if (&node == front_) {
         front_ = node.next();
         front_->setPrev(0);
      } else {
         node.prev()->setNext( node.next() );
         node.next()->setPrev( node.prev() );
      }
      contractBounds(node);
      node.clear();
      --size_;
   }

   /*
   * Insert one new node.
   */
   template <typename Data>
   void List<Data>::insert(Node<Data> &node)
   {
      //List<Data>& list = node.list();
      Node<Data>* current;
      Node<Data>* target;

      // Assert that node does not already belong to this list
      assert( &node.list() != this);

      // If list is empty
      if (front_ == 0) {
         node.setList(*this);
         front_ = &node;
         back_ = &node;
         lower_ = &node;
         upper_ = &node;
         size_  = 1;
         return;
      }

      // If node address is above highest address
      if (&node > upper_) {
         target = upper_->next();
         upper_->attachNext(node);
         upper_ = &node;
      } else // If node address is below lowest address
      if (&node < lower_) {
         target = lower_->prev();
         if (target) {
            target->attachNext(node);
         } else {
            node.setPrev(0);
            node.setList(*this);
            front_ = &node;
         }
         target = lower_;
         lower_ = &node;
      } else // If between lowest and highest addresses
      {
         current = &node;
         do {
            assert(current > lower_);
            --current;
         } while (&current->list() != this);
         target = current->next();
         current->attachNext(node);
      }

      // Set Next for new Node, or set to back_
      if (target != 0) {
         node.attachNext(*target);
      } else {
         node.setNext(0);
         back_ = &node;
      }

      ++size_;
   }

   /*
   * Set iterator to front of list.
   */
   template <typename Data>
   void List<Data>::begin(ListIterator<Data> &iterator) const
   {
      assert(front_ != 0);  
      iterator.setCurrent(front_); 
   }

   /*
   * Check validity of linked list.
   */
   template <typename Data>
   bool List<Data>::isValid() const
   {
      if (front_ == 0) {

         if (size_  != 0) {
            UTIL_THROW("List<Data>::isValid: front_ == 0 and size_ != 0.");
            return false;
         }

         if (back_  != 0) {
            UTIL_THROW("List<Data>::isValid:  front_ == 0 and back_ != 0.");
            return false;
         }

         if (lower_ != nodes_ + capacity_) {
            UTIL_THROW("front_ == 0 and lower_ != nodes_+capacity.");
            return false;
         }

         if (upper_ != nodes_ -1) {
            UTIL_THROW("List<Data>::isValid: front_ == 0 and upper_ != nodes_-1.");
            return false;
         }

      } else {

         int i=1;
         Node<Data>* prev = 0;
         Node<Data>* node = front_;
         while (node->next() != 0 ) {
            if (node->prev() != prev) {
               UTIL_THROW("List<Data>::isValid:  node->prev() != prev.");
               return false;
            }
            if (&node->list() != this) {
               UTIL_THROW("List<Data>::isValid: &node->list() != this.");
               return false;
            }
            if (node > upper_) {
               UTIL_THROW("List<Data>::isValid: node > upper_.");
               return false;
            }
            if (node < lower_) {
               UTIL_THROW("List<Data>::isValid: node < lower_.");
               return false;
            }
            prev = node;
            node = node->next();
            ++i;
         }

         if (i < size_) {
            UTIL_THROW("List<Data>::isValid: # elements < size_.");
            return false;
         }
         if (i > size_) {
            UTIL_THROW("List<Data>::isValid: # elements > size_.");
            return false;
         }
         if (node != back_) {
            UTIL_THROW("List<Data>::isValid: Last node != back_.");
            return false;
         }
         if (node->next() != 0) {
            UTIL_THROW("List<Data>::isValid: Last node, node->next() != 0.");
            return false;
         }

      }
      // If no errors were detected to this point, the list is valid.
      return true;
   }

   // Private methods
   
   /*
   * Set list to empty state.
   */
   template <typename Data>
   void List<Data>::setEmpty()
   {
      size_  = 0;
      front_ = 0;
      back_  = 0;

      // Set lower above all addresses and upper below all so that the
      // expandBounds method will reset both values when a node is added.
      lower_  = nodes_ + capacity_ ;
      upper_  = nodes_ - 1;

   }

   /*
   * Update upper_ and lower_ bounds for a new node.
   */
   template <typename Data>
   void List<Data>::expandBounds(Node<Data>& node)
   {
      if (&node > upper_) {
         upper_ = &node;
      }
      if (&node < lower_) {
         lower_ = &node;
      }
   }

   /*
   * Update upper_ and lower_ bounds to reflect removal of a node.
   */
   template <typename Data>
   void List<Data>::contractBounds(Node<Data>& node)
   {
      if (&node == upper_) {
         do {
            --upper_;
         } while ( &(upper_->list()) != this );
      }
      if (&node == lower_) {
         do {
            ++lower_;
         } while ( &(lower_->list()) != this );
      }
   }

} 
#endif

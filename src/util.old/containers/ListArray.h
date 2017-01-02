#ifndef UTIL_LIST_ARRAY_H
#define UTIL_LIST_ARRAY_H

/*
* Simpatico - Simulation Package for Polymeric and Molecular Liquids
*
* Copyright 2010 - 2014, The Regents of the University of Minnesota
* Distributed under the terms of the GNU General Public License.
*/

#include <util/containers/Array.h>
#include <util/containers/Node.h>
#include <util/containers/List.h>
#include <util/misc/Memory.h>
#include <util/global.h>

namespace Util
{

   /**
   * An array of objects that are accessible by one or more linked List objects.
   *
   * A ListArray is an allocatable array of data objects that also provides access
   * to some or all of its via one or more associated List objects. Each element 
   * of the array may be part of at most one List.
   *
   * \ingroup List_Module
   */
   template <typename Data>
   class ListArray 
   { 
   
   public: 
   
      /**
      * Constructor.
      */
      ListArray();
   
      /**
      * Destructor.
      *
      * Delete dynamically allocated arrays of Node and List objects.
      */
      virtual ~ListArray();
    
      /**
      * Allocate arrays of Node and List objects.
      *
      * \param capacity size of array Node<Data> objects
      * \param nList    size of array of List<Data> linked list objects
      */
      void allocate(int capacity, int nList);
   
      /**
      * Get the number of associated linked lists.
      *
      * \return size of array of associated List<Data> objects.
      */
      int nList() const;

      /**
      * Return allocated size of underlying array of nodes.
      *
      * \return Number of elements allocated in array.
      */
      int capacity() const
      { return capacity_; }

      /**
      * Return data for node element i.
      *
      * \param  i array index
      * \return reference to element i
      */
      Data& operator[] (int i)
      {
         assert(nodes_ != 0);
         assert(i < capacity_);
         assert(i >= 0);
         return nodes_[i].data();
      }

      /**
      * Return const refereence to Data in Node element number i.
      *
      * \param i array index
      * \return const reference to element i
      */
      const Data& operator[] (int i) const
      {
         assert(nodes_ != 0);
         assert(i < capacity_);
         assert(i >= 0 );
         return nodes_[i].data();
      }

      /**
      * Return a reference to a specific List.
      *
      * \param  i  array index
      * \return reference to List number i
      */
      List<Data>& list(int i);
   
      /**
      * Return a const reference to a specific List.
      *
      * \param  i  array index
      * \return reference to List number i
      */
      const List<Data>& list(int i) const;

      /**
      * Return reference to node number i.
      *
      * \param  i  array index
      * \return reference to Data object element number i
      */
      Node<Data>& node(int i);
   
      /**
      * Return true if the ListAray is valid, or throw an exception.
      */
      bool isValid() const;
   
   private:
  
      // C array of Node<Data> objects 
      Node<Data> *nodes_;
   
      // C array of List<Data> objects
      List<Data> *lists_;
   
      // Allocated size of nodes_ array.
      int         capacity_;
   
      // Number of lists.
      int         nList_;
   
   }; 


   /* 
   * Default constructor.
   */
   template <typename Data>
   ListArray<Data>::ListArray()
    : nodes_(0),    
      lists_(0),
      capacity_(0),
      nList_(0)
   {}

   /* 
   * Destructor.
   *
   * Delete dynamically allocated arrays of Data, Node, and List objects.
   */
   template <typename Data>
   ListArray<Data>::~ListArray()
   {
      if (nodes_) {
         Memory::deallocate< Node<Data> >(nodes_, capacity_);
      }
      if (lists_) {
         Memory::deallocate< List<Data> >(lists_, nList_);
      }
   }
 
   /* 
   * Allocate arrays of Data objects, list nodes, and linked lists.
   *
   * \param capacity size of arrays of Data and ListNode<Data> objects
   * \param nList    size of array of List<Data> linked list objects
   */
   template <typename Data>
   void ListArray<Data>::allocate(int capacity, int nList) 
   {
      int i;
      capacity_ = capacity;
      nList_    = nList;

      // Allocate array of nodes
      Memory::allocate< Node<Data> >(nodes_, capacity_);

      // Allocate and initialize array of lists
      Memory::allocate< List<Data> >(lists_, nList_);
      for (i=0; i < nList_; ++i) {
         lists_[i].initialize(nodes_, capacity_);
      }

   }

   /* 
   * Get the number of associated linked lists.
   *
   * \return size of array of associated List<Data> objects.
   */
   template <typename Data>
   inline int ListArray<Data>::nList() const
   { return nList_;}

   /* 
   * Return a reference to a specific List.
   *
   * \param  i  array index
   * \return reference to List number i
   */
   template <typename Data>
   List<Data>& ListArray<Data>::list(int i) 
   { 
      assert(lists_ != 0); 
      assert(i >= 0); 
      assert(i < nList_); 

      return *(lists_ + i); 
   }

   /* 
   * Return a const reference to a specific List.
   *
   * \param  i  array index
   * \return reference to List number i
   */
   template <typename Data>
   const List<Data>& ListArray<Data>::list(int i) const
   { 
      assert(lists_ != 0); 
      assert(i >= 0); 
      assert(i < nList_); 

      return *(lists_ + i); 
   }

   /* 
   * Return reference to node number i.
   *
   * \param  i  array index
   * \return reference to Data object element number i
   */
   template <typename Data>
   Node<Data>& ListArray<Data>::node(int i) 
   { 
      assert(nodes_ != 0); 
      assert(i >= 0); 
      assert(i < capacity_); 

      return *(nodes_ + i); 
   }

   /*
   * Return true if the ListAray is valid, or throw an exception.
   */
   template <typename Data>
   bool ListArray<Data>::isValid() const
   {
      if (nodes_ != 0) {
         if (lists_ == 0) {
            UTIL_THROW("nodes_ is allocated but lists_ is not");
         }

         // Check validity of all lists 
         for (int i=0; i < nList_ ; ++i) {
            if (!lists_[i].isValid()){
               UTIL_THROW("Invalid list in ListArray");
            }
         }

      } else {

         if (lists_ != 0) {
            UTIL_THROW("nodes_ is not allocated but lists_ is");
         }

      }
      return true;

   }

} 
#endif

#ifndef UTIL_R_ARRAY_H
#define UTIL_R_ARRAY_H

/*
* Simpatico - Simulation Package for Polymeric and Molecular Liquids
*
* Copyright 2010 - 2014, The Regents of the University of Minnesota
* Distributed under the terms of the GNU General Public License.
*/

#include <util/containers/Array.h>
#include <util/containers/FSArray.h>
#include <util/global.h>

namespace Util
{

   /**
   * An Array that acts as a reference to another Array or C array.
   *
   * An RArray is associated with a "target" DArray or C array by the 
   * associate() method. The RArray and its target array then wrap 
   * the same underlying C array, and so access the same data. The
   * associate() method simply copies the address and capacity of 
   * a C array. An RArray can be associated only once, 
   * after which it can be safely used as an alias for its target.
   *
   * An RArray can only be associated with a DArray after the target
   * DArray has been allocated. Because a DArray can be allocated 
   * only once, this association cannot be corrupted by re-allocation 
   * or re-sizing of the target DArray. 
   *
   * An RArray can be created from another RArray only after the
   * target RArray has already been associated with some other Array.
   *
   * An RArray differs from a C++ reference to an Array because a
   * C++ reference must be initialized when it is instantiated,
   * whereas an RArray is associated after it is instantiated.
   * Because association is implemented by copying the address and 
   * capacity of a shared C array, access through an RArray should 
   * be exactly as efficient as access through a DArray.
   *
   * \ingroup Array_Module
   */
   template <typename Data>
   class RArray : public Array<Data>
   { 

      using Array<Data>::data_;
      using Array<Data>::capacity_;
   
   public: 
   
      /**
      * Constructor.
      */
      RArray() : 
         Array<Data>()
      {}
   
      /**
      * Copy constructor.
      *
      * Shallow copy of another RArray
      *
      * \param other another RArray<Data> for which this is an alias.
      */
      RArray(const RArray<Data>& other) : 
         Array<Data>()
      {
         data_     = other.data_;
         capacity_ = other.capacity_;
      }
  
      /**
      * Associate this RArray with an existing Array object.
      *
      * The target (i.e., the parameter array) must be allocated when this
      * method is invoked, as discussed in the RArray class documentation.
      *
      * \param array the target Array
      */
      void associate(Array<Data> &array) 
      {
         if (data_ != 0) {
            UTIL_THROW("Attempt to re-associate an RArray");
         }
         if (array.capacity() <= 0) {
            UTIL_THROW("Unallocated target array: Capacity_ <= 0");
         }
         data_     =  &array[0];
         capacity_ =   array.capacity();
      }
 
      /**
      * Associate this RArray with an existing C array.
      *
      * \param array    the target C array
      * \param capacity the number of elements in the target array
      */
      void associate(Data* array, int capacity) 
      {
         if (data_ != 0) {
            UTIL_THROW("Attempt to re-associate an RArray");
         }
         data_     = array;
         capacity_ = capacity;
      }
 
   private:
 
      /**
      * Assignment, private and not implemented to prohibit operation.
      */
      RArray<Data>& operator=(const RArray<Data>& other);
   
   }; // end class RArray

} 
#endif

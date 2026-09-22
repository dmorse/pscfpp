#ifndef PSCF_CONST_HOST_ARRAY_CP_H
#define PSCF_CONST_HOST_ARRAY_CP_H

/*
* Util Package - C++ Utilities for Scientific Computation
*
* Copyright 2010 - 2017, The Regents of the University of Minnesota
* Distributed under the terms of the GNU General Public License.
*/

#include <util/containers/ConstArray.h>  // base class

// Forward declarations
namespace Pscf {
   template <typename Data, typename T> class ConstHostArray;
   template <typename Data, typename T> class DeviceArray;
   class CPT;
}

namespace Pscf {

   using namespace Util;

   /**
   * Read-only psuedo-"host" array for use with C++ backend.
   *
   * This class template partial specialization may be used in template
   * code in which data is copied from a psuedo-device array to a
   * pseudo-host array when the device array is const and/or if the host
   * array only requires read access. 
   *
   * Construction or assignment (operator =)from a DeviceArray<Data,CPT>
   * to a ConstHostArray<Data,CPT> creates a shallow read-only copy of a
   * C array that is owned by the device array. Assignment does nothing
   * if such an association already exists.  Because the relevant
   * constructor and assignment operator each take a const reference to
   * a DeviceArray<Data,CPT> as a parameter, they can be used in contexts
   * in which the device array is declared const. 
   *
   * The use of construction or an assignment (=) operator to create a
   * shallow copy allows this class to be used in template code to
   * imitate the syntax of an actual device-to-host data copy that
   * would be performed by a ConstHostArray<Data, CUT> in a template
   * specialization that is designed to use a GPU.
   *
   * An association with a device array can be released by calling the
   * dissociate function, or will be released upon destruction.  Such an 
   * association must be released by one of these two methods before the 
   * device array that owns the associated data is de-allocated or 
   * destroyed. A reference counting system detects and reports erroneous
   * de-allocation of a source array that is still referred to by one or 
   * more other array view containers.
   *
   * \ingroup Pscf_Backend_Cpp_Module
   */
   template <typename Data>
   class ConstHostArray<Data, CPT> : public ConstArray<Data>
   {

   public:

      /**
      * Default constructor.
      */
      ConstHostArray();

      /**
      * Copy construction from a DeviceArray.
      *
      * Create an association (shallow copy) with memory owned by
      * the pre-existing device array.
      *
      * \param other  array container
      */
      ConstHostArray(DeviceArray<Data,CPT> const & other);

      // Prohibit copy construction from another ConstHostArray.
      ConstHostArray(ConstHostArray<Data,CPT> const & other) = default;

      /**
      * Destructor.
      */
      ~ConstHostArray();

      // Prohibit assignment from another ConstHostArray.
      ConstHostArray<Data,CPT>&
      operator = (ConstHostArray<Data,CPT> const&) = delete;

      /**
      * Create read-only association with a DeviceArray.
      *
      * If this is already associated with other (the device array),
      * do nothing and return. Otherwise, create an association of this
      * with the other device array (i.e., create a shallow copy).
      *
      * \throw Exception if other array is not allocated
      * \throw Exception if this is associated with other with wrong size
      *
      * \param other  array container on RHS of assigment (input)
      */
      ConstHostArray<Data,CPT>& operator = (DeviceArray<Data,CPT> const & other);

      /**
      * Release association with a device array.
      */
      void dissociate();

   private:

      /// Reference to a device array that owns memory referenced by this.
      CountedReference ref_;

      /**
      * Associate this object with source array.
      *
      * This associates this object with a source array.
      *
      * \throw Exception if this array is already associated.
      * \throw Exception if source array is not allocated on entry.
      *
      * \param source  array that owns the data
      */
      void associate(DeviceArray<Data,CPT> const & source);

      using ConstArray<Data>::data_;
      using ConstArray<Data>::capacity_;

   };

} // namespace Pscf

#include <pscf/backend/cpp/DeviceArray.h>
#include <util/global.h>

namespace Pscf {

   // Member functions

   /*
   * Default constructor.
   */
   template <typename Data>
   ConstHostArray<Data,CPT>::ConstHostArray()
    : ConstArray<Data>(),
      ref_()
   {}

   /*
   * Destructor.
   */
   template <typename Data>
   ConstHostArray<Data,CPT>::~ConstHostArray()
   {
      if (ref_.isAssociated()) {
         ref_.dissociate(); // decrements counter of source array
      }
   }

   /*
   * Copy construction from a DeviceArray.
   *
   * Creates an association with a DeviceArray, which must own the data.
   */
   template <typename Data>
   ConstHostArray<Data,CPT>::ConstHostArray(
                                    DeviceArray<Data,CPT> const & other)
    : ConstArray<Data>()
   {  associate(other); }

   /*
   * Assignment from a DeviceArray.
   *
   * Creates an association with a DeviceArray, which must own the data.
   */
   template <typename Data>
   ConstHostArray<Data,CPT>&
   ConstHostArray<Data,CPT>::operator = (DeviceArray<Data,CPT> const & other)
   {
      UTIL_CHECK(other.isAllocated());
      Data const * data = ConstArray<Data>::cArray();
      if ((bool)data && other.cArray() == data) {
         UTIL_CHECK(other.capacity() == ConstArray<Data>::capacity());
      } else {
         associate(other);
      }
      return *this;
   }

   /*
   * Release association with a device array.
   */
   template <typename Data>
   void ConstHostArray<Data,CPT>::dissociate()
   {
      UTIL_CHECK(data_);
      UTIL_CHECK(ref_.isAssociated());

      data_ = nullptr;
      capacity_ = 0;
      ref_.dissociate(); // decrements reference counter of source array
   }

   // Private function

   /*
   * Associate this object with a device array.
   */
   template <typename Data>
   void ConstHostArray<Data,CPT>::associate(DeviceArray<Data,CPT> const & source)
   {
      UTIL_CHECK(source.isAllocated());
      UTIL_CHECK(!ref_.isAssociated());

      // Copy data pointer and size
      data_ = const_cast<Data*>( source.cArray() );
      capacity_ = source.capacity();

      // Note: const_cast to non-const pointer is permissible because
      // the ConstArray<Data> class public interface is designed to
      // prevent modification of individual array elements.

      // Associate ReferencecCounter base class of the source array with
      // the CountedReference ref_ member variable of this data user.
      ref_.associate(source);

      // On exit, the ReferenceCounter of the data source is incremented
      // and the ref_ CountedReference member variable of this object
      // holds a pointer to that ReferenceCounter.
   }

}
#endif

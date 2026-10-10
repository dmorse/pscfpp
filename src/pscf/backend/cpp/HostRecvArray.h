#ifndef PSCF_HOST_RECV_ARRAY_CP_H
#define PSCF_HOST_RECV_ARRAY_CP_H

/*
* Util Package - C++ Utilities for Scientific Computation
*
* Copyright 2010 - 2017, The Regents of the University of Minnesota
* Distributed under the terms of the GNU General Public License.
*/

#include <util/containers/ConstArray.h>  // base class
#include <util/misc/CountedReference.h>  // member
#include <pscf/backend/cpp/CPT.h>        // template argument

// Forward declarations
namespace Pscf {
   template <typename Data, typename T> class DeviceArray;
}

namespace Pscf {

   using namespace Util;

   // Declare primary template
   template <typename Data, typename T> class HostRecvArray;

   /**
   * Array for pseudo device-to-host array transfer with C++ backend.
   *
   * This template specialization should be used in backend-independent 
   * template code as the host array in contexts in which the specialization
   * for the CUDA backend (T=Pscf:CUT) would copy an array from the GPU
   * device to the CPU host. This specialization for the C++ backend
   * (T=Pscf::CPT) instead makes the host array a shallow copy of the
   * device array, thus avoiding the cost of an actual deep copy. This
   * specialization supports random read-only access to array elements 
   * via a subscript ([]) operator that returns a const reference, or 
   * via a ConstArrayIterator. 
   *
   * Conversion construction and assignment (operator =) from an instance 
   * of DeviceArray<Data,CPT> to a HostRecvArray<Data,CPT> each create a 
   * shallow read-only copy of a C array that is owned by the associated
   * device array, via a private pointer to this shared array.  These two 
   * functions each take a const reference to DeviceArray<Data,CPT> as a 
   * parameter, and so can be used in contexts in which the device array
   * must be passed as a const reference. Conversion construction is
   * equivalent to default construction followed by assignment.
   *
   * An association with a device array can be released by calling the
   * dissociate member function, or will be released upon destruction of
   * this host array. Such an association must be released by one of these 
   * two methods before the device array that owns the associated data is 
   * de-allocated or destroyed. A reference counting system is used to 
   * detect and report dangerous de-allocation of a device array that 
   * is still referred to by one or more other containers, which would
   * create a dangling reference.
   *
   * <b> Usage: </b>
   *
   * Typical usage is shown below for a pseudo device-to-host transfer 
   * from a longer lived instance of DeviceArray<Data,CPT> named 
   * deviceArray to a shorter lived instance of HostRecvArray<Data,CPT> 
   * named hostArray. The alias Data denotes the data type of each array 
   * element.
   *
   * \code
   *    HostRecvArray<Data,CPT> hostArray;
   *    hostArray = deviceArray
   *
   *    // Read the data in hostArray, and possibly peform a computation 
   *
   *    hostArray.dissociate();
   * \endcode
   * The default constructor and assignment operations may be combined 
   * into a single call of the conversion constructor, giving the shorter 
   * version:
   * \code
   *    HostRecvArray<Data,CPT> hostArray(deviceArray)
   *
   *    // Read the data in hostArray, and possibly peform a computation 
   *
   *    hostArray.dissociate();
   * \endcode
   *
   * Comments:
   * 
   *  - The assignment (=) operator or conversion constructor each make 
   *    hostArray a shallow copy of deviceArray, which references an array
   *    that is owned by deviceArray.
   *
   *  - The dissociate member function releases the association with the 
   *    device array, by nullifying the private Data* pointer member of 
   *    host array that is used to point to an array.
   *
   *  - A call to the dissociate function call must occur after all 
   *    statements that access elements of the host array. 
   *
   *  - The call to the dissociate function may be omitted if hostArray is 
   *    a local variable that is defined within a function but deviceArray 
   *    lives beyond the end of the function. In this case, any remaining
   *    association will be released by the host array destructor when 
   *    function ends and the host array is destroyed.
   *    association.
   *
   *  - The HostRecvArray template does not define an public 
   *    "associate" member function like that defined by the HostSendArray
   *    template. By convention, an association may only be created by
   *    the conversion constructor or the assignment operator. 
   *
   * \see Pscf::HostRecvArray<Data, CUT>
   * \ingroup Pscf_Backend_Cpp_Module
   */
   template <typename Data>
   class HostRecvArray<Data, CPT> final : public ConstArray<Data>
   {

   public:

      /**
      * Backend identifier class typename alias.
      */
      using BackendIdClass = CPT;

      // Default constructor.
      HostRecvArray() = default;

      // Copy construction (delete)
      HostRecvArray(HostRecvArray<Data,CPT> const & other) = delete;

      /**
      * Destructor.
      *
      * Releases any remaining association.
      */
      ~HostRecvArray();

      // Assignment from another HostRecvArray (delete)
      HostRecvArray<Data,CPT>&
      operator = (HostRecvArray<Data,CPT> const&) = delete;

      /**
      * Conversion construction from a DeviceArray.
      *
      * Create a read-only association with the device array, thus
      * creating a shallow copy that points to the memory owned by
      * the device array.
      *
      * \param other  array container
      */
      HostRecvArray(DeviceArray<Data,CPT> const & other);

      /**
      * Create a read-only association with a DeviceArray.
      *
      * Create an association with the other device array, thus
      * creating a shallow copy that points to the memory owned by
      * the device array. This array may not be allocated on entry.
      *
      * \throw Exception if other array is not allocated
      * \throw Exception if this is already associated with an array
      *
      * \param other  array container on RHS of assigment (input)
      */
      HostRecvArray<Data,CPT>& 
      operator = (DeviceArray<Data,CPT> const & other);

      /**
      * Release association with a device array.
      *
      * \throw Exception if no such association exists.
      */
      void dissociate();

   private:

      /// Reference to a device array that owns memory referenced by this.
      CountedReference ref_;

      /**
      * Associate this object with source array.
      *
      * This private function associates this object with a source array.
      * It is called by the conversion constructor and assignment
      * operator.
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

   // Explicit instantiation declarations
   extern template class HostRecvArray<double,CPT>;
   extern template class HostRecvArray<fftw_complex,CPT>;

} // namespace Pscf

#include <pscf/backend/cpp/DeviceArray.h>
#include <util/global.h>

namespace Pscf {

   // Public member function definitions

   /*
   * Destructor.
   */
   template <typename Data>
   HostRecvArray<Data,CPT>::~HostRecvArray()
   {
      if (ref_.isAssociated()) {
         ref_.dissociate(); // decrements counter of source array
      }
   }

   /*
   * Copy construction from a DeviceArray.
   *
   * Creates an association with a DeviceArray that owns the data.
   */
   template <typename Data>
   HostRecvArray<Data,CPT>::HostRecvArray(
                                    DeviceArray<Data,CPT> const & other)
    : ConstArray<Data>()
   {  associate(other); }

   /*
   * Assignment from a DeviceArray.
   *
   * Creates an association with a DeviceArray that owns the data.
   */
   template <typename Data>
   HostRecvArray<Data,CPT>&
   HostRecvArray<Data,CPT>::operator = (DeviceArray<Data,CPT> const& other)
   {
      associate(other);
      return *this;
   }

   /*
   * Release association with a device array.
   */
   template <typename Data>
   void HostRecvArray<Data,CPT>::dissociate()
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
   void 
   HostRecvArray<Data,CPT>::associate(DeviceArray<Data,CPT> const & source)
   {
      UTIL_CHECK(source.isAllocated());
      UTIL_CHECK(!ref_.isAssociated());

      // Copy data pointer and size
      data_ = const_cast<Data*>( source.cArray() );
      capacity_ = source.capacity();

      // Note: The const_cast to non-const pointer is permissible because
      // the public interface of the ConstArray<Data> base class does
      // not allow modification of array elements.

      // Associate ReferenceCounter of the source array with this 
      ref_.associate(source);

      // On exit, the ReferenceCounter sub-object of the source array is 
      // incremented and the ref_ CountedReference member variable of 
      // this object holds a pointer to that ReferenceCounter.
   }

}
#endif

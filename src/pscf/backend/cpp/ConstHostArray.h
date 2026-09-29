#ifndef PSCF_CONST_HOST_ARRAY_CP_H
#define PSCF_CONST_HOST_ARRAY_CP_H

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
   template <typename Data, typename T> class ConstHostArray;

   /**
   * Read-only pseudo-"host" array for use with C++ backend.
   *
   * This template specialization should be used in backend-independent 
   * template code in which data is copied from a device array to a
   * a host array. This is derived from ConstArray<Data>, and thus 
   * provides read-only access to the underlying data.
   *
   * Conversion construction and assignment (operator =) from an instance 
   * of DeviceArray<Data,CPT> to a ConstHostArray<Data,CPT> each create a 
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
   * The use of construction or assignment to create a shallow read-only 
   * copy allows this class to be used in backend-independent template 
   * code to imitate the syntax of an actual device-to-host data copy 
   * that would be performed by analogous functions in the specialization
   * for a CUDA backend, without performing an unnecessary deep copy. 
   *
   * <b> Usage: </b>
   *
   * Typical usage is shown below for device-to-host transfer from a longer
   * lived instance of DeviceArray<Data,CPT> named dArray (the device array)
   * to a shorter lived instance of ConstHostArray<Data,CPT> named hArray
   * (the host array). The alias Data denotes the data type of each array 
   * element.
   *
   * \code
   *    ConstHostArray<Data,CPT> hArray;
   *    hArray = dArray
   *
   *    // (Read and use the data in dArray)
   *
   *    hArray.dissociate();
   * \endcode
   * The default constructor and assignment operations may be combined 
   * into a single call of the conversion constructor, giving the shorter 
   * version:
   * \code
   *    ConstHostArray<Data,CPT> hArray(dArray)
   *
   *    // Read and use the data in dArray
   *
   *    hArray.dissociate();
   * \endcode
   *
   * Comments:
   * 
   *  - In this specialization for a C++ backend (T=CPT), the assignment
   *    operator or conversion constructor creates a shallow copy. In 
   *    the corresponding specialization for a GPU backend (T=CUT), these
   *    functions would instead allocate the host array if not already
   *    allocated, and then perform a deep copy of data from GPU device 
   *    memory to CPU host memory.
   *
   *  - In this specialization for a C++ backend (T=CPT), the dissociate
   *    member function release the association with the device array, by 
   *    nullifying a pointer owned by the host array. In the specialization
   *    for a CUDA backend (T=CUT), the dissociate function does nothing.
   *
   *  - The dissociate function call must occur after all statements 
   *    that access elements of the host array. In the specialization for
   *    a C++ backend (T=CPT), the host array has no access to data after 
   *    the dissociate function is called.
   *
   *  - The call to the dissociate function may be omitted if hArray is 
   *    a local variable that is defined within a function but dArray 
   *    lives beyond the end of the function. Any remaining association
   *    will be deleted by the host array destructor when the host array 
   *    goes out of scope. There is, however, also no harm in calling 
   *    dissociate explicitly to document the need to release the 
   *    association.
   *
   *  - If data that is copied from device to host needs to be modified
   *    on the host, it must first be copied to a container that provides
   *    read-write access. A ConstHostArray container only provides read
   *    access.
   * 
   *  - The ConstHostArray template does not define an public 
   *    "associate" member function like that defined by the HostArray
   *    template. By convention, an association may only be created by
   *    the conversion constructor or the assignment operator. 
   *
   * \see Pscf::ConstHostArray<Data, CUT>
   * \ingroup Pscf_Backend_Cpp_Module
   */
   template <typename Data>
   class ConstHostArray<Data, CPT> final : public ConstArray<Data>
   {

   public:

      /**
      * Backend identifier class typename alias.
      */
      using BackendIdClass = CPT;

      // Default constructor.
      ConstHostArray() = default;

      // Copy construction (delete)
      ConstHostArray(ConstHostArray<Data,CPT> const & other) = delete;

      /**
      * Destructor.
      *
      * Releases any remaining association.
      */
      ~ConstHostArray();

      // Assignment from another ConstHostArray (delete)
      ConstHostArray<Data,CPT>&
      operator = (ConstHostArray<Data,CPT> const&) = delete;

      /**
      * Conversion construction from a DeviceArray.
      *
      * Create an association the other device array, thus
      * creating a shallow copy that points to the memory owned by
      * the device array.
      *
      * \param other  array container
      */
      ConstHostArray(DeviceArray<Data,CPT> const & other);

      /**
      * Create a read-only association with a DeviceArray.
      *
      * Create an association with the other device array, thus
      * creating a shallow copy that points to the memory owned by
      * the device array. This host array may not be allocated on
      * entry.
      *
      * \throw Exception if other array is not allocated
      * \throw Exception if this is already associated with an array
      *
      * \param other  array container on RHS of assigment (input)
      */
      ConstHostArray<Data,CPT>& 
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
   extern template class ConstHostArray<double,CPT>;
   extern template class ConstHostArray<fftw_complex,CPT>;

} // namespace Pscf

#include <pscf/backend/cpp/DeviceArray.h>
#include <util/global.h>

namespace Pscf {

   // Public member function definitions

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
   * Creates an association with a DeviceArray that owns the data.
   */
   template <typename Data>
   ConstHostArray<Data,CPT>::ConstHostArray(
                                    DeviceArray<Data,CPT> const & other)
    : ConstArray<Data>()
   {  associate(other); }

   /*
   * Assignment from a DeviceArray.
   *
   * Creates an association with a DeviceArray that owns the data.
   */
   template <typename Data>
   ConstHostArray<Data,CPT>&
   ConstHostArray<Data,CPT>::operator = (DeviceArray<Data,CPT> const& other)
   {
      associate(other);
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
   void 
   ConstHostArray<Data,CPT>::associate(DeviceArray<Data,CPT> const & source)
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

#ifndef PSCF_HOST_ARRAY_CP_H
#define PSCF_HOST_ARRAY_CP_H

/*
* Util Package - C++ Utilities for Scientific Computation
*
* Copyright 2010 - 2017, The Regents of the University of Minnesota
* Distributed under the terms of the GNU General Public License.
*/

#include <util/containers/Array.h>       // base class template
#include <util/misc/CountedReference.h>  // member
#include <pscf/backend/cpp/CPT.h>        // template argument

// Forward declaration
namespace Pscf {
   template <typename Data, typename T> class DeviceArray;
}

namespace Pscf {

   using namespace Util;

   // Declare primary template
   template <typename Data, typename T> class HostArray;

   /**
   * Psuedo-"host" array used for host-to-device data transfer.
   *
   * This class template should be used in device-independent template
   * code in which an array is copied from host to device.  The
   * DeviceArray<Data,CPT> template defines an assignment operator
   * that assigns from a host array to a device array. 
   *
   * The "associate" member function creates a shallow copy of a C 
   * array that is owned by an associated device array.  This function 
   * takes a non-const reference to the associated DeviceArray<Data,CPT> 
   * its only parameter. This parameter is non-const because the shallow 
   * copy provides write access to the underlying shared array.
   *
   * The lifetime of an association with a device array should not be
   * allowed to extend beyond the function in which the association is
   * created. This is necessary to guarantee that no such association
   * remains when the device array is deallocated or destroyed. Such an
   * association may be released by calling the dissociate member function 
   * of the host array, or will be released by the host array destructor.
   * A reference counting system is used to detect and report erroneous 
   * de-allocation of a device array when it is still referred to by 
   * one or more other containers.
   *
   * <b> Usage </b>: 
   *
   * Typical usage is shown below for host-to-device transfer to a long
   * lived instance of DeviceArray<Data,CPT> named dArray from a shorter
   * lived instance of HostArray<Data,CPT> named hArray, where Data
   * denotes the type of each array element:
   * \code
   *    HostArray<Data,CPT> hArray;
   *    hArray.associate(dArray);
   *
   *    \\ ( Initialize data in hArray )
   *
   *    dArray = hArray
   *    hArray.dissociate();
   * \endcode
   *
   * Comments:
   *
   *   - In this specialization for a CPU backend (T=CPT), the 
   *     assignment (=) operator that assigns a RHS HostArray<Data,CPT> 
   *     to a LHS DeviceArray<Data,CPT> template normally does nothing.
   *     It is provided for compatability with the syntax required for
   *     GPU code, with T=CUT, for which the analogous assignment operator
   *     actually copies data from the CPU host memory to GPU device 
   *     global memory. 
   * 
   *   - In this specialization for a CPU backend (T=CPT), the dissociate
   *     function destroys the shallow copy owned by the host array, by 
   *     nullifying a pointer held by the host array. In the corresponding 
   *     specialization for a Cuda GPU backend, with T=CUT, the analogous 
   *     function does nothing.
   *
   *   - The host array should never be used to modify data after the 
   *     assignment operator and before the dissociate function is
   *     invoked. Doing so would modify data owned by the device array
   *     in CPU code (T=CPT) but would have no effect on data owned by
   *     the device array in GPU code (T=CUT), causing inconsistent
   *     behavior. To enforce this, it is good practice to invoke the 
   *     dissociate member function immediately after host-to-device 
   *     assignment, as shown above.
   * 
   *   - The invocation of dissociate may sometimes be omitted when
   *     the host array is a local variable that is defined within 
   *     a function and the assignment operator occurs near the end
   *     of that function, because any remaining association with a
   *     device array will be destroyed by the host array destructor 
   *     function when the host array goes out of scope.
   *
   * \ingroup Pscf_Backend_Cpp_Module
   */
   template <typename Data>
   class HostArray<Data, CPT> final : public Array<Data> 
   {

   public:

      /**
      * Backend identifier class typename alias.
      */
      using BackendIdClass = CPT;

      // Default constructor.
      HostArray() = default;

      // Copy construction (delete).
      HostArray(HostArray<Data,CPT> const &) = delete;

      /**
      * Destructor.
      */
      ~HostArray();

      // Assignment (delete).
      HostArray<Data,CPT>& 
      operator = (HostArray<Data,CPT> const&) = delete;

      /**
      * Associate this object with a device array.
      *
      * Associates this object with a device array, thus creating a
      * shallow copy that points to the same memory. After successful
      * return, isAllocated() will return true.
      *
      * \throw Exception if source array is not allocated on entry
      * \throw Exception if this array is already associated
      *
      * \param other  device array that owns the data
      */
      void associate(DeviceArray<Data,CPT> & other);

      /**
      * Release association with a device array.
      *
      * Upon successful return, isAllocated() will return false.
      *
      * \throw Exception if this is not associated with another array.
      */
      void dissociate();

   private:

      /// Reference to a device array that has an array referenced by this.
      CountedReference ref_;

      // Hide inherited protected members by making them private
      using Array<Data>::data_;
      using Array<Data>::capacity_;

   };

   // Explicit instantiation declarations
   extern template class HostArray<double,CPT>;
   extern template class HostArray<fftw_complex,CPT>;

} // namespace Pscf

#include <pscf/backend/cpp/DeviceArray.h>
#include <util/global.h>

namespace Pscf {

   // Member functions

   /*
   * Destructor.
   */
   template <typename Data>
   HostArray<Data,CPT>::~HostArray()
   {
      if (ref_.isAssociated()) {
         ref_.dissociate(); // decrements counter of source array
      }
   }

   /*
   * Associate this object with a device array.
   */
   template <typename Data>
   void
   HostArray<Data,CPT>::associate(DeviceArray<Data,CPT> & source)
   {
      UTIL_CHECK(source.isAllocated());
      UTIL_CHECK(!ref_.isAssociated());

      // Copy data pointer and size
      data_ = source.cArray();
      capacity_ = source.capacity();

      // Associate ReferenceCounter of the source array with this
      ref_.associate(source);

      // On exit, the ReferenceCounter sub-object of the data source is
      // incremented and the ref_ CountedReference member variable of
      // this object holds a pointer to that ReferenceCounter.
   }

   /*
   * Release pre-existing association with a device array.
   */
   template <typename Data>
   void HostArray<Data,CPT>::dissociate()
   {
      UTIL_CHECK(data_);
      UTIL_CHECK(ref_.isAssociated());

      data_ = nullptr;
      capacity_ = 0;
      ref_.dissociate(); // decrements reference counter of source array
   }

}
#endif

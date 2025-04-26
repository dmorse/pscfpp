#ifndef PRDC_CPU_FFTW_D_ARRAY_H
#define PRDC_CPU_FFTW_D_ARRAY_H

/*
* PSCF Package 
*
* Copyright 2015 - 2025, The Regents of the University of Minnesota
* Distributed under the terms of the GNU General Public License.
*/

#include <util/containers/Array.h>
#include <util/global.h>

#include <fftw3.h>

namespace Pscf {
namespace Prdc {
namespace Cpu {

   using namespace Util;

   /**
   * Dynamic array with data aligned for use with FFTW library.
   *
   * The allocate and deallocate functions of this class use functions
   * provided by the FFTW library to allocate and free aligned memory. 
   * The class is otherwise similar in most respects to a Util::DArray.
   *
   * \ingroup Prdc_Cpu_Module
   */
   template <typename Data>
   class FftwDArray : public Array<Data>
   {

   public:

      /**
      * Default constructor.
      */
      FftwDArray();

      /**
      * Destructor.
      *
      * Deletes underlying C array, if allocated previously.
      */
      virtual ~FftwDArray();

      /**
      * Allocate the underlying C array.
      *
      * \throw Exception if the FftwDArray is already allocated.
      *
      * \param capacity number of elements to allocate.
      */
      void allocate(int capacity);

      /**
      * Dellocate the underlying C array.
      *
      * \throw Exception if the FftwDArray is not allocated.
      */
      virtual void deallocate();

      /**
      * Return true if the FftwDArray has been allocated, false otherwise.
      */
      bool isAllocated() const;

      /**
      * Serialize a FftwDArray to/from an Archive.
      *
      * \param ar       archive
      * \param version  archive version id
      */
      template <class Archive>
      void serialize(Archive& ar, const unsigned int version);

   protected:

       using Array<Data>:: data_;
       using Array<Data>:: capacity_;

   private:

      /**
      * Copy constructor (private and not implemented to prohibit).
      */
      FftwDArray(FftwDArray<Data> const & other);

      /**
      * Assignment operator (private and non implemented to prohibit).
      */
      FftwDArray<Data>& operator = (FftwDArray<Data> const & other);

   };

   /*
   * Return true if the FftwDArray has been allocated, false otherwise.
   */
   template <typename Data>
   inline bool FftwDArray<Data>::isAllocated() const
   {  return (bool) data_; }

   /*
   * Serialize a FftwDArray to/from an Archive.
   */
   template <typename Data>
   template <class Archive>
   void FftwDArray<Data>::serialize(Archive& ar, const unsigned int version)
   {
      int capacity;
      if (Archive::is_saving()) {
         capacity = capacity_;
      }
      ar & capacity;
      if (Archive::is_loading()) {
         if (!isAllocated()) {
            if (capacity > 0) {
               allocate(capacity);
            }
         } else {
            if (capacity != capacity_) {
               UTIL_THROW("Inconsistent FftwDArray capacities");
            }
         }
      }
      if (isAllocated()) {
         for (int i = 0; i < capacity_; ++i) {
            ar & data_[i];
         }
      }
   }

}
}
}
#include "FftwDArray.tpp"
#endif

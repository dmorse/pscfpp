#ifndef PRDC_CPU_C_FIELD_H
#define PRDC_CPU_C_FIELD_H

/*
* PSCF Package 
*
* Copyright 2016 - 2022, The Regents of the University of Minnesota
* Distributed under the terms of the GNU General Public License.
*/

#include <prdc/cpu/FftwDArray.h>
#include <prdc/cpu/complex.h>
#include <pscf/math/IntVec.h>
#include <util/global.h>

namespace Pscf {
namespace Prdc {
namespace Cpu {

   using namespace Util;
   using namespace Pscf;

   /**
   * Field of complex double precision values on an FFT mesh.
   * 
   * \ingroup Prdc_Cpu_Module 
   */
   template <int D>
   class CField : public FftwDArray<fftw_complex>
   {

   public:

      // Typedefs

      /**
      * Type of each element.
      */
      typedef fftw_complex ElementType;

      /**
      * Complex number type.
      */
      typedef fftw_complex Complex;

      /**
      * Real and imaginary parts of a Complex number.
      */
      typedef double Real;

      // Member functions

      /**
      * Default constructor.
      */
      CField();

      /**
      * Copy constructor.
      *
      * Allocates new memory and copies all elements by value.
      *
      *\param other the CField to be copied.
      */
      CField(const CField& other);

      /**
      * Destructor.
      *
      * Deletes underlying C array, if allocated previously.
      */
      virtual ~CField();

      /**
      * Assignment operator.
      *
      * If this Field is not allocated, allocates and copies all elements.
      *
      * If this and the other Field are both allocated, the capacities must
      * be exactly equal. If so, this method copies all elements.
      *
      * \param other the RHS CField
      */
      CField& operator = (const CField& other);

      /**
      * Allocate the underlying C array for an FFT grid.
      *
      * \throw Exception if the CField is already allocated.
      *
      * \param meshDimensions vector of numbers of grid points per direction
      */
      void allocate(const IntVec<D>& meshDimensions);

      /**
      * Deallocate underlying C array and clear mesh dimensions.
      */
      virtual void deallocate();

      /**
      * Return mesh dimensions by constant reference.
      */
      const IntVec<D>& meshDimensions() const;

      /**
      * Serialize a Field to/from an Archive.
      *
      * \param ar       archive
      * \param version  archive version id
      */
      template <class Archive>
      void serialize(Archive& ar, const unsigned int version);

   private:

      // Vector containing number of grid points in each direction.
      IntVec<D> meshDimensions_;

      using FftwDArray<fftw_complex>::allocate;

   };

   /*
   * Return mesh dimensions by constant reference.
   */
   template <int D>
   inline const IntVec<D>& CField<D>::meshDimensions() const
   {  return meshDimensions_; }

   /*
   * Serialize a Field to/from an Archive.
   */
   template <int D>
   template <class Archive>
   void CField<D>::serialize(Archive& ar, const unsigned int version)
   {
      FftwDArray<fftw_complex>::serialize(ar, version);
      ar & meshDimensions_;
   }

   #ifndef PRDC_CPU_C_FIELD_TPP
   extern template class CField<1>;
   extern template class CField<2>;
   extern template class CField<3>;
   #endif

}
}
}
#endif

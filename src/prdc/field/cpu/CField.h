#ifndef PRDC_CPU_C_FIELD_H
#define PRDC_CPU_C_FIELD_H

/*
* PSCF - Polymer Self-Consistent Field 
*
* Copyright 2015 - 2025, The Regents of the University of Minnesota
* Distributed under the terms of the GNU General Public License.
*/

#include <pscf/cpu/CppTp.h>         // backend 
#include <pscf/cpu/FftwDRArray.h>   // base class
#include <pscf/math/IntVec.h>       // member

#include <fftw3.h>

namespace Pscf {
namespace Prdc {

   using namespace Util;

   // Declare primary template
   template <int D, class T> class CField;

   /**
   * Field of complex double precision values on a mesh (CPU).
   * 
   * \ingroup Prdc_Cpu_Module 
   */
   template <int D>
   class CField<D, CppTp<D> > 
    : public FftwDRArray<fftw_complex>
   {

   public:

      // Type aliases

      using FftwDRArray<fftw_complex>::ValueType;

      /**
      * Type of real and imaginary parts of a complex value
      */
      using RealType = double;

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
      * If this field is not allocated, this function allocates the field
      * and copies all elements.
      *
      * If this and the other field are both allocated, the capacities must
      * be equal. If so, this functions copies all elements.
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

      using FftwDRArray<fftw_complex>::allocate;

   };

   /*
   * Return mesh dimensions by constant reference.
   */
   template <int D> inline 
   IntVec<D> const & CField<D, CppTp<D> >::meshDimensions() const
   {  return meshDimensions_; }

   /*
   * Serialize a Field to/from an Archive.
   */
   template <int D>
   template <class Archive>
   void CField<D, CppTp<D> >::serialize(
                        Archive& ar, 
                        const unsigned int version)
   {
      FftwDRArray<fftw_complex>::serialize(ar, version);
      ar & meshDimensions_;
   }

   // Explicit instantiation declarations
   extern template class CField<1, CppTp<1> >;
   extern template class CField<2, CppTp<2> >;
   extern template class CField<3, CppTp<3> >;

} // namespace Prdc
} // namespace Pscf
#endif

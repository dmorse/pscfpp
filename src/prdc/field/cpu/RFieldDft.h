#ifndef PRDC_CPU_R_FIELD_DFT_H
#define PRDC_CPU_R_FIELD_DFT_H

/*
* PSCF - Polymer Self-Consistent Field 
*
* Copyright 2015 - 2025, The Regents of the University of Minnesota
* Distributed under the terms of the GNU General Public License.
*/

#include <pscf/cpu/CppTp.h>         // class template argument
#include <pscf/cpu/FftwDRArray.h>   // base class
#include <pscf/math/IntVec.h>       // member

#include <fftw3.h>

namespace Pscf {
namespace Prdc {

   using namespace Util;

   // Declaration of primary template
   template <int D, class T> class RFieldDft;

   /**
   * Fourier transform of a real field on an FFT mesh (CPU version).
   *
   * \ingroup Prdc_Cpu_Module
   */
   template <int D>
   class RFieldDft<D, CppTp<D> >
    : public FftwDRArray<fftw_complex>
   {

   public:

      // Type aliases

      /**
      * Complex type of an array element.
      */
      using FftwDRArray<fftw_complex>::ValueType;

      /**
      * Type of real and imaginary parts of a complex element value.
      */
      using RealType = double;

      // Member functions

      /**
      * Default constructor.
      */
      RFieldDft();

      /**
      * Copy constructor.
      *
      * Allocates new memory and copies all elements by value.
      *
      *\param other the RFieldDft to be copied.
      */
      RFieldDft(RFieldDft<D, CppTp<D> > const & other);

      /**
      * Destructor.
      *
      * Deletes underlying C array, if allocated previously.
      */
      virtual ~RFieldDft();

      /**
      * Assignment operator.
      *
      * If this Field is not allocated, allocates and copies all elements.
      *
      * If this and the other Field are both allocated, the capacities must
      * be exactly equal. If so, this method copies all elements.
      *
      * \param other the RHS Field
      */
      RFieldDft<D, CppTp<D> >& 
      operator = (RFieldDft<D, CppTp<D> > const & other);

      /**
      * Allocate the underlying C array and set mesh dimensions.
      *
      * \throw Exception if the RFieldDft is already allocated.
      *
      * \param meshDimensions vector of grid points in each direction
      */
      void allocate(IntVec<D> const & meshDimensions);

      /**
      * Deallocate underlying C array and clear mesh dimensions.
      */
      virtual void deallocate();

      /**
      * Return vector of spatial mesh dimensions by const reference.
      */
      IntVec<D> const & meshDimensions() const;

      /**
      * Return vector of dft (Fourier) grid dimensions by const reference.
      *  
      * The last element of dftDimensions() and meshDimensions() differ by
      * about a factor of two: dftDimension()[D-1] = meshDimensions()/2 + 1.
      * For D > 1, other elements are equal. 
      */
      IntVec<D> const & dftDimensions() const;

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

      // Vector containing dimensions of dft (Fourier) grid.
      IntVec<D> dftDimensions_;

      using FftwDRArray<fftw_complex>::allocate;

   };

   /*
   * Return mesh dimensions by const reference.
   */
   template <int D> inline 
   IntVec<D> const & RFieldDft<D, CppTp<D> >::meshDimensions() const
   {  return meshDimensions_; }

   /*
   * Return dimensions of dft grid by const reference.
   */
   template <int D> inline 
   IntVec<D> const & RFieldDft<D, CppTp<D> >::dftDimensions() const
   {  return dftDimensions_; }

   /*
   * Serialize a Field to/from an Archive.
   */
   template <int D>
   template <class Archive>
   void RFieldDft<D, CppTp<D> >::serialize(Archive& ar, 
                                           const unsigned int version)
   {
      FftwDRArray<fftw_complex>::serialize(ar, version);
      ar & meshDimensions_;
      ar & dftDimensions_;
   }

   // Explicit instantiation declarations
   extern template class RFieldDft<1, CppTp<1> >;
   extern template class RFieldDft<2, CppTp<2> >;
   extern template class RFieldDft<3, CppTp<3> >;

} // namespace Prdc
} // namespace Pscf
#endif

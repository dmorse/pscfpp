#ifndef PRDC_CPU_R_FIELD_H
#define PRDC_CPU_R_FIELD_H

/*
* PSCF Package 
*
* Copyright 2016 - 2022, The Regents of the University of Minnesota
* Distributed under the terms of the GNU General Public License.
*/

#include <prdc/cpu/Field.h>
#include <pscf/math/IntVec.h>
#include <util/global.h>

namespace Pscf {
namespace Prdc {
namespace Cpu {

   using namespace Util;
   using namespace Pscf;

   /**
   * Field of real double precision values on an FFT mesh.
   * 
   * \ingroup Prdc_Cpu_Module 
   */
   template <int D>
   class RField : public Field<double>
   {

   public:

      /**
      * Default constructor.
      */
      RField();

      /**
      * Copy constructor.
      *
      * Allocates new memory and copies all elements by value.
      *
      *\param other the RField to be copied.
      */
      RField(const RField& other);

      /**
      * Destructor.
      *
      * Deletes underlying C array, if allocated previously.
      */
      virtual ~RField();

      /**
      * Assignment operator.
      *
      * If this Field is not allocated, allocates and copies all elements.
      *
      * If this and the other Field are both allocated, the capacities must
      * be exactly equal. If so, this method copies all elements.
      *
      * \param other the RHS RField
      */
      RField& operator = (const RField& other);

      /**
      * Allocate the underlying C array for an FFT grid.
      *
      * \throw Exception if the RField is already allocated.
      *
      * \param meshDimensions vector of numbers of grid points per direction
      */
      void allocate(IntVec<D> const & meshDimensions);

      /**
      * Deallocate memory and return to empty state.
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

      using Field<double>::allocate;

   };

   /*
   * Return mesh dimensions by constant reference.
   */
   template <int D>
   inline const IntVec<D>& RField<D>::meshDimensions() const
   {  return meshDimensions_; }

   /*
   * Serialize a Field to/from an Archive.
   */
   template <int D>
   template <class Archive>
   void RField<D>::serialize(Archive& ar, const unsigned int version)
   {
      Field<double>::serialize(ar, version);
      ar & meshDimensions_;
   }

   #ifndef PRDC_CPU_R_FIELD_TPP
   extern template class RField<1>;
   extern template class RField<2>;
   extern template class RField<3>;
   #endif

}
}
}
#endif

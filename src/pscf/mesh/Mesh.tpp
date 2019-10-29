#ifndef PSCF_MESH_TPP
#define PSCF_MESH_TPP

/*
* PSCF - Polymer Self-Consistent Field Theory
*
* Copyright 2016 - 2019, The Regents of the University of Minnesota
* Distributed under the terms of the GNU General Public License.
*/

#include "Mesh.h"
#include <util/format/Int.h>
#include <util/global.h>

namespace Pscf
{

   using namespace Util;

   template <int D>
   Mesh<D>::Mesh()
    : dimensions_(0),
      offsets_(0),
      size_(0)
   {
      IntVec<D> dimensions(1);
      setDimensions(dimensions);
   }

   template <int D>
   Mesh<D>::Mesh(const IntVec<D>& dimensions)
    : dimensions_(0),
      offsets_(0),
      size_(0)
   {
      setDimensions(dimensions);
   }

   template <int D>
   void Mesh<D>::setDimensions(const IntVec<D>& dimensions)
   {
      int i;
      for (i = 0; i < D; ++i) {
         if (dimensions[i] <= 0) {
            UTIL_THROW("Mesh dimensions must be positive");
         }
      }

      dimensions_ = dimensions;
      offsets_[D -1] = 1;
      for (i = D - 1; i > 0; --i) {
         offsets_[i-1] = offsets_[i]*dimensions_[i];
      }
      size_ = offsets_[0]*dimensions_[0];
   }

   template <int D>
   int Mesh<D>::rank(const IntVec<D>& position) const
   {
      int i;
      int result = 0;
      for (i = 0; i < D - 1; ++i) {
         assert(position[i] >= 0);
         assert(position[i] < dimensions_[i]);
         result += position[i]*offsets_[i];
      }
      assert(position[i] >= 0);
      assert(position[i] < dimensions_[i]);
      result += position[i];
      return result;
   }

   template <int D>
   IntVec<D> Mesh<D>::position(int id) const
   {
      IntVec<D> position;
      int       remainder = id;

      int i;
      for (i = 0; i < D - 1; ++i) {
         position[i] = remainder/offsets_[i];
         remainder -= position[i]*offsets_[i];
      }
      position[i] = remainder;
      return position;
   }

   template <int D>
   bool Mesh<D>::isInMesh(int coordinate, int i) const
   {
      bool result = true;
      if (coordinate <  0)
         result = false;
      if (coordinate >= dimensions_[i])
         result = false;
      return result;
   }

   template <int D>
   bool Mesh<D>::isInMesh(IntVec<D>& position) const
   {
      bool result = true;
      for (int i = 0; i < D; ++i) {
         if (position[i] <  0)
            result = false;
         if (position[i] >= dimensions_[i])
            result = false;
      }
      return result;
   }

   template <int D>
   int Mesh<D>::shift(int& coordinate, int i) const
   {
      int shift;
      if (coordinate >= 0) {
         shift = coordinate/dimensions_[i];
      } else {
         shift = -1 + ((coordinate+1)/dimensions_[i]);
      }
      coordinate -= shift*dimensions_[i];
      return shift;
   }

   template <int D>
   IntVec<D> Mesh<D>::shift(IntVec<D>& position) const
   {
      IntVec<D> shifts;
      for (int i = 0; i < D; ++i) {
         shifts[i] = shift(position[i], i);
      }
      return shifts;
   }

   template <int D>
   std::istream& operator >> (std::istream& in, 
                              Mesh<D>& mesh)
   {
      IntVec<D> dimensions;
      in >> dimensions;
      for (int i = 0; i < D; ++i) {
         UTIL_CHECK(dimensions[i] > 0);
      }
      mesh.setDimensions(dimensions);
      return in;
   }

   template <int D>
   std::ostream& operator << (std::ostream& out,
                              Mesh<D>& mesh)
   {
      for (int i = 0; i < D; ++i) {
         out << " " << Int(mesh.dimensions_[i], 6);
      }
      return out;
   }

}
#endif

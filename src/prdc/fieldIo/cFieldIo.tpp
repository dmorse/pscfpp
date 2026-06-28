#ifndef PRDC_CFIELD_IO_TPP
#define PRDC_CFIELD_IO_TPP

/*
* PSCF - Polymer Self-Consistent Field
*
* Copyright 2015 - 2025, The Regents of the University of Minnesota
* Distributed under the terms of the GNU General Public License.
*/

#include <pscf/mesh/MeshIterator.h>
#include <pscf/math/IntVec.h>
#include <pscf/math/arithmetic.h>

#include <util/containers/DArray.h>
#include <util/format/Dbl.h>

namespace Pscf {
namespace Prdc {

   template <int D, class AT>
   void readCFieldsData(std::istream& in, 
                        DArray< AT > & fields,
                        IntVec<D> const& dimensions)
   {
      double x, y;
      MeshIterator<D> iter(dimensions);
      int nMonomer = fields.capacity();
      UTIL_CHECK(nMonomer > 0);
      int rank, j;
      for (j = 0; j < nMonomer; ++j) {
         for (iter.begin(); !iter.atEnd(); ++iter) {
            rank = iter.rank();
            in >> x;
            UTIL_ASSERT(in.good());
            in >> y;
            UTIL_ASSERT(in.good());
            assign(fields[j][rank], x, y);
         }
      }
      UTIL_CHECK(in.good());
   }

   template <int D, class AT>
   void readCFieldData(std::istream& in, 
                      AT& field,
                      IntVec<D> const& dimensions)
   {
      double x, y;
      MeshIterator<D> iter(dimensions);
      int rank;
      for (iter.begin(); !iter.atEnd(); ++iter) {
         rank = iter.rank();
         in >> x;
         UTIL_ASSERT(in.good());
         in >> y;
         UTIL_ASSERT(in.good());
         assign(field[rank], x, y);
      }
      UTIL_CHECK(in.good());
   }

   template <int D, class AT, typename CT, typename RT>
   void writeCFieldsData(std::ostream& out,
                       DArray< AT > const& fields,
                       IntVec<D> const& dimensions)
   {
      RT x, y;
      MeshIterator<D> iter(dimensions);
      int nMonomer = fields.capacity();
      UTIL_CHECK(nMonomer > 0);
      int rank, j;
      for (j = 0; j < nMonomer; ++j) {
         for (iter.begin(); !iter.atEnd(); ++iter) {
            rank = iter.rank();
            x = real( fields[j][rank] );
            y = imag( fields[j][rank] );
            out << " " << Dbl(x, 21, 13)
                << " " << Dbl(y, 21, 13);
         }
         out << std::endl;
      }
   }

   template <int D, class AT, typename CT, typename RT>
   void writeCFieldData(std::ostream& out,
                       AT const& field,
                       IntVec<D> const& dimensions)
   {
      RT x, y;
      MeshIterator<D> iter(dimensions);
      int rank;
      for (iter.begin(); !iter.atEnd(); ++iter) {
         rank = iter.rank();
         x = real( field[rank] );
         y = imag( field[rank] );
         out << " " << Dbl(x, 21, 13)
             << " " << Dbl(y, 21, 13)
             << std::endl;
      }
   }

} // namespace Prdc
} // namespace Pscf
#endif

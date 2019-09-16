#ifndef PSCF_SPACE_GROUP_H
#define PSCF_SPACE_GROUP_H

/*
* Simpatico - Simulation Package for Polymeric and Molecular Liquids
*
* Copyright 2016 - 2019, The Regents of the University of Minnesota
* Distributed under the terms of the GNU General Public License.
*/

#include <pscf/crystal/SpaceSymmetry.h>
#include <pscf/crystal/SymmetryGroup.h>
#include <pscf/math/IntVec.h>
#include <util/containers/FSArray.h>
#include <util/param/Label.h>
#include <iostream>

namespace Pscf
{

   using namespace Util;

   /**
   * Crystallographic space group.
   *
   * \ingroup Pscf_Crystal_Module
   */
   template <int D>
   class SpaceGroup : public SymmetryGroup< SpaceSymmetry <D> >
   {};

   // Template function definition

   /**
   * Output stream inserter operator for a SpaceGroup<D>.
   *
   * \param out  output stream
   * \param g  space group
   * \ingroup Pscf_Crystal_Module
   */ 
   template <int D>
   std::ostream& operator << (std::ostream& out, const SpaceGroup<D>& g)
   {
      int size = g.size();
      out << "dim  " << D << std::endl;
      out << "size " << size << std::endl;
      for (int i = 0; i < size; ++i) {
         out << std::endl;
         out << g[i];
      }
      return out;
   }

   /**
   * Input stream extractor operator for a SpaceGroup<D>.
   *
   * \param in  input stream
   * \param g  space group
   *
   * \ingroup Pscf_Crystal_Module
   */ 
   template <int D>
   std::istream& operator >> (std::istream& in, SpaceGroup<D>& g)
   {
      int dim, size;
      in >> Label("dim") >> dim;
      UTIL_CHECK(D == dim);
      in >> Label("size") >> size;

      SpaceSymmetry<D> s;
      g.clear();
      for (int i = 0; i < size; ++i) {
         in >> s;
         g.add(s);
      }
      return in;
   }

   #ifndef PSCF_SPACE_GROUP_CPP
   extern template class SpaceGroup<1>;
   extern template class SpaceGroup<2>;
   extern template class SpaceGroup<3>;
   #endif

}
#endif

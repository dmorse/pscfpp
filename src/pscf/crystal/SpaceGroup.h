#ifndef PSCF_SPACE_GROUP_H
#define PSCF_SPACE_GROUP_H

/*
* Simpatico - Simulation Package for Polymeric and Molecular Liquids
*
* Copyright 2010 - 2017, The Regents of the University of Minnesota
* Distributed under the terms of the GNU General Public License.
*/

#include <pscf/crystal/SpaceSymmetry.h>
#include <pscf/crystal/SymmetryGroup.h>
#include <pscf/math/IntVec.h>
#include <util/containers/FSArray.h>
#include <iostream>

namespace Pscf
{

   using namespace Util;

   /**
   * Group of crystal symmetries with no translations.
   *
   * \ingroup Pscf_Crystal_Module
   */
   template <int D>
   class SpaceGroup : public SymmetryGroup< SpaceSymmetry <D> >
   {};

   // Template function definition

   /*
   * Output stream inserter operator for a SpaceGroup<D>.
   */ 
   template <int D>
   std::ostream& operator << (std::ostream& out, const SpaceGroup<D>& g)
   {
      int i, size;
      size = g.size();
      out << size << std::endl;
      for (i = 0; i < size; ++i) {
         out << std::endl;
         out << g[i];
      }
      return out;
   }

   /*
   * Output stream inserter operator for a SpaceGroup<D>.
   */ 
   template <int D>
   std::istream& operator >> (std::istream& in, SpaceGroup<D>& g)
   {
      SpaceSymmetry<D> s;
      int size;
      in >> size;
      for (int i = 0; i < size; ++i) {
         in >> s;
         g.add(s);
      }
      return in;
   }
}
#endif

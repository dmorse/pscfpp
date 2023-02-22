#ifndef PSCF_SPACE_GROUP_H
#define PSCF_SPACE_GROUP_H

/*
* Simpatico - Simulation Package for Polymeric and Molecular Liquids
*
* Copyright 2016 - 2022, The Regents of the University of Minnesota
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
   class SpaceGroup : public SymmetryGroup< SpaceSymmetry<D> >
   {

   public:

      /**
      * Determines if a symmetric structure has an inversion center.
      *
      * Returns true if an inversion center exists, and false otherwise.
      * If an inversion center exists, its location is returned as the
      * output value of output argument "center". 
      *
      * \param center location of inversion center, if any (output)
      */
      bool 
      hasInversionCenter(typename SpaceSymmetry<D>::Translation& center) 
      const;

      /**
      * Shift the origin of space used in the coordinate system.
      *
      * This function modifies each symmetry elements in the group so as
      * to refer to an equivalent symmetry defined using a new coordinate 
      * system with a shifted origin. The argument gives the coordinates 
      * of the origin of the new coordinate system as defined in the old
      * coordinate system.
      *
      * \param origin  location of origin of the new coordinate system
      */
      void 
      shiftOrigin(typename SpaceSymmetry<D>::Translation const & origin);

      // Using declarations for some inherited functions
      using SymmetryGroup< SpaceSymmetry <D> >::size;

   };

   // Template function definition

   /**
   * Output stream inserter operator for a SpaceGroup<D>.
   *
   * \param out  output stream
   * \param g  space group
   * \ingroup Pscf_Crystal_Module
   */ 
   template <int D>
   std::ostream& operator << (std::ostream& out, SpaceGroup<D> const & g)
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

   /**
   * Open and read a group file.
   * 
   * \param groupName  name of group, or group file (input)
   * \param group  space group (output)
   *
   * \ingroup Pscf_Crystal_Module
   */
   template <int D>
   void readGroup(std::string groupName, SpaceGroup<D>& group);

   /**
   * Open and write a group file.
   * 
   * \param filename  output file name
   * \param group  space group 
   *
   * \ingroup Pscf_Crystal_Module
   */
   template <int D>
   void writeGroup(std::string filename, SpaceGroup<D> const & group);

   #ifndef PSCF_SPACE_GROUP_TPP
   extern template class SpaceGroup<1>;
   extern template class SpaceGroup<2>;
   extern template class SpaceGroup<3>;
   extern template void readGroup(std::string, SpaceGroup<1>& );
   extern template void readGroup(std::string, SpaceGroup<2>& );
   extern template void readGroup(std::string, SpaceGroup<3>& );
   extern template void writeGroup(std::string, SpaceGroup<1> const &);
   extern template void writeGroup(std::string, SpaceGroup<2> const &);
   extern template void writeGroup(std::string, SpaceGroup<3> const &);
   #endif

}
#endif

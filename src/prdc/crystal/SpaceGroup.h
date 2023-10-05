#ifndef PSCF_SPACE_GROUP_H
#define PSCF_SPACE_GROUP_H

/*
* Simpatico - Simulation Package for Polymeric and Molecular Liquids
*
* Copyright 2016 - 2022, The Regents of the University of Minnesota
* Distributed under the terms of the GNU General Public License.
*/

#include <prdc/crystal/SpaceSymmetry.h>
#include <prdc/crystal/SymmetryGroup.h>
#include <pscf/math/IntVec.h>
#include <util/containers/FSArray.h>
#include <util/param/Label.h>
#include <iostream>

namespace Pscf {
namespace Prdc {

   using namespace Util;

   /**
   * Crystallographic space group.
   *
   * \ingroup Pscf_Prdc_Crystal_Module
   */
   template <int D>
   class SpaceGroup : public SymmetryGroup< SpaceSymmetry<D> >
   {

   public:

      /**
      * Determines if this space group has an inversion center.
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

      /**
      * Check if input mesh dimensions are compatible with space group.
      *
      * This function checks if a mesh with the specified dimensions is
      * invariant under all operations of this space group, i.e., whether
      * each crystal symmetry operation maps the position of every node 
      * of the mesh onto the position of another node. It is only possible
      * define how a symmetry operation transforms a function that is 
      * defined only on the nodes of mesh if the mesh is invariant under 
      * the symmetry operation, in this sense. An invariant mesh must thus 
      * be used necessary to describe a function whose values on the mesh 
      * nodes are invariant under all operations in the space group. 
      *
      * If the mesh is not invariant under all operations of the space
      * group, an explanatory error message is printed and an Exception 
      * is thrown to halt execution.
      *
      * The mesh for a unit cell within a Bravais lattice is assumed to 
      * be a regular orthogonal mesh in a space of reduced coordinates, 
      * which are the components of position defined using a Bravais
      * basis (i.e., a basis of Bravais lattice basis vectors). Each
      * element of the dimensions vector is equal to the number of grid 
      * points along a direction corresponding to a Bravais lattice vector.
      * A Bravais basis is also used to define elements of the matrix 
      * representation of the point group operation and the translation
      * vector in the representation of a crystal symmetry operation as
      * an instance of class Pscf::SpaceSymmetry<D>. 
      *
      * \param dimensions vector of mesh dimensions
      */
      void checkMeshDimensions(IntVec<D> const & dimensions) const;

      // Using declarations for some inherited functions
      using SymmetryGroup< SpaceSymmetry <D> >::size;

   };

   // Template function definition

   /**
   * Output stream inserter operator for a SpaceGroup<D>.
   *
   * \param out  output stream
   * \param g  space group
   * \ingroup Pscf_Prdc_Crystal_Module
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
   * \ingroup Pscf_Prdc_Crystal_Module
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
   * \ingroup Pscf_Prdc_Crystal_Module
   */
   template <int D>
   void readGroup(std::string groupName, SpaceGroup<D>& group);

   /**
   * Open and write a group file.
   * 
   * \param filename  output file name
   * \param group  space group 
   *
   * \ingroup Pscf_Prdc_Crystal_Module
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
}
#endif

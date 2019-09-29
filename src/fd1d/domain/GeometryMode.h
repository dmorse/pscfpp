#ifndef FD1D_GEOMETRY_MODE_H
#define FD1D_GEOMETRY_MODE_H

/*
* Simpatico - Simulation Package for Polymeric and Molecular Liquids
*
* Copyright 2010 - 2014, The Regents of the University of Minnesota
* Distributed under the terms of the GNU General Public License.
*/

#include <util/archives/serialize.h>
#include <iostream>

namespace Pscf{
namespace Fd1d
{

   /**
   * Enumeration of geometrical modes for functions of one coordinate.
   *
   * Allowed values are: Planar, Cylindrical Spherical.
   *
   * \ingroup Fd1d_Domain_Module
   */
   enum GeometryMode {Planar, Cylindrical, Spherical};

   /**
   * istream extractor for a GeometryMode.
   *
   * \param  in       input stream
   * \param  lattice  GeometryMode to be read
   * \return modified input stream
   */
   std::istream& operator >> (std::istream& in, GeometryMode& lattice);

   /**
   * ostream inserter for an GeometryMode.
   *
   * \param  out      output stream
   * \param  lattice  GeometryMode to be written
   * \return modified output stream
   */
   std::ostream& operator << (std::ostream& out, GeometryMode lattice);

   /**
   * Serialize a GeometryMode value.
   *
   * \param ar      archive object
   * \param lattice value to be serialized
   * \param version archive version id
   */
   template <class Archive>
   void serialize(Archive& ar, GeometryMode& lattice, const unsigned int version)
   {  serializeEnum(ar, lattice, version); }

}
}
#endif

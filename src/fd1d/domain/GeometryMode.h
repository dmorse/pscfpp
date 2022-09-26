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
   * GeometryMode values can be read from or written to iostream using
   * overloaded extractor (>>) and inserter (<<) operators. The text
   * representations of the three values are "planar", "cylindrical" and
   * "spherical".
   *
   * \ingroup Fd1d_Domain_Module
   */
   enum GeometryMode {Planar, Cylindrical, Spherical};

   /**
   * istream extractor for a GeometryMode.
   *
   * \param  in       input stream
   * \param  mode  GeometryMode to be read
   * \return modified input stream
   */
   std::istream& operator >> (std::istream& in, GeometryMode& mode);

   /**
   * ostream inserter for an GeometryMode.
   *
   * \param  out      output stream
   * \param  mode  GeometryMode to be written
   * \return modified output stream
   */
   std::ostream& operator << (std::ostream& out, GeometryMode mode);

   /**
   * Serialize a GeometryMode value.
   *
   * \param ar      archive object
   * \param mode value to be serialized
   * \param version archive version id
   */
   template <class Archive>
   void serialize(Archive& ar, GeometryMode& mode, const unsigned int version)
   {  serializeEnum(ar, mode, version); }

}
}
#endif

#ifndef PRDC_FIELD_HEADER_H
#define PRDC_FIELD_HEADER_H

/*
* PSCF - Polymer Self-Consistent Field Theory
*
* Copyright 2015 - 2025, The Regents of the University of Minnesota
* Distributed under the terms of the GNU General Public License.
*/

#include "UnitCell.h"
#include <iostream>
#include <iomanip>

namespace Pscf {
namespace Prdc {

   using namespace Util;

   /**
   * Read common part of field header (fortran PSCF format).
   *
   * If the group_name label and value are absent, then this function
   * returns an empty groupName string, groupName == "". 
   *
   * \param in  input stream
   * \param ver1  major file format version number (output)
   * \param ver2  major file format version number (output)
   * \param cell  UnitCell<D> object (output)
   * \param groupName  string identifier for space group (output)
   * \param nMonomer  number of monomers (output)
   * \ingroup Prdc_Crystal_Module
   */
   template <int D>
   void readFieldHeader(std::istream& in, int& ver1, int& ver2, 
                        UnitCell<D>& cell, std::string& groupName,
                        int& nMonomer);

   /**
   * Write common part of field header (fortran PSCF format).
   *
   * If the groupName parameter is an empty string (groupName == ""), 
   * then this function does not write the group_name label or value. 
   *
   * If nMonomer == 0, then this function does not write the N_monomer
   * label or value. 
   *
   * \param out  output stream
   * \param ver1  major file format version number (input)
   * \param ver2  major file format version number (input)
   * \param cell  UnitCell<D> object (input)
   * \param groupName  string identifier for space group (input)
   * \param nMonomer  number of monomers (input)
   * \ingroup Prdc_Crystal_Module
   */
   template <int D>
   void writeFieldHeader(std::ostream &out, int ver1, int ver2,
                         UnitCell<D> const & cell,
                         std::string const & groupName,
                         int nMonomer);

}
}
#include "fieldHeader.tpp"
#endif

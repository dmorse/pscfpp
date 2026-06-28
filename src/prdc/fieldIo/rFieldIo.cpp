/*
* PSCF - Polymer Self-Consistent Field
*
* Copyright 2015 - 2025, The Regents of the University of Minnesota
* Distributed under the terms of the GNU General Public License.
*/

#include "rFieldIo.h"

namespace Pscf {
namespace Prdc {

   /*
   * Read the number of basis functions from a field file header.
   */
   int readNBasis(std::istream& in)
   {
   
      // Read the label, which can be N_star or N_basis
      std::string label;
      in >> label;
      UTIL_ASSERT(in.good());
      if (label != "N_star" && label != "N_basis") {
         std::string msg =  "\n";
         msg += "Error reading field file:\n";
         msg += "Expected N_basis or N_star, but found [";
         msg += label;
         msg += "]";
         UTIL_THROW(msg.c_str());
      }
 
      // Read the value of nBasis
      int nBasis;
      in >> nBasis;
      UTIL_CHECK(in.good());
      UTIL_CHECK(nBasis > 0);

      return nBasis;
   }

   /*
   * Write the number of basis functions to a field file header.
   */
   void writeNBasis(std::ostream& out, int nBasis)
   {
      out << "N_basis      " << std::endl
          << "             " << nBasis << std::endl;
   }

} // namespace Prdc
} // namespace Pscf

/*
* PFTS - Polymer Field Theory Simulator
*
* Copyright 2013, David Morse (morse012@.umn.edu)
* Distributed under the terms of the GNU General Public License.
*/

#include "Monomer.h"

namespace Chem{ 

   Monomer::Monomer()
    : id_(-1),
      step_(0.0),
      name_()
   {}

   /* 
   * Extract a Monomer from an istream.
   */
   std::istream& operator >> (std::istream& in, Monomer& monomer)
   {
      in >> monomer.id_;
      in >> monomer.name_;
      in >> monomer.step_;
      return in;
   }
   
   /* 
   * Output a Monomer to an ostream, without line breaks.
   */
   std::ostream& operator << (std::ostream& out, const Monomer& monomer) 
   {
      out << monomer.id_;
      out << "  " << monomer.name_;
      out << "  ";
      out.setf(std::ios::scientific);
      out.width(15);
      out.precision(8);
      out << monomer.step_;
      return out;
   }

} 

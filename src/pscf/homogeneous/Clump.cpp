/*
* PSCF - Polymer Self-Consistent Field Theory
*
* Copyright 2016 - 2019, The Regents of the University of Minnesota
* Distributed under the terms of the GNU General Public License.
*/

#include "Clump.h"

namespace Pscf { 
namespace Homogeneous { 

   /*
   * Constructor.
   */ 
   Clump::Clump()
    : monomerId_(-1),
      size_(0.0)
   {}
  
   /*
   * Set the monomer id.
   */ 
   void Clump::setMonomerId(int monomerId)
   {  monomerId_ = monomerId; }
  
   /*
   * Set the size of this block.
   */ 
   void Clump::setSize(double size)
   {  size_ = size; }
  
   /* 
   * Extract a Clump from an istream.
   */
   std::istream& operator>>(std::istream& in, Clump &block)
   {
      in >> block.monomerId_;
      in >> block.size_;
      return in;
   }
   
   /* 
   * Output a Clump to an ostream, without line breaks.
   */
   std::ostream& operator<<(std::ostream& out, const Clump &block) 
   {
      out << "  " << block.monomerId_;
      out << "  ";
      out.setf(std::ios::scientific);
      out.width(16);
      out.precision(8);
      out << block.size_;
      return out;
   }

} 
} 

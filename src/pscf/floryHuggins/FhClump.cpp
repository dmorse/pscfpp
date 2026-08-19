/*
* PSCF - Polymer Self-Consistent Field
*
* Copyright 2015 - 2026, The Regents of the University of Minnesota
* Distributed under the terms of the GNU General Public License.
*/

#include "FhClump.h"

namespace Pscf { 

   /*
   * Constructor.
   */ 
   FhClump::FhClump()
    : monomerId_(-1),
      size_(0.0)
   {}
  
   /*
   * Set the monomer id.
   */ 
   void FhClump::setMonomerId(int monomerId)
   {  monomerId_ = monomerId; }
  
   /*
   * Set the size of this block.
   */ 
   void FhClump::setSize(double size)
   {  size_ = size; }
  
   /* 
   * Extract a FhClump from an istream.
   */
   std::istream& operator>>(std::istream& in, FhClump &block)
   {
      in >> block.monomerId_;
      in >> block.size_;
      return in;
   }
   
   /* 
   * Output a FhClump to an ostream, without line breaks.
   */
   std::ostream& operator<<(std::ostream& out, const FhClump &block) 
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

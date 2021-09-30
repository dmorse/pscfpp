/*
* PSCF - Polymer Self-Consistent Field Theory
*
* Copyright 2016 - 2019, The Regents of the University of Minnesota
* Distributed under the terms of the GNU General Public License.
*/

#include "SolventDescriptor.h"

namespace Pscf
{ 

   /*
   * Constructor.
   */ 
   SolventDescriptor::SolventDescriptor()
    : monomerId_(-1),
      size_()
   {}

   SolventDescriptor::~SolventDescriptor()
   {}

   /*
   * Set the id for this solvent.
   */ 
   void SolventDescriptor::setMonomerId(int monomerId)
   {  monomerId_ = monomerId; }
  
   /*
   * Set the id for this solvent.
   */ 
   void SolventDescriptor::setSize(double size)
   {  size_ = size; }
  
   /* 
   * Extract a SolventDescriptor from an istream.
   */
   std::istream& operator >> (std::istream& in, SolventDescriptor &solvent)
   {
      in >> solvent.monomerId_;
      in >> solvent.size_;
      return in;
   }
   
   /* 
   * Output a SolventDescriptor to an ostream, without line breaks.
   */
   std::ostream& 
   operator<<(std::ostream& out, const SolventDescriptor &solvent) 
   {
      out << solvent.monomerId_;
      out << "  ";
      out.setf(std::ios::scientific);
      out.width(16);
      out.precision(8);
      out << solvent.size_;
      return out;
   }

} 

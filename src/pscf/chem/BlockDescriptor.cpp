/*
* PSCF - Polymer Self-Consistent Field Theory
*
* Copyright 2016 - 2019, The Regents of the University of Minnesota
* Distributed under the terms of the GNU General Public License.
*/

#include "BlockDescriptor.h"

namespace Pscf
{ 

   /*
   * Constructor.
   */ 
   BlockDescriptor::BlockDescriptor()
    : id_(-1),
      monomerId_(-1),
      length_(-1.0),
      vertexIds_(),
      polymerType_(PolymerType::Branched)
   {}

   /*
   * Set the id for this block.
   */ 
   void BlockDescriptor::setId(int id)
   {  id_ = id; }
  
   /*
   * Set indices of associated vertices.
   */ 
   void BlockDescriptor::setVertexIds(int vertexId0, int vertexId1)
   { 
      vertexIds_[0] = vertexId0; 
      vertexIds_[1] = vertexId1; 
   }
  
   /*
   * Set the monomer id.
   */ 
   void BlockDescriptor::setMonomerId(int monomerId)
   {  monomerId_ = monomerId; }
  
   /*
   * Set the length of this block.
   */ 
   void BlockDescriptor::setLength(double length)
   {  length_ = length; }
  
   /*
   * Set the type of the polymer containing this block.
   */ 
   void BlockDescriptor::setPolymerType(PolymerType::Enum type)
   {  polymerType_ = type; }
  
   /* 
   * Extract a BlockDescriptor from an istream.
   */
   std::istream& operator>>(std::istream& in, BlockDescriptor &block)
   {
      in >> block.monomerId_;
      in >> block.length_;
      if (block.polymerType_ == PolymerType::Branched) {
         in >> block.vertexIds_[0];
         in >> block.vertexIds_[1];
      }
      return in;
   }
   
   /* 
   * Output a BlockDescriptor to an ostream, without line breaks.
   */
   std::ostream& operator<<(std::ostream& out, const BlockDescriptor &block) 
   {
      out << "  " << block.monomerId_;
      out << "  ";
      out.setf(std::ios::scientific);
      out.width(20);
      out.precision(12);
      out << block.length_;
      if (block.polymerType_ == PolymerType::Branched) {
         out << "  " << block.vertexIds_[0];
         out << "  " << block.vertexIds_[1];
      }
      return out;
   }

} 

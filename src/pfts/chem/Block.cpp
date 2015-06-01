/*
* PFTS - Polymer Field Theory Simulator
*
* Copyright 2013, David Morse (morse012@.umn.edu)
* Distributed under the terms of the GNU General Public License.
*/

#include "Block.h"

namespace Pfts{ 

   /*
   * Constructor.
   */ 
   Block::Block()
    : id_(-1),
      monomerId_(-1),
      vertexIds_(),
      length_(0.0)
   {}

   /*
   * Set the id for this block.
   */ 
   void Block::setId(int id)
   {  id_ = id; }
  
   /*
   * Set indices of associated vertices.
   */ 
   void Block::setVertexIds(int vertexId0, int vertexId1)
   { 
      vertexIds_[0] = vertexId0; 
      vertexIds_[1] = vertexId1; 
   }
  
   /*
   * Set the monomer id.
   */ 
   void Block::setMonomerId(int monomerId)
   {  monomerId_ = monomerId; }
  
   /*
   * Set the length of this block.
   */ 
   void Block::setLength(double length)
   {  length_ = length; }
  
   /* 
   * Extract a Block from an istream.
   */
   std::istream& operator>>(std::istream& in, Block &block)
   {
      in >> block.id_;
      in >> block.monomerId_;
      in >> block.vertexIds_[0];
      in >> block.vertexIds_[1];
      in >> block.length_;
      return in;
   }
   
   /* 
   * Output a Block to an ostream, without line breaks.
   */
   std::ostream& operator<<(std::ostream& out, const Block &block) 
   {
      // out.setf(std::ios::scientific);
      // out.width(Block::Width);
      // out.precision(Block::Precision);
      out << block.id_;
      out << block.monomerId_;
      out << block.vertexIds_[0];
      out << block.vertexIds_[1];
      out << block.length_;
      return out;
   }

} 

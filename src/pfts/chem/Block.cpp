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
   void Block::setVertexIds(int vertexAId, int vertexBId)
   { 
      vertexIds_[0] = vertexAId; 
      vertexIds_[1] = vertexBId; 
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
  
} 

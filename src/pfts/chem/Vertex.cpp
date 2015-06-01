/*
* PFTS - Polymer Field Theory Simulator
*
* Copyright 2013, David Morse (morse012@.umn.edu)
* Distributed under the terms of the GNU General Public License.
*/

#include "Vertex.h"
#include "Block.h"
#include <util/global.h>

namespace Pfts{ 

   Vertex::Vertex()
    : inSolverIds_(),
      outSolverIds_(),
      id_(-1)
   {}

   Vertex::~Vertex()
   {}

   void Vertex::setId(int id)
   {  id_ = id; }

   void Vertex::addBlock(const Block& block)
   {
      // Preconditions
      if (id_ < 0) {
         UTIL_THROW("Negative vertex id");
      }
      if (block.id() < 0) {
         UTIL_THROW("Negative block id");
      }
      if (block.vertexId(0) == block.vertexId(1)) {
         UTIL_THROW("Error: Equal vertex indices in block");
      }

      Pair<int> solverId;
      solverId[0] = block.id();
      if (block.vertexId(0) == id_) {
         solverId[1] = 0;
         outSolverIds_.append(solverId);
         solverId[1] = 1;
         inSolverIds_.append(solverId);
      } else
      if (block.vertexId(1) == id_) {
         solverId[1] = 1;
         outSolverIds_.append(solverId);
         solverId[1] = 0;
         inSolverIds_.append(solverId);
      } else {
         UTIL_THROW("Neither block vertex id matches this vertex");
      }
   }

} 

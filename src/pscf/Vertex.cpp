/*
* PFTS - Polymer Field Theory Simulator
*
* Copyright 2013, David Morse (morse012@.umn.edu)
* Distributed under the terms of the GNU General Public License.
*/

#include "Vertex.h"
#include "Block.h"
#include <util/global.h>

namespace Pscf
{ 

   Vertex::Vertex()
    : inPropagatorIds_(),
      outPropagatorIds_(),
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

      Pair<int> propagatorId;
      propagatorId[0] = block.id();
      if (block.vertexId(0) == id_) {
         propagatorId[1] = 0;
         outPropagatorIds_.append(propagatorId);
         propagatorId[1] = 1;
         inPropagatorIds_.append(propagatorId);
      } else
      if (block.vertexId(1) == id_) {
         propagatorId[1] = 1;
         outPropagatorIds_.append(propagatorId);
         propagatorId[1] = 0;
         inPropagatorIds_.append(propagatorId);
      } else {
         UTIL_THROW("Neither block vertex id matches this vertex");
      }
   }

} 

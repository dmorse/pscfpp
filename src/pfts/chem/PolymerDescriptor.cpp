/*
* PFTS - Polymer Field Theory Simulator
*
* Copyright 2013, David Morse (morse012@.umn.edu)
* Distributed under the terms of the GNU General Public License.
*/

#include "PolymerDescriptor.h"

namespace Pfts{ 

   /*
   * Constructor.
   */
   PolymerDescriptor::PolymerDescriptor()
    : blocks_(),
      vertices_(),
      plan_(),
      volume_(0.0),
      nBlock_(0),
      nVertex_()
   {
      setClassName("PolymerDescriptor");
   }

   void PolymerDescriptor::readParameters(std::istream& in)
   {
      read<int>(in, "nBlock", nBlock_);
      read<int>(in, "nVertex", nVertex_);
      blocks_.allocate(nBlock_);
      readDArray<Block>(in, "blocks", blocks_, nBlock_);

      // Allocate array of vertices and set vertex indices
      vertices_.allocate(nVertex_);
      for (int vertexId = 0; vertexId < nVertex_; ++vertexId) {
         vertices_[vertexId].setId(vertexId);
      }

      // Add blocks to vertices
      int vertexId0, vertexId1;
      Block* blockPtr;
      for (int blockId = 0; blockId < nBlock_; ++blockId) {
          blockPtr = &(blocks_[blockId]);
          vertexId0 = blockPtr->vertexId(0);
          vertexId1 = blockPtr->vertexId(1);
          vertices_[vertexId0].addBlock(*blockPtr);
          vertices_[vertexId1].addBlock(*blockPtr);
      }

   }

} 

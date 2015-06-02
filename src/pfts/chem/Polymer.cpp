/*
* PFTS - Polymer Field Theory Simulator
*
* Copyright 2013, David Morse (morse012@.umn.edu)
* Distributed under the terms of the GNU General Public License.
*/

#include "Polymer.h"

namespace Pfts{ 

   /*
   * Constructor.
   */
   Polymer::Polymer()
    : blocks_(),
      vertices_(),
      plan_(),
      volume_(0.0),
      nBlock_(0),
      nVertex_()
   {}

   void Polymer::readParameters(std::istream& in)
   {
      read<int>(in, "nBlock", nBlock_);
      read<int>(in, "nVertex", nVertex_);
      readDArray<Block>(in, "blocks", blocks_, nBlock_);

      // Allocate array of vertices and set vertex indices
      vertices_.allocate(nVertex_);
      for (int vertexId = 0; vertexId < nVertex_; ++vertexId) {
         vertices_[vertexId].setId(vertexId);
      }

      // Add blocks to vertices
      int vertexId0, vertexId1;
      for (int blockId = 0; blockId < nBlock_; ++blockId) {
          vertexId0 = blocks_[blockId].vertexId(0);
          vertexId1 = blocks_[blockId].vertexId(1);
          vertices_[vertexId0].addBlock(blocks_[blockId]);
          vertices_[vertexId1].addBlock(blocks_[blockId]);
      }

   }

} 

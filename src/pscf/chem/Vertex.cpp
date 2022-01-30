/*
* PSCF - Polymer Self-Consistent Field Theory
*
* Copyright 2016 - 2019, The Regents of the University of Minnesota
* Distributed under the terms of the GNU General Public License.
*/

#include "Vertex.h"
#include "BlockDescriptor.h"
#include <util/global.h>

namespace Pscf
{ 

   /*
   * Constructor.
   */
   Vertex::Vertex()
    : inPropagatorIds_(),
      outPropagatorIds_(),
      id_(-1)
   {}

   /*
   * Destructor.
   */
   Vertex::~Vertex()
   {}


   /*
   * Set integer id.
   */
   void Vertex::setId(int id)
   {  id_ = id; }

   /*
   * Add this block to the list.
   */
   void Vertex::addBlock(const BlockDescriptor& block)
   {
      // Preconditions
      if (id_ < 0) {
         UTIL_THROW("Negative vertex id");
      }
      if (block.id() < 0) {
         UTIL_THROW("Negative block id");
      }
      if (block.vertexId(0) < 0) {
         UTIL_THROW("Error: Negative block vertexId 0");
      }
      if (block.vertexId(1) < 0) {
         UTIL_THROW("Error: Negative block vertexId 1");
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

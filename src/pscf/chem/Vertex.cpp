/*
* PSCF - Polymer Self-Consistent Field Theory
*
* Copyright 2015 - 2025, The Regents of the University of Minnesota
* Distributed under the terms of the GNU General Public License.
*/

#include "Vertex.h"
#include "Edge.h"
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
   * Add a edge to the vertex.
   */
   void Vertex::addEdge(const Edge& edge)
   {
      // Preconditions
      if (id_ < 0) {
         UTIL_THROW("Negative vertex id");
      }
      if (edge.id() < 0) {
         UTIL_THROW("Negative edge id");
      }
      if (edge.vertexId(0) < 0) {
         UTIL_THROW("Error: Negative edge vertexId 0");
      }
      if (edge.vertexId(1) < 0) {
         UTIL_THROW("Error: Negative edge vertexId 1");
      }
      if (edge.vertexId(0) == edge.vertexId(1)) {
         UTIL_THROW("Error: Equal vertex indices in edge");
      }

      Pair<int> propagatorId;
      propagatorId[0] = edge.id();
      if (edge.vertexId(0) == id_) {
         propagatorId[1] = 0;
         outPropagatorIds_.append(propagatorId);
         propagatorId[1] = 1;
         inPropagatorIds_.append(propagatorId);
      } else
      if (edge.vertexId(1) == id_) {
         propagatorId[1] = 1;
         outPropagatorIds_.append(propagatorId);
         propagatorId[1] = 0;
         inPropagatorIds_.append(propagatorId);
      } else {
         UTIL_THROW("Neither edge vertex id matches this vertex");
      }
   }

} 

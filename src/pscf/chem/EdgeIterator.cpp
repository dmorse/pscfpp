/*
* PSCF - Polymer Self-Consistent Field Theory
*
* Copyright 2016 - 2025, The Regents of the University of Minnesota
* Distributed under the terms of the GNU General Public License.
*/

#include "EdgeIterator.h"
#include <pscf/chem/PolymerSpecies.h>
#include <pscf/chem/Edge.h>
#include <util/containers/Pair.h>

namespace Pscf {

   /*
   * Constructor
   */
   EdgeIterator::EdgeIterator(PolymerSpecies const & polymer)
    : currentVertexId_(-1),
      targetVertexId_(-1),
      polymerPtr_(&polymer)
   {}

   /*
   * Destructor
   */
   EdgeIterator::~EdgeIterator()
   {}

   /*
   * Initialize iterator.
   */
   void EdgeIterator::begin(int sourceId, int targetId)
   {
      UTIL_CHECK(sourceId != targetId);
      currentEdgeId_ = sourceId;
      targetEdgeId_  = targetId;

      // Identify ids for vertices of source edge
      int s0 = polymerPtr_->edge(sourceId).vertexId(0);
      int s1 = polymerPtr_->edge(sourceId).vertexId(1);

      // Identify ids for vertices of target edge
      int t0 = polymerPtr_->edge(targetId).vertexId(0);
      int t1 = polymerPtr_->edge(targetId).vertexId(1);

      // Paths from source vertices towards target vertices
      Pair<int> p00, p01, p10, p11;
      p00 = polymerPtr_->path(s0, t0);
      p01 = polymerPtr_->path(s0, t1);
      p10 = polymerPtr_->path(s1, t0);
      p11 = polymerPtr_->path(s1, t1);
      // Find paths that go through the source edge
      if (p00[0] == currentEdgeId_) {
         UTIL_CHECK(p01[0] == currentEdgeId_);
         UTIL_CHECK(p10[0] != currentEdgeId_);
         UTIL_CHECK(p11[0] != currentEdgeId_);
         currentVertexId_ = s1;
      } else 
      if (p10[0] == currentEdgeId_) {
         UTIL_CHECK(p11[0] == currentEdgeId_);
         UTIL_CHECK(p00[0] != currentEdgeId_);
         UTIL_CHECK(p01[0] != currentEdgeId_);
         currentVertexId_ = s0;
      } else {
         UTIL_THROW("Error in finding leading vertex for source");
      }

      // Paths from target vertices towards source vertices
      p00 = polymerPtr_->path(t0, s0);
      p01 = polymerPtr_->path(t0, s1);
      p10 = polymerPtr_->path(t1, s0);
      p11 = polymerPtr_->path(t1, s1);
      // Find paths that go through the target edge
      if (p00[0] == targetEdgeId_) {
         UTIL_CHECK(p01[0] == targetEdgeId_);
         UTIL_CHECK(p10[0] != targetEdgeId_);
         UTIL_CHECK(p11[0] != targetEdgeId_);
         targetVertexId_ = t0;
      } else 
      if (p10[0] == targetEdgeId_) {
         UTIL_CHECK(p11[0] == targetEdgeId_);
         UTIL_CHECK(p00[0] != targetEdgeId_);
         UTIL_CHECK(p01[0] != targetEdgeId_);
         targetVertexId_ = t1;
      } else {
         UTIL_THROW("Error in finding target vertex");
      }

   }

   /*
   * Increment vertex - update current vertex to next one in path.
   */
   EdgeIterator& EdgeIterator::operator ++ ()
   {
      UTIL_CHECK(notEnd());
      Pair<int> propId = polymerPtr_->path(currentVertexId_, targetVertexId_);
      int edgeId = propId[0];
      int dirId = propId[1];
      UTIL_CHECK(edgeId >= 0);
      UTIL_CHECK(edgeId < polymerPtr_->nBlock());
      UTIL_CHECK(dirId >= 0);
      UTIL_CHECK(dirId < 2);
      Edge const & edge = polymerPtr_->edge(edgeId);
      UTIL_CHECK(edge.vertexId(dirId) == currentVertexId_);
      if (dirId == 0) {
         currentEdgeId_ = edgeId;
         currentVertexId_ = edge.vertexId(1);
      } else {
         currentEdgeId_ = edgeId;
         currentVertexId_ = edge.vertexId(0);
      }
      return *this;
   }

   /*
   * Get the current edge id.
   */ 
   int EdgeIterator::currentEdgeId() const
   {  return currentEdgeId_; }

   /*
   * Get the current vertex id.
   */ 
   int EdgeIterator::currentVertexId() const
   {  return currentVertexId_; }

   /*
   * Is the current vertex equal to the target:
   */
   bool EdgeIterator::isEnd() const
   {
      return (bool)( currentVertexId_ == targetVertexId_ );
   }

   /*
   * Is the current vertex not equal to the target.
   */
   bool EdgeIterator::notEnd() const
   {
      return (bool)( currentVertexId_ != targetVertexId_ );
   }

}

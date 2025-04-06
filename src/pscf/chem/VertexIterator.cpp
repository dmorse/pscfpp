/*
* PSCF - Polymer Self-Consistent Field Theory
*
* Copyright 2016 - 2025, The Regents of the University of Minnesota
* Distributed under the terms of the GNU General Public License.
*/

#include "VertexIterator.h"
#include <pscf/chem/PolymerSpecies.h>
#include <pscf/chem/Edge.h>

namespace Pscf {

   /*
   * Constructor
   */
   VertexIterator::VertexIterator(PolymerSpecies const & polymer)
    : currentId_(-1),
      targetId_(-1),
      polymerPtr_(&polymer)
   {}

   /*
   * Destructor
   */
   VertexIterator::~VertexIterator()
   {}

   /*
   * Initialize iterator.
   */
   void VertexIterator::begin(int sourceId, int targetId)
   {
      currentId_ = sourceId;
      targetId_  = targetId;
   }

   /*
   * Increment vertex - update current vertex to next one in path.
   */
   VertexIterator& VertexIterator::operator ++ ()
   {
      UTIL_CHECK(notEnd());
      Pair<int> propId = polymerPtr_->path(currentId_, targetId_);
      int edgeId = propId[0];
      int dirId = propId[1];
      UTIL_CHECK(edgeId >= 0);
      UTIL_CHECK(edgeId < polymerPtr_->nBlock());
      UTIL_CHECK(dirId >= 0);
      UTIL_CHECK(dirId < 2);
      Edge const & edge = polymerPtr_->edge(edgeId);
      if (dirId == 0) {
         currentId_ = edge.vertexId(1);
      } else {
         currentId_ = edge.vertexId(0);
      }
      return *this;
   }

   /*
   * Get the current vertex id.
   */ 
   int VertexIterator::currentId() const
   {  return currentId_; }

   /*
   * Is the current vertex equal to the target:
   */
   bool VertexIterator::isEnd() const
   {
      return (bool)( currentId_ == targetId_ );
   }

   /*
   * Is the current vertex not equal to the target.
   */
   bool VertexIterator::notEnd() const
   {
      return (bool)( currentId_ != targetId_ );
   }

}

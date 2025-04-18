#ifndef PSCF_EDGE_ITERATOR_H
#define PSCF_EDGE_ITERATOR_H

/*
* PSCF - Polymer Self-Consistent Field Theory
*
* Copyright 2016 - 2022, The Regents of the University of Minnesota
* Distributed under the terms of the GNU General Public License.
*/

#include <pscf/chem/PolymerSpecies.h>   

namespace Pscf { 

   // Forward reference
   class PolymerSpecies;

   using namespace Util;

   /**
   * Edge iterator for graph associated with a polymer.
   *
   * Usage: Suppose that object p is an instance of a subclass of 
   * PolymerSpecies. The following snippet illustrates how to iterate 
   * from edge (or block) iSource to edge iTarget of the associated 
   * polymer:
   * \code
   *    EdgeIterator = iter(p);
   *    int iEdge;
   *    for (iter.begin(iSource, iTarget); iter.notEnd(); ++iter) {
   *       iEdge = iter.currentEdgeId();
   *       // do something with edge iEdge
   *    }
   * \endcode
   *
   * \ingroup Pscf_Chem_Module
   */
   class EdgeIterator 
   {

   public:

      /**
      * Constructor.
      *
      * \param polymer  associated PolymerSpecies object
      */
      EdgeIterator(PolymerSpecies const & polymer);
   
      /**
      * Destructor.
      */
      ~EdgeIterator();

      /**
      * Initialize iterator. 
      *
      * \param sourceId  index of the initial (or source) edge
      * \param targetId  index of the final (or target) edge
      */
      void begin(int sourceId, int targetId);

      /**
      * Increment operator - move to next vertex.
      */
      EdgeIterator& operator ++ ();

      /**
      * Get index of the current edge.
      */
      int currentEdgeId() const;

      /**
      * Get direction index for the path within the current edge.
      */
      int currentDirectionId() const;

      /**
      * Get index of the current vertex.
      *
      * When the current edge is not the target edge, the current vertex
      * is the vertex of the current edge that is closer to the closest
      * vertex of the target edge.
      *
      * When the current edge is also the target edge, the current vertex
      * is the vertex of the target edge that is farther from the initial
      * source edge.
      */
      int currentVertexId() const;

      /**
      * Return true iff currentId == targetId.
      */
      bool isEnd() const;

      /**
      * Return true iff currentId != targetId.
      */
      bool notEnd() const;

   private:

      // Index of current edge.
      int currentEdgeId_;

      // Direction index for the current edge.
      int currentDirectionId_;

      // Index of current vertex.
      int currentVertexId_;

      // Index of target edge.
      int targetEdgeId_;

      // Index of target vertex.
      int targetVertexId_;

      // Pointer to associated PolymerSpecies object.
      PolymerSpecies const * polymerPtr_;

   };

}
#endif 

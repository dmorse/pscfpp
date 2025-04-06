#ifndef PSCF_VERTEX_ITERATOR_H
#define PSCF_VERTEX_ITERATOR_H

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
   * Vertex iterator for graph associated with a polymer.
   *
   * \ingroup Pscf_Chem_Module
   */
   class VertexIterator 
   {

   public:

      /**
      * Constructor.
      *
      * \param polymer  associated PolymerSpecies object
      */
      VertexIterator(PolymerSpecies const & polymer);
   
      /**
      * Destructor.
      */
      ~VertexIterator();

      /**
      * Initialize iterator. 
      */
      void begin(int sourceId, int targetId);

      /**
      * Increment operator - move to next vertex.
      */
      VertexIterator& operator ++ ();

      /**
      * Get index of the current vertex.
      */
      int currentId() const;

      /**
      * Return true iff currentId == targetId.
      */
      bool isEnd() const;

      /**
      * Return true iff currentId != targetId.
      */
      bool notEnd() const;

   private:

      // Index of current vertex.
      int currentId_;

      // Index of target vertex.
      int targetId_;

      // Pointer to associated PolymerSpecies object.
      PolymerSpecies const * polymerPtr_;

   };

}
#endif 

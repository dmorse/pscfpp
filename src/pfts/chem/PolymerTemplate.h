#ifndef PFTS_TY_H
#define PFTS_TY_H

/*
* PFTS - Polymer Field Theory Simulator
*
* Copyright 2013, David Morse (morse012@.umn.edu)
* Distributed under the terms of the GNU General Public License.
*/

namespace Pfts{

   using namespace Util;

   template <class TPropagator, class TCField>
   class PolymerTemplate : public Species, public PolymerDescriptor
   {
   public:
   
      /**
      * Read parameters and initialize.
      */
      virtual void readParameters();
   
      /**
      * Get monomer concentration field for specific block.
      *
      * \param blockId integer index of associated block
      */
      const TCField& blockCField(int blockId) const;
   
      /**
      * Get propagator for a specific block and direction.
      *
      * Convention: The directionId for a propagagor is based on the
      * indexing of vertex ids for the associated Block, as returned 
      * by Block::vertexId(int). Let p be a PolymerTemplate object,
      * let i be a block index, and let v0 = p.block.vertexId(0) and 
      * v1 = p.block(i).vertexId(1) be indices of associated vertices.
      * Then, p.propagator(i,v0) propagates vertex v0 to v1, while
      * p.propagator(i, v1) propagates from vertex v1 to v0.
      *
      * \param blockId integer index of associated block
      * \param directionId integer index for direction (0 or 1)
      */
      const TPropagator& propagator(int blockId, int directionId) const;
   
   protected:

      /**
      * Create connections among propagators.
      */
      void init();
   
      /**
      * Array of propagators, indexed by block and direction.
      */ 
      DArray<Pair<TPropagator>> propagators_;
   
      /**
      * Array of block concentration fields, indexed by block.
      */ 
      DArray<TCField> blockCFields_;
   
   };

} 
#endif 

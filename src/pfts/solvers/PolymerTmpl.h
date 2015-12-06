#ifndef PFTS_POLYMER_TMPL_H
#define PFTS_POLYMER_TMPL_H

/*
* PFTS - Polymer Field Theory Simulator
*
* Copyright 2013, David Morse (morse012@.umn.edu)
* Distributed under the terms of the GNU General Public License.
*/

#include <pfts/solvers/Species.h>
#include <pfts/chem/PolymerDescriptor.h>
#include <util/containers/DArray.h>
#include <util/containers/DMatrix.h>

namespace Pfts
{

   using namespace Util;

   /**
   * Base class template for classes that represent polymer species.
   */
   template <class TProp, class TCField>
   class PolymerTmpl : public PolymerDescriptor, public Species
   {
   public:
   
      /**
      * Read parameters and initialize.
      */
      virtual void readParameters(std::istream& in);
   
      /**
      * Get the monomer concentration field for a specific block.
      *
      * \param blockId integer index of associated block
      */
      TCField& blockCField(int blockId);
   
      /**
      * Get propagator for a specific block and direction.
      *
      * Suppose p is a PolymerTmpl, and b = p.block(i) is a block 
      * terminating at vertices with indices v0 = b.vertexId(0) and 
      * v1 = b.vertexId(1), then p.propagator(i, 0) is a propagator 
      * that from v0 to v1, while p.propagator(i, 1) propagates
      * from v1 to v0.
      *
      * \param blockId integer index of associated block
      * \param directionId integer index for direction (0 or 1)
      */
      TProp& propagator(int blockId, int directionId);
   
   private:

      /**
      * Array of propagators, indexed by block and direction.
      */ 
      DMatrix<TProp> propagators_;
   
      /**
      * Array of block concentration fields, indexed by block.
      */ 
      DArray<TCField> blockCFields_;
   
   };

   /*
   * Read parameters and allocate arrays.
   */
   template <class TProp, class TCField>
   void 
   PolymerTmpl<TProp, TCField>::readParameters(std::istream& in)
   {
      PolymerDescriptor::readParameters(in);
      blockCFields_.allocate(nBlock());
      propagators_.allocate(nBlock(), 2);

      // Add sources to all propagators
      Vertex const * vertexPtr = 0;
      TProp const * sourcePtr = 0;
      TProp* propagatorPtr = 0;
      Pair<int> propagatorId;
      int blockId, directionId, vertexId, i;
      for (blockId = 0; blockId < nBlock(); ++blockId) {
         for (directionId = 0; directionId < 2; ++directionId) {
            vertexId = block(blockId).vertexId(directionId);
            vertexPtr = &vertex(vertexId);
            propagatorPtr = &propagator(blockId, directionId);
            for (i = 0; i < vertexPtr->size(); ++i) {
               propagatorId = vertexPtr->inSolverId(i);
               if (propagatorId[0] != blockId) {
                   sourcePtr = & propagator(propagatorId[0], propagatorId[1]);
                   propagatorPtr->addSource(*sourcePtr);
               }
            }
         }
      }


   }

   /*
   * Get a block concentration.
   */
   template <class TProp, class TCField>
   TCField& 
   PolymerTmpl<TProp, TCField>::blockCField(int blockId)
   {  return blockCFields_[blockId]; }

   /*
   * Get a propagator.
   */
   template <class TProp, class TCField>
   TProp& 
   PolymerTmpl<TProp, TCField>::propagator(int blockId, int directionId)
   {  return propagators_(blockId, directionId); }

} 
#endif

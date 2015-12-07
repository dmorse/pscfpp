#ifndef PFTS_POLYMER_TMPL_H
#define PFTS_POLYMER_TMPL_H

/*
* PFTS - Polymer Field Theory Simulator
*
* Copyright 2013, David Morse (morse012@.umn.edu)
* Distributed under the terms of the GNU General Public License.
*/

#include <pfts/Species.h>
#include <pfts/PolymerDescriptor.h>
#include <util/containers/DArray.h>
#include <util/containers/DMatrix.h>

namespace Pfts
{

   using namespace Util;

   /**
   * Base class template for classes that represent polymer species.
   */
   template <class TP, class TC>
   class PolymerTmpl : public PolymerDescriptor, public Species
   {
   public:
 
      // Modified diffusion equation propagator for one block.
      typedef TP Propagator;

      // Monomer concentration field.
      typedef TC CField;
 
      /**
      * Read parameters and initialize.
      */
      virtual void readParameters(std::istream& in);
   
      /**
      * Get the monomer concentration field for a specific block.
      *
      * \param blockId integer index of associated block
      */
      CField& blockCField(int blockId);
   
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
      Propagator& propagator(int blockId, int directionId);
   
      /**
      * Get propagator indexed in order of computation.
      *
      * \param id integer index, in order of computation plan
      */
      Propagator& propagator(int id);
   
   private:

      /**
      * Array of propagators, indexed by block and direction.
      */ 
      DMatrix<Propagator> propagators_;
   
      /**
      * Array of block concentration fields, indexed by block.
      */ 
      DArray<CField> blockCFields_;
   
   };

   /*
   * Read parameters and allocate arrays.
   */
   template <class TP, class TC>
   void PolymerTmpl<TP, TC>::readParameters(std::istream& in)
   {
      PolymerDescriptor::readParameters(in);
      blockCFields_.allocate(nBlock());
      propagators_.allocate(nBlock(), 2);

      // Associate propagators with blocks and directions
      Propagator* propagatorPtr = 0;
      int blockId, directionId;
      for (blockId = 0; blockId < nBlock(); ++blockId) {
         for (directionId = 0; directionId < 2; ++directionId) {
            propagatorPtr = &propagator(blockId, directionId);
            propagatorPtr->setBlock(block(blockId), directionId);
         }
      }
      
      // Add sources to all propagators
      Vertex const * vertexPtr = 0;
      Propagator const * sourcePtr = 0;
      Pair<int> propagatorId;
      int vertexId, i;
      for (blockId = 0; blockId < nBlock(); ++blockId) {
         for (directionId = 0; directionId < 2; ++directionId) {
            vertexId = block(blockId).vertexId(directionId);
            vertexPtr = &vertex(vertexId);
            propagatorPtr = &propagator(blockId, directionId);
            for (i = 0; i < vertexPtr->size(); ++i) {
               propagatorId = vertexPtr->inPropagatorId(i);
               if (propagatorId[0] != blockId) {
                   sourcePtr = & propagator(propagatorId[0], propagatorId[1]);
                   propagatorPtr->addSource(*sourcePtr);
               }
            }
         }
      }

   }

   /*
   * Get a block monomer concentration.
   */
   template <class TP, class TC>
   TC& PolymerTmpl<TP, TC>::blockCField(int blockId)
   {  return blockCFields_[blockId]; }

   /*
   * Get a propagator indexed by block and direction.
   */
   template <class TP, class TC>
   TP& PolymerTmpl<TP, TC>::propagator(int blockId, int directionId)
   {  return propagators_(blockId, directionId); }

   /*
   * Get a propagator indexed in order of computation.
   */
   template <class TP, class TC>
   TP& PolymerTmpl<TP, TC>::propagator(int id)
   {
      UTIL_CHECK(id >= 0);
      UTIL_CHECK(id < nPropagator());
      Pair<int> propId = propagatorId(id);
      UTIL_CHECK(propId[0] >= 0);
      UTIL_CHECK(propId[0] < nBlock());
      UTIL_CHECK(propId[1] >= 0);
      UTIL_CHECK(propId[1] < 2);
      return propagators_(propId[0], propId[1]); 
   }

}
#endif

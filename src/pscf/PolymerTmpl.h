#ifndef PSCF_POLYMER_TMPL_H
#define PSCF_POLYMER_TMPL_H

/*
* PFTS - Polymer Field Theory Simulator
*
* Copyright 2013, David Morse (morse012@.umn.edu)
* Distributed under the terms of the GNU General Public License.
*/

#include <pscf/Species.h>
#include <pscf/PolymerDescriptor.h>
#include <util/containers/DArray.h>
#include <util/containers/DMatrix.h>

#include <cmath>

namespace Pscf
{

   using namespace Util;

   /**
   * Base class template for classes that represent polymer species.
   */
   template <class TProp>
   class PolymerTmpl : public PolymerDescriptor, public Species
   {
   public:
 
      // Modified diffusion equation propagator for one block.
      typedef TProp Propagator;

      // Monomer concentration field.
      typedef typename TProp::CField CField;
 
      // Chemical potential field.
      typedef typename TProp::WField WField;

      /**
      * Constructor.
      */
      PolymerTmpl();
 
      /**
      * Destructor.
      */
      ~PolymerTmpl();

      /**
      * Read parameters and initialize.
      */
      virtual void readParameters(std::istream& in);
 
      /**
      * Compute solution to modified diffusion equation.
      */ 
      virtual void compute(const DArray<WField>& wFields);
 
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

   // Inline functions

   /*
   * Get a block monomer concentration.
   */
   template <class TProp>
   inline
   typename PolymerTmpl<TProp>::CField& PolymerTmpl<TProp>::blockCField(int blockId)
   {  return blockCFields_[blockId]; }

   /*
   * Get a propagator indexed by block and direction.
   */
   template <class TProp>
   inline TProp& PolymerTmpl<TProp>::propagator(int blockId, int directionId)
   {  return propagators_(blockId, directionId); }

   /*
   * Get a propagator indexed in order of computation.
   */
   template <class TProp>
   inline TProp& PolymerTmpl<TProp>::propagator(int id)
   {
      Pair<int> propId = propagatorId(id);
      return propagators_(propId[0], propId[1]); 
   }

   // Non-inline functions

   /*
   * Constructor.
   */
   template <class TProp>
   PolymerTmpl<TProp>::PolymerTmpl()
   {}

   /*
   * Destructor.
   */
   template <class TProp>
   PolymerTmpl<TProp>::~PolymerTmpl()
   {}

   /*
   * Read parameters and allocate arrays.
   */
   template <class TProp>
   void PolymerTmpl<TProp>::readParameters(std::istream& in)
   {
      // Read polymer structure
      PolymerDescriptor::readParameters(in);

      // Read ensemble and phi or mu
      ensemble_ = Species::Closed;
      readOptional<Species::Ensemble>(in, "ensemble", ensemble_);
      if (ensemble_ == Species::Closed) {
         read(in, "phi", phi_);
      } else {
         read(in, "mu", mu_);
      }

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
      
      // Set sources and partners for all propagators
      Vertex const * vertexPtr = 0;
      Propagator const * sourcePtr = 0;
      Pair<int> propagatorId;
      int vertexId, i;
      for (blockId = 0; blockId < nBlock(); ++blockId) {
         // Add sources
         for (directionId = 0; directionId < 2; ++directionId) {
            vertexId = block(blockId).vertexId(directionId);
            vertexPtr = &vertex(vertexId);
            propagatorPtr = &propagator(blockId, directionId);
            for (i = 0; i < vertexPtr->size(); ++i) {
               propagatorId = vertexPtr->inPropagatorId(i);
               if (propagatorId[0] == blockId) {
                  UTIL_CHECK(propagatorId[1] != directionId);
               } else {
                  sourcePtr = & propagator(propagatorId[0], propagatorId[1]);
                  propagatorPtr->addSource(*sourcePtr);
               }
            }
         }
         // Set partners
         propagator(blockId, 0).setPartner(propagator(blockId, 1));
         propagator(blockId, 1).setPartner(propagator(blockId, 0));
      }
   }

   /*
   * Compute solution to modified diffusion equation, concentrations, etc.
   */ 
   template <class TProp>
   void PolymerTmpl<TProp>::compute(const DArray<WField>& wFields)
   {

      // Clear all propagators
      for (int j = 0; j < nPropagator(); ++j) {
         propagator(j).setIsSolved(false);
      }

      // Solve modified diffusion equation for all propagators
      int monomerId;
      for (int j = 0; j < nPropagator(); ++j) {
         if (!propagator(j).isReady()) {
            UTIL_THROW("Propagator not ready");
         }
         monomerId = propagator(j).block().monomerId();
         propagator(j).solve(wFields[monomerId]);
      }

      // Compute molecular partition function
      double q = propagator(0,0).computeQ();
      if (ensemble() == Species::Closed) {
         mu_ = log(phi_/q);
      } 
      else if (ensemble() == Species::Open) {
         phi_ = exp(mu_)*q;
      }

      // Compute block concentration fields
      double prefactor = phi_ / (q*length());
      for (int i = 0; i < nBlock(); ++i) {
         propagator(i, 0).computeConcentration(prefactor, blockCFields_[i]);
      }

   }
 
}
#endif

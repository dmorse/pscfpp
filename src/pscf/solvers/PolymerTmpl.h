#ifndef PSCF_POLYMER_TMPL_H
#define PSCF_POLYMER_TMPL_H

/*
* PSCF - Polymer Self-Consistent Field Theory
*
* Copyright 2016 - 2022, The Regents of the University of Minnesota
* Distributed under the terms of the GNU General Public License.
*/

#include <pscf/chem/Species.h>           // base class
#include <util/param/ParamComposite.h>   // base class

#include <pscf/chem/Monomer.h>           // member template argument
#include <pscf/chem/Vertex.h>            // member template argument
#include <pscf/chem/PolymerType.h>       // member
#include <util/containers/Pair.h>        // member template
#include <util/containers/DArray.h>      // member template

#include <pscf/chem/PolymerModel.h>
#include <pscf/chem/Edge.h>  
#include <util/containers/GArray.h>
#include <util/containers/FArray.h>
#include <util/containers/DMatrix.h>

#include <cmath>

namespace Pscf
{

   class Block;
   using namespace Util;

   /**
   * Descriptor and MDE solver for a block polymer.
   *
   * A PolymerTmpl<Block> object has arrays of Block and Vertex objects.
   * Each Block has two propagator MDE solver objects.  The solve() 
   * member function solves the modified diffusion equation (MDE) for 
   * all propagators in molecule.
   *
   * \ingroup Pscf_Solver_Module
   */
   template <class Block>
   class PolymerTmpl : public Species, public ParamComposite
   {

   public:

      // Modified diffusion equation solver for one block.
      typedef typename Block::Propagator Propagator;

      // Monomer concentration field.
      typedef typename Propagator::CField CField;

      // Chemical potential field.
      typedef typename Propagator::WField WField;

      /**
      * Constructor.
      */
      PolymerTmpl();

      /**
      * Destructor.
      */
      ~PolymerTmpl();

      /**
      * Read and initialize.
      *
      * \param in input parameter stream
      */
      virtual void readParameters(std::istream& in);

      /**
      * Solve modified diffusion equation.
      *
      * Upon return, propagators for all blocks, molecular partition
      * function q, phi and mu are all set.  The implementation of 
      * PolymerTmpl::solve() the calls the solve() function for all 
      * propagators in the molecule in a predeternined order, then
      * computes q, and finally computes mu from phi or phi from mu,
      * depending on the ensemble.
      *
      * This function should be called within a member function named
      * "compute" of each concrete subclass.  This compute function must
      * take an array of chemical potential fields (w-fields) as an
      * argument. Before calling the solve() function declared here,
      * the compute() function of each such subclass must setup the 
      * solvers by storing information required to solve the MDE within 
      * each block within a set of implementation dependent private 
      * members of the class that represents a block of a block polymer. 
      * After calling this solve() function, the compute function must
      * loop over blocks to compute the monomer concentration fields 
      * for all blocks.
      *
      * The optional parameter phiTot is only relevant to problems
      * such as thin films in which the material is excluded from part
      * of the unit cell by an inhogeneous constraint on the sum of
      * monomer concentration (i.e., a "mask").
      *
      * \param phiTot  fraction of volume occupied by material
      */
      virtual void solve(double phiTot = 1.0);

      /// \name Accessors (objects, by reference)
      ///@{

      /**
      * Get a specified Edge (block descriptor) by non-const reference.
      *
      * \param id block index, 0 <= id < nBlock
      */
      Edge& edge(int id);

      /**
      * Get a specified Edge (block descriptor) by const reference.
      *
      * \param id block index, 0 <= id < nBlock
      */
      Edge const& edge(int id) const;

      /**
      * Get a specified Block.
      *
      * \param id block index, 0 <= id < nBlock
      */
      Block& block(int id);

      /**
      * Get a specified Block by const reference.
      *
      * \param id block index, 0 <= id < nBlock
      */
      Block const& block(int id) const ;

      /**
      * Get a specified Vertex by const reference.
      *
      * Both chain ends and junctions are vertices.
      *
      * \param id vertex index, 0 <= id < nVertex
      */
      const Vertex& vertex(int id) const;

      /**
      * Get propagator for a specific block and direction (non-const).
      *
      * \param blockId integer index of associated block
      * \param directionId integer index for direction (0 or 1)
      */
      Propagator& propagator(int blockId, int directionId);

      /**
      * Get a propagator for a specific block and direction (const).
      *
      * \param blockId integer index of associated block
      * \param directionId integer index for direction (0 or 1)
      */
      Propagator const & propagator(int blockId, int directionId) const;

      /**
      * Get propagator indexed in order of computation (non-const).
      *
      * The propagator index must satisfy 0 <= id < 2*nBlock.
      *
      * \param id propagator index, in order of computation plan
      */
      Propagator& propagator(int id);

      /**
      * Get propagator identifier, indexed by order of computation.
      *
      * The return value is a pair of integers. The first integer
      * is a block index between 0 and nBlock - 1, and the second
      * is a propagator direction id, which must be 0 or 1.
      *
      * \param id  propagator index, in order of computation plan
      */
      Pair<int> const & propagatorId(int id) const;

      /**
      * Get propagator id from a source vertex towards a target.
      *
      * For is != it, the return value is an identifier for an outgoing
      * propagator that begins at the source vertex (vertex index is) 
      * and is part of the directed path that leads towards the target
      * vertex. In this case, the return value is a propagator identifier
      * similar to that returned by the propagatorId(int) member function,
      * which is a pair of integers in which the first element is block 
      * id and the second is a direction id for an outgoing propagator.
      *
      * For is == it, the return value is a pair [-1, -1]. 
      *
      * \param is  vertex index of the source vertex
      * \param it  vertex index of the target vertex
      */
      Pair<int> const & path(int is, int it) const;

      ///@}
      /// \name Accessors (by value)
      ///@{

      /**
      * Number of blocks.
      */
      int nBlock() const;

      /**
      * Number of vertices (junctions and chain ends).
      *
      * A theorem of graph theory tells us that, for any linear or
      * acyclic branched polymer, nVertex = nBlock + 1.
      */
      int nVertex() const;

      /**
      * Number of propagators (nBlock*2).
      */
      int nPropagator() const;  //

      /**
      * Sum of the lengths of all blocks in the polymer (thread model).
      *
      * Precondition: PolymerModel::isThread()
      */
      double length() const;

      /**
      * Total number of beads in the polymer (bead model).
      *
      * Precondition: PolymerModel::isBead()
      */
      int nBead() const;

      /**
      * Get Polymer type (Branched or Linear)
      */
      PolymerType::Enum type() const;

      ///@}

   protected:

      /**
      * Make a plan for order in which propagators should be computed.
      *
      * The algorithm creates a plan for computing propagators in an
      * that guarantees that the inital conditions required for each
      * propagator are known before it is processed. The algorithm is
      * works for any acyclic branched block copolymer. This function
      * is called in the default implementation of readParameters,
      * and must be called the readParameters method of any subclass.
      */
      virtual void makePlan();

      /**
      * Construct the paths data structure.
      */
      void makePaths();

      /**
      * Check the validity of the polymer graph.
      */
      void isValid();

   private:

      /// Array of Block objects in this polymer.
      DArray<Block> blocks_;

      /// Array of Vertex objects in this polymer.
      DArray<Vertex> vertices_;

      /// Propagator ids, indexed in order of computation.
      DArray< Pair<int> > propagatorIds_;

      /// Container for path-to-vertex signposts.
      DArray< DArray< Pair<int> > > paths_;

      /// Number of blocks in this polymer
      int nBlock_;

      /// Number of vertices (ends or junctions) in this polymer
      int nVertex_;

      /// Number of propagators (two per block).
      int nPropagator_;

      /// Polymer type (Branched or Linear)
      PolymerType::Enum type_;

   };

   // Inline functions

   /*
   * Number of vertices (ends and/or junctions)
   */
   template <class Block>
   inline int PolymerTmpl<Block>::nVertex() const
   {  return nVertex_; }

   /*
   * Number of blocks (or edges of the graph).
   */
   template <class Block>
   inline int PolymerTmpl<Block>::nBlock() const
   {  return nBlock_; }

   /*
   * Number of propagators.
   */
   template <class Block>
   inline int PolymerTmpl<Block>::nPropagator() const
   {  return nPropagator_; }

   /*
   * Total length of all blocks = volume / reference volume
   */
   template <class Block>
   inline double PolymerTmpl<Block>::length() const
   {
      UTIL_CHECK(PolymerModel::isThread());
      double value = 0.0;
      for (int blockId = 0; blockId < nBlock_; ++blockId) {
         value += edge(blockId).length();
      }
      return value;
   }

   /*
   * Total number of beads in all blocks.
   */
   template <class Block>
   inline int PolymerTmpl<Block>::nBead() const
   {
      UTIL_CHECK(PolymerModel::isBead());
      int value = 0.0;
      for (int blockId = 0; blockId < nBlock_; ++blockId) {
         value += edge(blockId).nBead();
      }
      return value;
   }

   /*
   * Get a specified Edge (block descriptor) by non-const reference.
   */
   template <class Block>
   inline Edge& PolymerTmpl<Block>::edge(int id)
   {  return blocks_[id]; }

   /*
   * Get a specified Edge (block descriptor) by const reference.
   */
   template <class Block>
   inline Edge const & PolymerTmpl<Block>::edge(int id) const
   {  return blocks_[id]; }

   /*
   * Get a specified Block solver by non-const reference.
   */
   template <class Block>
   inline Block& PolymerTmpl<Block>::block(int id)
   {  return blocks_[id]; }

   /*
   * Get a specified Block solver by const reference.
   */
   template <class Block>
   inline Block const & PolymerTmpl<Block>::block(int id) const
   {  return blocks_[id]; }

   /*
   * Get a specified Vertex by const reference.
   */
   template <class Block>
   inline
   Vertex const & PolymerTmpl<Block>::vertex(int id) const
   {  return vertices_[id]; }

   /*
   * Get a Propagator, indexed by block and direction ids (non-const).
   */
   template <class Block>
   inline
   typename Block::Propagator&
   PolymerTmpl<Block>::propagator(int blockId, int directionId)
   {  return block(blockId).propagator(directionId); }

   /*
   * Get a Propagator, indexed by block and direction ids (const).
   */
   template <class Block>
   inline
   typename Block::Propagator const &
   PolymerTmpl<Block>::propagator(int blockId, int directionId) const
   {  return block(blockId).propagator(directionId); }

   /*
   * Get a propagator indexed in order of computation (non-const).
   */
   template <class Block>
   inline
   typename Block::Propagator& PolymerTmpl<Block>::propagator(int id)
   {
      Pair<int> propId = propagatorId(id);
      return propagator(propId[0], propId[1]);
   }

   /*
   * Get a propagator id, indexed in order of computation.
   */
   template <class Block>
   inline
   Pair<int> const & PolymerTmpl<Block>::propagatorId(int id) const
   {
      UTIL_CHECK(id >= 0);
      UTIL_CHECK(id < nPropagator_);
      return propagatorIds_[id];
   }

   /*
   * Get a propagator id that leads from a source vertex towards a target.
   */
   template <class Block>
   inline
   Pair<int> const & PolymerTmpl<Block>::path(int is, int it) const
   {
      UTIL_CHECK(is >= 0);
      UTIL_CHECK(is < nVertex_);
      UTIL_CHECK(it >= 0);
      UTIL_CHECK(it < nVertex_);
      return paths_[is][it];
   }

   /*
   * Get the polymer type enumeration value (Branched or Linear).
   */
   template <class Block>
   inline PolymerType::Enum PolymerTmpl<Block>::type() const
   {  return type_; }

   // Non-inline functions

   /*
   * Constructor.
   */
   template <class Block>
   PolymerTmpl<Block>::PolymerTmpl()
    : Species(),
      blocks_(),
      vertices_(),
      propagatorIds_(),
      nBlock_(0),
      nVertex_(0),
      nPropagator_(0)
   {  setClassName("PolymerTmpl"); }

   /*
   * Destructor.
   */
   template <class Block>
   PolymerTmpl<Block>::~PolymerTmpl()
   {}

   template <class Block>
   void PolymerTmpl<Block>::readParameters(std::istream& in)
   {
      // Read polymer type (linear by default)
      type_ = PolymerType::Linear;
      readOptional<PolymerType::Enum>(in, "type", type_);

      read<int>(in, "nBlock", nBlock_);

      // Note: For any acyclic graph, nVertex = nBlock + 1
      nVertex_ = nBlock_ + 1;

      // Allocate all arrays (blocks_, vertices_ and propagatorIds_)
      blocks_.allocate(nBlock_);
      vertices_.allocate(nVertex_);
      propagatorIds_.allocate(2*nBlock_);

      // Set block id and polymerType for all blocks
      for (int blockId = 0; blockId < nBlock_; ++blockId) {
         edge(blockId).setId(blockId);
         edge(blockId).setPolymerType(type_);
      }

      // Set all vertex ids
      for (int vertexId = 0; vertexId < nVertex_; ++vertexId) {
         vertices_[vertexId].setId(vertexId);
      }

      // If polymer is linear polymer, set all block vertex Ids:
      if (type_ == PolymerType::Linear) {
         // In a linear chain, block i connects vertices i and i+1.
         for (int blockId = 0; blockId < nBlock_; ++blockId) {
            edge(blockId).setVertexIds(blockId, blockId + 1);
         }
         if (PolymerModel::isBead()) {
            // For bead model, set vertex ownership. For a linear chain:
            //    - block i owns vertex i+1.
            //    - block i owns vertex i iff i == 0.
            bool own0, own1;
            for (int blockId = 0; blockId < nBlock_; ++blockId) {
               own1 = true;
               own0 = false;
               if (blockId == 0) own0 = true;
               edge(blockId).setVertexOwnership(own0, own1);
            }
         }
      }

      // Read all other required block data
      readDArray<Block>(in, "blocks", blocks_, nBlock_);

      /*
      * The parameter file format for each block in the array blocks_
      * is different for branched and linear polymer, and different for
      * thread and bead models. These formats are defined by the >> and
      * << stream io operators for a Pscf::Block and are discussed in the
      * web manual
      *
      * For a linear polymer in the thread model, blocks must be entered 
      * sequentially, in the order they appear along the chain, and block 
      * vertex id values are set automatically. In this case, the format 
      * for one block is:
      *
      *   monomerId length
      *
      * By convention block i connects vertex i to vertex i+1. In the 
      * thread model, the block length is a real (i.e., floating point)
      * number. For a linear polymer in the bead model, the floating
      * point length parameter is replaced in this format by an integer
      * nBead that gives the number of beads in the block. 
      *
      * For a branched polymer in the thread model, the text format for 
      * each block is:
      *
      *   monomerId length vertexId(0) vertexId(1)
      *
      * where monomerId is the index of the block monomer type, length
      * is the block length, and vertexId(0) and vertexId(1) are the 
      * integer indices of the two vertices to which the block is 
      * attached.
      */

      // Read phi or mu (but not both)
      bool hasPhi = readOptional(in, "phi", phi_).isActive();
      if (hasPhi) {
         ensemble_ = Species::Closed;
      } else {
         ensemble_ = Species::Open;
         read(in, "mu", mu_);
      }

      // Reading of parameter file is now complete

      // Add edges to attached vertices
      int vertexId0, vertexId1;
      Edge* edgePtr;
      for (int blockId = 0; blockId < nBlock_; ++blockId) {
         edgePtr = &(edge(blockId));
         vertexId0 = edgePtr->vertexId(0);
         vertexId1 = edgePtr->vertexId(1);
         vertices_[vertexId0].addEdge(*edgePtr);
         vertices_[vertexId1].addEdge(*edgePtr);
      }

      // Polymer graph topology is now fully specified.

      // Construct a plan for the order in which block propagators
      // should be computed when solving the MDE.
      makePlan();

      // Construct path signposts
      makePaths();

      // Check internal consistency of edge and vertex data
      isValid();

      // Remainder of function - set and check propagator data

      // Set sources for all propagators
      Vertex const * vertexPtr = 0;
      Propagator const * sourcePtr = 0;
      Propagator * propagatorPtr = 0;
      Pair<int> propagatorId;
      int blockId, directionId, vertexId, i;
      for (blockId = 0; blockId < nBlock(); ++blockId) {
         // Add sources
         for (directionId = 0; directionId < 2; ++directionId) {
            vertexId = block(blockId).vertexId(directionId);
            vertexPtr = &vertex(vertexId);
            propagatorPtr = &block(blockId).propagator(directionId);
            for (i = 0; i < vertexPtr->size(); ++i) {
               propagatorId = vertexPtr->inPropagatorId(i);
               if (propagatorId[0] == blockId) {
                  UTIL_CHECK(propagatorId[1] != directionId);
               } else {
                  sourcePtr =
                     &block(propagatorId[0]).propagator(propagatorId[1]);
                  propagatorPtr->addSource(*sourcePtr);
               }
            } // end loop over inputs to head vertex
         } // end loop over directionId
      } // end loop over blockId

      // If using bead model, set propagator vertex ownership
      Block* blockPtr;
      if (PolymerModel::isBead()) {
         bool own0, own1;
         for (int blockId = 0; blockId < nBlock_; ++blockId) {
            blockPtr = &(blocks_[blockId]);
            own0 = blockPtr->ownsVertex(0);
            own1 = blockPtr->ownsVertex(1);
            blockPtr->propagator(0).setVertexOwnership(own0, own1);
            blockPtr->propagator(1).setVertexOwnership(own1, own0);
         }
      }

      // If using bead model, check ownership of vertex beads
      if (PolymerModel::isBead()) {
         int vSize, ib, id, nOwner;
         Pair<int> propId;

         // Loop over vertices
         for (int vertexId = 0; vertexId < nVertex_; ++vertexId) {
            vSize = vertices_[vertexId].size();

            // Incoming propagators
            nOwner = 0;
            for (int ip = 0; ip < vSize; ++ip) {
               propId = vertices_[vertexId].inPropagatorId(ip);
               ib = propId[0];
               id = propId[1];
               if (blocks_[ib].propagator(id).ownsTail()) {
                 ++nOwner;
               }
            }
            UTIL_CHECK(nOwner == 1);

            // Outgoing propagators
            nOwner = 0;
            for (int i = 0; i < vSize; ++i) {
               propId = vertices_[vertexId].outPropagatorId(i);
               ib = propId[0];
               id = propId[1];
               if (blocks_[ib].propagator(id).ownsHead()) {
                 ++nOwner;
               }
            }
            UTIL_CHECK(nOwner == 1);

         } // end loop over vertices

         // Check consistency of block and propagator ownership flags
         bool own0, own1;
         for (int blockId = 0; blockId < nBlock_; ++blockId) {
            blockPtr = &(blocks_[blockId]);
            own0 = blockPtr->ownsVertex(0);
            own1 = blockPtr->ownsVertex(1);
            UTIL_CHECK(blockPtr->propagator(0).ownsHead() == own0);
            UTIL_CHECK(blockPtr->propagator(1).ownsTail() == own0);
            UTIL_CHECK(blockPtr->propagator(0).ownsTail() == own1);
            UTIL_CHECK(blockPtr->propagator(1).ownsHead() == own1);
         } 

      } // end if (PolymerModel::isBead())

   }

   /*
   * Make a plan for the order of computation of block propagators.
   */
   template <class Block>
   void PolymerTmpl<Block>::makePlan()
   {
      if (nPropagator_ != 0) {
         UTIL_THROW("nPropagator !=0 on entry");
      }

      // Allocate and initialize isFinished matrix
      DMatrix<bool> isFinished;
      isFinished.allocate(nBlock_, 2);
      for (int iBlock = 0; iBlock < nBlock_; ++iBlock) {
         for (int iDirection = 0; iDirection < 2; ++iDirection) {
            isFinished(iBlock, iDirection) = false;
         }
      }

      Pair<int> propagatorId;
      Vertex* inVertexPtr = 0;
      int inVertexId = -1;
      bool isReady;
      while (nPropagator_ < nBlock_*2) {
         for (int iBlock = 0; iBlock < nBlock_; ++iBlock) {
            for (int iDirection = 0; iDirection < 2; ++iDirection) {
               if (isFinished(iBlock, iDirection) == false) {
                  inVertexId = edge(iBlock).vertexId(iDirection);
                  inVertexPtr = &vertices_[inVertexId];
                  isReady = true;
                  for (int j = 0; j < inVertexPtr->size(); ++j) {
                     propagatorId = inVertexPtr->inPropagatorId(j);
                     if (propagatorId[0] != iBlock) {
                        if (!isFinished(propagatorId[0], propagatorId[1])){
                           isReady = false;
                           break;
                        }
                     }
                  }
                  if (isReady) {
                     propagatorIds_[nPropagator_][0] = iBlock;
                     propagatorIds_[nPropagator_][1] = iDirection;
                     isFinished(iBlock, iDirection) = true;
                     ++nPropagator_;
                  }
               }
            }
         }
      }

   }

   /*
   * Compute a solution to the MDE and block concentrations.
   */
   template <class Block>
   void PolymerTmpl<Block>::solve(double phiTot)
   {

      // Clear all propagators
      for (int j = 0; j < nPropagator(); ++j) {
         propagator(j).setIsSolved(false);
      }

      // Solve modified diffusion equation for all propagators in
      // the order specified by function makePlan.
      for (int j = 0; j < nPropagator(); ++j) {
         UTIL_CHECK(propagator(j).isReady());
         propagator(j).solve();
      }

      // Compute molecular partition function q_
      q_ = block(0).propagator(0).computeQ();

      // The Propagator::computeQ function returns a spatial average.
      // Correct for partial occupation of the unit cell.
      q_ = q_/phiTot;

      // Compute mu_ or phi_, depending on ensemble
      if (ensemble() == Species::Closed) {
         mu_ = log(phi_/q_);
      } else
      if (ensemble() == Species::Open) {
         phi_ = exp(mu_)*q_;
      }

      #if 0
      // Compute block concentration fields
      double prefactor;
      if (PolymerModel::isThread()) {
         prefactor = phi_ / ( q_ * length() );
      } else 
      if (PolymerModel::isBead()) {
         prefactor = phi_ / ( q_ * (double)nBead() );
      }
      for (int i = 0; i < nBlock(); ++i) {
         block(i).computeConcentration(prefactor);
      }
      #endif

   }

   /*
   * Identify paths between vertices.
   */
   template <class Block>
   void PolymerTmpl<Block>::makePaths()
   {
      UTIL_CHECK(nVertex_ > 0);

      // Local variables
      DArray< GArray< FArray<int, 3> > >  sendPaths;
      DArray< GArray< FArray<int, 3> > >  recvPaths;
      DArray< GArray< FArray<int, 3> > >  oldPaths;
      FArray<int, 3> path;
      Pair<int> pId;
      int is; // sender vertex id
      int ir; // receiver vertex id
      int ib; // bond id
      int io; // outgoing (source) direction id
      int ii; // incoming (receiver) direction id
      int j, k, n;

      // Allocate and clear path containers
      sendPaths.allocate(nVertex_);
      recvPaths.allocate(nVertex_);
      oldPaths.allocate(nVertex_);
      for (is = 0; is < nVertex_; ++is) {
         sendPaths[is].clear();
         recvPaths[is].clear();
         oldPaths[is].clear();
      }

      // Initialize sendPaths with path-to-self entries
      for (is = 0; is < nVertex_; ++is) {
         path[0] = is;
         path[1] = -1;
         path[2] = -1;
         sendPaths[is].append(path);
      }

      // While loop to completion
      bool done = false;
      int iter = 0;
      while (!done) {
         UTIL_CHECK(iter < nBlock_ + 2);

         // Check that recvPaths container is empty
         for (ir = 0; is < nVertex_; ++is) {
            UTIL_CHECK(recvPaths[ir].size() == 0);
         }

         // Loop over source vertices for sending
         done = true;
         for (is = 0; is < nVertex_; ++is) {

            // Send any sendPaths data to neighbors
            int n = sendPaths[is].size();
            if (n > 0) {
               done = false;
               Vertex const & sender = vertex(is);

               // Send sendPaths data to all new neighbors
               for (j = 0; j < sender.size(); ++j) {
                  pId = sender.outPropagatorId(j);
                  ib = pId[0];   // block identifier
                  io = pId[1];   // outgoing direction identifier
                  if (io == 0) {
                     UTIL_CHECK(edge(ib).vertexId(0) == is);
                     ii = 1;
                  } else {
                     UTIL_CHECK(edge(ib).vertexId(1) == is);
                     ii = 0;
                  }
                  ir = edge(ib).vertexId(ii);
                  for (k = 0; k < n; ++k) {
                     path = sendPaths[is][k];
                     // Send unless just received along same bond
                     if (ib != path[1]) {
                        path[1] = ib;
                        path[2] = ii;
                        recvPaths[ir].append(path);
                     }
                  }
               }

            } // if (n > 0)
         } // Loop over source vertices

         // Loop over vertices to transfer data structures
         for (is = 0; is < nVertex_; ++is) {

            // Move sendPaths to oldPaths, clear sendPaths
            n = sendPaths[is].size();
            if (n > 0) {
               for (k = 0; k < n; ++k) {
                  path = sendPaths[is][k];
                  oldPaths[is].append(path);
               }
            }
            sendPaths[is].clear();

            // Move recvPaths to sendPaths, clear recvPaths
            n = recvPaths[is].size();
            if (n > 0) {
               for (k = 0; k < n; ++k) {
                  sendPaths[is].append(recvPaths[is][k]);
               }
            }
            recvPaths[is].clear();

         }

         ++iter;
      } // while not done

      // Allocate and initialize member variable paths_
      paths_.allocate(nVertex_);
      for (is = 0; is < nVertex_; ++is) {
         paths_[is].allocate(nVertex_);
         for (ir = 0; ir < nVertex_; ++ir) {
            paths_[is][ir][0] = -1;
            paths_[is][ir][1] = -1;
         }
      }

      // Assign values to all elements of paths_ container
      for (is = 0; is < nVertex_; ++is) {
         n = oldPaths[is].size();
         UTIL_CHECK(n == nVertex_);
         for (k = 0; k < n; ++k) {
            path = oldPaths[is][k];
            ir = path[0]; // id of target vertex
            UTIL_CHECK(ir >= 0);
            UTIL_CHECK(ir < nVertex_);
            // Check that element was not set previously
            UTIL_CHECK(paths_[is][ir][0] == -1);
            UTIL_CHECK(paths_[is][ir][1] == -1);
	    //std::cout << std::endl << is << " " << ir << " " 
            //          << path[1] << " " << path[2];
            if (ir == is) {
               UTIL_CHECK(path[1] == -1);
               UTIL_CHECK(path[2] == -1);
            } else {
               UTIL_CHECK(path[1] >= 0);
               UTIL_CHECK(path[1] < nBlock_);
               UTIL_CHECK(path[2] >= 0);
               UTIL_CHECK(path[2] < 2);
            }
            paths_[is][ir][0] = path[1];
            paths_[is][ir][1] = path[2];
         }
      }

      #if 0
      std::cout << std::endl << "Paths:" ;
      for (is = 0; is < nVertex_; ++is) {
         for (ir = 0; ir < nVertex_; ++ir) {
            pId = paths_[is][ir];
            std::cout << std::endl;
            std::cout << is << ir << pId[0] << pId[1];
         }
      }
      std::cout << std::endl;
      #endif

   }

   /*
   * Check consistency of data structures.
   */
   template <class Block>
   void PolymerTmpl<Block>::isValid()
   {
      Pair<int> pair;
      int ib, iv, ip, id, iv0, iv1, n;

      // Check validity of ids owned by blocks
      for (ib = 0; ib < nBlock_; ++ib) {
         UTIL_CHECK(edge(ib).id() == ib);
         iv0 = edge(ib).vertexId(0);
         iv1 = edge(ib).vertexId(1);
         UTIL_CHECK(iv0 != iv1);
         UTIL_CHECK(iv0 >= 0);
         UTIL_CHECK(iv0 < nVertex_);
         UTIL_CHECK(iv1 >= 0);
         UTIL_CHECK(iv1 < nVertex_);
      }

      // Check consistency of vertex::outPropagatorId
      for (iv = 0; iv < nVertex_; ++iv) {
         UTIL_CHECK(vertex(iv).id() == iv);
         n = vertex(iv).size();
         for (ip = 0; ip < n; ++ip) {
            pair = vertex(iv).outPropagatorId(ip);
            ib = pair[0];
            id = pair[1];
            UTIL_CHECK(ib >= 0);
            UTIL_CHECK(ib < nBlock_);
            UTIL_CHECK(id >= 0);
            UTIL_CHECK(id < 2);
            UTIL_CHECK(edge(ib).vertexId(id) == iv);
         }
      }

      // Check consistency of vertex::inPropagatorId
      for (iv = 0; iv < nVertex_; ++iv) {
         UTIL_CHECK(vertex(iv).id() == iv);
         n = vertex(iv).size();
         for (ip = 0; ip < n; ++ip) {
            pair = vertex(iv).inPropagatorId(ip);
            ib = pair[0];
            id = pair[1];
            UTIL_CHECK(ib >= 0);
            UTIL_CHECK(ib < nBlock_);
            UTIL_CHECK(id >= 0);
            UTIL_CHECK(id < 2);
            if (id == 0) {
               UTIL_CHECK(edge(ib).vertexId(1) == iv);
            } else {
               UTIL_CHECK(edge(ib).vertexId(0) == iv);
            }
         }
      }

      // Check consistency of vertex ids in paths
      for (iv0 = 0; iv0 < nVertex_; ++iv0) {
         for (iv1 = 0; iv1 < nVertex_; ++iv1) {
            pair = path(iv0, iv1);
            ib = pair[0];
            id = pair[1];
            if (iv0 == iv1) {
               UTIL_CHECK(ib == -1);
               UTIL_CHECK(id == -1);
            } else {
               UTIL_CHECK(edge(ib).vertexId(id) == iv0);
            }
         }
      }
     
   }

}
#endif

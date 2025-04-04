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

#include <util/containers/GArray.h>
#include <util/containers/FArray.h>
#include <util/containers/DMatrix.h>

#include <cmath>

namespace Pscf
{

   class Block;
   using namespace Util;

   /**
   * Descriptor and MDE solver for an acyclic block polymer.
   *
   * A PolymerTmpl<Block> object has arrays of Block and Vertex
   * objects. Each Block has two propagator MDE solver objects.
   * The solve() member function solves the modified diffusion
   * equation (MDE) for the entire molecule and computes monomer
   * concentration fields for all blocks.
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
      * function q, phi and mu, and the block concentration fields for
      * all blocks are set.
      *
      * This function should be called within a member function named
      * "compute" of each concrete subclass.  This compute function must
      * take an array of chemical potential fields (w-fields) as an
      * argument. Before calling the solve() function declared here,
      * the compute() function of each such subclass must store
      * information required to solve the MDE within each block within
      * a set of implementation dependent private members of the class
      * class that represents a block of a block polymer. The
      * implementation of PolymerTmpl::solve() the calls the solve()
      * function for all propagators in the molecule in a predeternined
      * order.
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
      * Get propagator indexed in order of computation (non-const).
      *
      * The propagator index must satisfy 0 <= id < 2*nBlock.
      *
      * \param id propagator index, in order of computation plan
      */
      Propagator&  propagator(int id);

      /**
      * Get a propagator for a specific block and direction (const).
      *
      * \param blockId integer index of associated block
      * \param directionId integer index for direction (0 or 1)
      */
      Propagator const & propagator(int blockId, int directionId) const;

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
      * Get propagator identifier from one vertex towards a target.
      *
      * For is != it, the return value is an identifier for an outgoing
      * propagator that begins at the source vertex (vertex index is) 
      * and is part of the directed path that leads towards the 
      * target vertex. In this case, the return value is a propagator
      * identifier similar to that returned by the propagatorId(int)
      * member function, for which the first element is an int block id
      * and the second is a direction id for an outgoing propagator.
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
      * Sum of the lengths of all blocks in the polymer.
      */
      double length() const;

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

      void makePaths();

      void isValid();

   private:

      /// Array of Block objects in this polymer.
      DArray<Block> blocks_;

      /// Array of Vertex objects in this polymer.
      DArray<Vertex> vertices_;

      /// Propagator ids, indexed in order of computation.
      DArray< Pair<int> > propagatorIds_;

      /// Container for path-to-vertex signposts
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

   /*
   * Number of vertices (ends and/or junctions)
   */
   template <class Block>
   inline int PolymerTmpl<Block>::nVertex() const
   {  return nVertex_; }

   /*
   * Number of blocks.
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
      double value = 0.0;
      for (int blockId = 0; blockId < nBlock_; ++blockId) {
         value += blocks_[blockId].length();
      }
      return value;
   }

   /*
   * Get a specified Vertex.
   */
   template <class Block>
   inline
   const Vertex& PolymerTmpl<Block>::vertex(int id) const
   {  return vertices_[id]; }

   /*
   * Get a specified Block.
   */
   template <class Block>
   inline Block& PolymerTmpl<Block>::block(int id)
   {  return blocks_[id]; }

   /*
   * Get a specified Block by const reference.
   */
   template <class Block>
   inline Block const & PolymerTmpl<Block>::block(int id) const
   {  return blocks_[id]; }

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
   * Get a propagator indexed by block and direction.
   */
   template <class Block>
   inline
   typename Block::Propagator&
   PolymerTmpl<Block>::propagator(int blockId, int directionId)
   {  return block(blockId).propagator(directionId); }

   /*
   * Get a const propagator indexed by block and direction.
   */
   template <class Block>
   inline
   typename Block::Propagator const &
   PolymerTmpl<Block>::propagator(int blockId, int directionId) const
   {  return block(blockId).propagator(directionId); }

   /*
   * Get a propagator indexed in order of computation.
   */
   template <class Block>
   inline
   typename Block::Propagator& PolymerTmpl<Block>::propagator(int id)
   {
      Pair<int> propId = propagatorId(id);
      return propagator(propId[0], propId[1]);
   }

   /*
   * Get a propagator id that leads from one vertex towards another.
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
   * Get Polymer type (Branched or Linear).
   */
   template <class Block>
   inline PolymerType::Enum PolymerTmpl<Block>::type() const
   { return type_; }

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
         blocks_[blockId].setId(blockId);
         blocks_[blockId].setPolymerType(type_);
      }

      // Set all vertex ids
      for (int vertexId = 0; vertexId < nVertex_; ++vertexId) {
         vertices_[vertexId].setId(vertexId);
      }

      // If polymer is linear polymer, set all block vertex Ids:
      // In a linear chain, block i connects vertex i and vertex i+1.
      if (type_ == PolymerType::Linear) {
         for (int blockId = 0; blockId < nBlock_; ++blockId) {
            blocks_[blockId].setVertexIds(blockId, blockId + 1);
         }
      }

      // Read all other required block data
      readDArray<Block>(in, "blocks", blocks_, nBlock_);

      /*
      * The parameter file format for each block in the array blocks_
      * is different for branched and linear polymer. These formats
      * are defined in >> and << stream io operators for a Pscf::Block.
      * The choice of format is controlled by the Block::polymerType.
      *
      * For a branched polymer, the text format for each block is:
      *
      *   monomerId length vertexId(0) vertexId(1)
      *
      * where monomerId is the index of the block monomer type, length
      * is the block length, and vertexId(0) and vertexId(1) are the
      * indices of the two vertices to which the block is attached.
      *
      * For a linear polymer, blocks must be entered sequentially, in
      * their order along the chain, and block vertex id values are
      * set automatically. In this case, the format for one block is:
      *
      *   monomerId length
      *
      * with no need for explicit vertex ids.
      */

      // Add blocks to attached vertices
      int vertexId0, vertexId1;
      Block* blockPtr;
      for (int blockId = 0; blockId < nBlock_; ++blockId) {
         blockPtr = &(blocks_[blockId]);
         vertexId0 = blockPtr->vertexId(0);
         vertexId1 = blockPtr->vertexId(1);
         vertices_[vertexId0].addEdge(*blockPtr);
         vertices_[vertexId1].addEdge(*blockPtr);
      }

      // Polymer topology is now fully specified.

      // Construct a plan for the order in which block propagators
      // should be computed when solving the MDE.
      makePlan();

      // Construct paths
      makePaths();

      // Read phi or mu (but not both)
      bool hasPhi = readOptional(in, "phi", phi_).isActive();
      if (hasPhi) {
         ensemble_ = Species::Closed;
      } else {
         ensemble_ = Species::Open;
         read(in, "mu", mu_);
      }

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
            }
         }

      }

      // Check internal consistency of data structures
      isValid();

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
                  inVertexId = blocks_[iBlock].vertexId(iDirection);
                  inVertexPtr = &vertices_[inVertexId];
                  isReady = true;
                  for (int j = 0; j < inVertexPtr->size(); ++j) {
                     propagatorId = inVertexPtr->inPropagatorId(j);
                     if (propagatorId[0] != iBlock) {
                        if (!isFinished(propagatorId[0], propagatorId[1])) {
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

      // Compute block concentration fields
      double prefactor = phi_ / ( q_ * length() );
      for (int i = 0; i < nBlock(); ++i) {
         block(i).computeConcentration(prefactor);
      }

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
                     UTIL_CHECK(block(ib).vertexId(0) == is);
                     ii = 1;
                  } else {
                     UTIL_CHECK(block(ib).vertexId(1) == is);
                     ii = 0;
                  }
                  ir = block(ib).vertexId(ii);
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
         UTIL_CHECK(block(ib).id() == ib);
         iv0 = block(ib).vertexId(0);
         iv1 = block(ib).vertexId(1);
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
            UTIL_CHECK(block(ib).vertexId(id) == iv);
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
               UTIL_CHECK(block(ib).vertexId(1) == iv);
            } else {
               UTIL_CHECK(block(ib).vertexId(0) == iv);
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
               UTIL_CHECK(block(ib).vertexId(id) == iv0);
            }
         }
      }
     
   }

}
#endif

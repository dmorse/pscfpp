/*
* PSCF - Polymer Self-Consistent Field Theory
*
* Copyright 2016 - 2022, The Regents of the University of Minnesota
* Distributed under the terms of the GNU General Public License.
*/

#include "PolymerSpecies.h"

#include <pscf/chem/Edge.h>  

#include <util/containers/GArray.h>
#include <util/containers/FArray.h>
#include <util/containers/DMatrix.h>

#include <cmath>

namespace Pscf
{

   /*
   * Constructor.
   */
   PolymerSpecies::PolymerSpecies()
    : Species(),
      vertices_(),
      propagatorIds_(),
      paths_(),
      nBlock_(0),
      nVertex_(0),
      nPropagator_(0),
      type_(PolymerType::Linear)
   {  setClassName("PolymerSpecies"); }

   /*
   * Destructor.
   */
   PolymerSpecies::~PolymerSpecies()
   {}

   /*
   * Read parameter file block.
   */
   void PolymerSpecies::readParameters(std::istream& in)
   {
      // Read polymer type (linear by default)
      type_ = PolymerType::Linear;
      readOptional<PolymerType::Enum>(in, "type", type_);

      read<int>(in, "nBlock", nBlock_);

      // Note: For any acyclic graph, nVertex = nBlock + 1
      nVertex_ = nBlock_ + 1;

      // Allocate all arrays
      allocateBlocks();
      vertices_.allocate(nVertex_);
      propagatorIds_.allocate(2*nBlock_);

      // Set block id and polymerType for all blocks
      for (int edgeId = 0; edgeId < nBlock_; ++edgeId) {
         edge(edgeId).setId(edgeId);
         edge(edgeId).setPolymerType(type_);
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

      // Read array of block data from parameter file
      readBlocks(in);

      // Read phi or mu (but not both) and set ensemble accordingly
      Species::readParameters(in);

      // Reading of parameter file is now complete

      // Add edges to attached vertices
      int vertexId0, vertexId1;
      Edge* edgePtr;
      for (int edgeId = 0; edgeId < nBlock_; ++edgeId) {
         edgePtr = &(edge(edgeId));
         vertexId0 = edgePtr->vertexId(0);
         vertexId1 = edgePtr->vertexId(1);
         vertices_[vertexId0].addEdge(*edgePtr);
         vertices_[vertexId1].addEdge(*edgePtr);
      }

      // Polymer graph topology is now fully specified.

      // Construct a plan for the order in which block propagators
      // should be computed when solving the MDE.
      makePlan();

      // Construct vertex-to-vertex path signposts
      makePaths();

      // Check internal consistency of edge and vertex data
      isValid();

   }

   /*
   * Make a plan for the order of computation of block propagators.
   */
   void PolymerSpecies::makePlan()
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

      Pair<int> propId;
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
                     propId = inVertexPtr->inPropagatorId(j);
                     if (propId[0] != iBlock) {
                        if (!isFinished(propId[0], propId[1])){
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
   * Identify paths between vertices.
   */
   void PolymerSpecies::makePaths()
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
   void PolymerSpecies::isValid()
   {
      Pair<int> pair;
      int ib, iv, ip, id, iv0, iv1, n;

      // Check validity of ids owned by edges
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

   /*
   * Total length of all blocks = volume / reference volume
   */
   double PolymerSpecies::length() const
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
   int PolymerSpecies::nBead() const
   {
      UTIL_CHECK(PolymerModel::isBead());
      int value = 0.0;
      for (int blockId = 0; blockId < nBlock_; ++blockId) {
         value += edge(blockId).nBead();
      }
      return value;
   }

}

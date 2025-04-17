#ifndef PSCF_POLYMER_SPECIES_H
#define PSCF_POLYMER_SPECIES_H

/*
* PSCF - Polymer Self-Consistent Field Theory
*
* Copyright 2016 - 2022, The Regents of the University of Minnesota
* Distributed under the terms of the GNU General Public License.
*/

#include <pscf/chem/Species.h>           // base class

#include <pscf/chem/Vertex.h>            // member
#include <pscf/chem/PolymerType.h>       // member
#include <pscf/chem/PolymerModel.h>      // member
#include <util/containers/Pair.h>        // member
#include <util/containers/DArray.h>      // member

namespace Pscf
{

   // Forward declaration
   class Edge;

   using namespace Util;

   /**
   * Descriptor for a linear or acyclic branched block polymer.
   *
   * A PolymerSpecies has an array of Vertex objects and an associated
   * array of Edge objects. Together, these provide a description of the
   * connectivity graph, block lengths and block monomer types of blocks 
   * within a block polymer.  A PolymerSpecies object does not provide
   * function required to solve the modified diffusion equation (MDE), 
   * and so is a descriptor but not a solver for a polymer species. 
   *
   * PolymerSpecies is an abstract base class for classes that act as
   * MDE solvers as well as descriptors.  Each implementation level 
   * sub-namespace of Pscf (i.e., R1d, Rpc or Rpg) contains a subclass 
   * of PolymerSpecies named Polymer that defines a function named 
   * compute() that solves the MDE for all blocks of a block polymer.  
   * Each such implementation level sub-namespace also contains a 
   * subclass of Edge named Block that is an MDE solver and descriptor 
   * for a single block within a block polymer.
   *
   * The Polymer class in each implementation-level namespace is derived 
   * directly from an template instantiation Pscf::PolymerTmpl<Block> of 
   * the class template Pscf::PolymerTmpl, while PolymerTmpl<Block> is 
   * derived directly from PolymerSpecies. A PolymerTmpl<Block> has a 
   * private member variable that is an array of Block solver objects. 
   * The PolymerTmpl template defines implementations of edge(int id) 
   * functions that are declared as pure virtual functions of 
   * PolymerSpecies, each of which which returns a specific Block solver 
   * object as a reference to an Edge.
   *
   * The PolymerSpecies class provides a generic description of polymer
   * structure that only uses names of generic descriptor types (such as
   * as Edge and Vertex) that are defined directly in the Pscf namespace.
   * The PolymerSpecies interface does not use any of the specialized 
   * solver classes that are defined in each implementation level 
   * sub-namespaces of Pscf, such as the different Block and Polymer 
   * classes defined in different sub-namespaces for use in different 
   * implementations. A reference to such a specialized Polymer solver
   * object as a generic PolymerSpecies can thus be used as an argument
   * to a generic function that requires a description of the polymer
   * molecular structure, but that does not actually solve the MDE.
   * 
   * \ingroup Pscf_Chem_Module
   */
   class PolymerSpecies : public Species
   {

   public:

      /**
      * Constructor.
      */
      PolymerSpecies();

      /**
      * Destructor.
      */
      ~PolymerSpecies();

      /**
      * Read parameters and initialize.
      *
      * The default implementation of this function calls the virtual
      * protected functions allocateBlocks() and readBlocks() to access
      * data that is owned by the PolymerTmpl<Block> subclass.
      *
      * \param in input parameter stream
      */
      virtual void readParameters(std::istream& in);

      /// \name Accessors (return objects by reference)
      ///@{

      /**
      * Get a specified Edge (block descriptor) by non-const reference.
      *
      * This function is implemented by the class template 
      * Pscf::PolymerTmpl<class Block>, in which Block must be a subclass
      * of Pscf::Edge.
      *
      * \param id block index, 0 <= id < nBlock
      */
      virtual Edge& edge(int id) = 0;

      /**
      * Get a specified Edge (block descriptor) by const reference.
      *
      * \param id block index, 0 <= id < nBlock
      */
      virtual Edge const& edge(int id) const = 0;

      /**
      * Get a specified Vertex by const reference.
      *
      * Chain ends and block junctions are all vertices.
      *
      * \param id vertex index, 0 <= id < nVertex
      */
      const Vertex& vertex(int id) const;

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
      * Number of propagators (2*nBlock).
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
      /// \name Directed Edges and Connectivity
      ///@{

      /**
      * Get a propagator identifier, indexed by order of computation.
      *
      * The return value is a pair of integers that identifies a
      * directed edge, or a propagator. The first integer is a block
      * index between 0 and nBlock - 1, and the second is a propagator
      * direction id, which must be 0 or 1. By convention, direction 0
      * of Edge i propagates from vertex edge(i).vertexId(0) to vertex
      * edge(i).vertexId(1), while direction 1 propagates in the other
      * direction, from edge(i).vertexId(1) to edge(i).vertexId(0).
      *
      * The propagators are indexed in an order that is computed by the 
      * protected function PolymerSpecies::makePlan. The indexing is 
      * chosen such that, if MDEs for different propagators are solved in 
      * sequential order by propagator index, then information required 
      * to construct the initial condition for the MDE within each 
      * propagator (the source slice) is guaranteed to be available from 
      * the solution of previous calculations.
      *
      * \param id  propagator index, in order of computation plan
      */
      Pair<int> const & propagatorId(int id) const;

      /**
      * Get an id for a propagator from one vertex towards a target.
      *
      * For is != it, the return value is an identifier for an outgoing
      * propagator that begins at the source vertex (vertex index is) and
      * is the first edge of the directed path that leads from the source
      * vertex to the target vertex (vertex id it). The return value is a
      * pair of integers analogous to that returned by the propagatorId(int)
      * member function, for which the first element is a block index and 
      * the second element is a direction id (0 or 1) for the propagator 
      * direction for that block that is outgoing from from the source 
      * vertex along the unique path leading to the target vertex. The
      * requirement that the propagator be outgoing from vertex is implies 
      * that if path(is, it) == [ib, id], then edge(ib).vertexId(id) == is.
      *
      * \param is  index of the source vertex
      * \param it  index of the target vertex
      */
      Pair<int> const & path(int is, int it) const;

      ///@}

   protected:

      /**
      * Allocate array of blocks.
      */
      virtual void allocateBlocks() = 0;

      /**
      * Read array of blocks from parameter file.
      *
      * \param in  parameter input stream
      */
      virtual void readBlocks(std::istream& in) = 0;

      /**
      * Make a plan for order in which propagators should be computed.
      *
      * The algorithm creates a plan for processing propagators in an order
      * that guarantees that the inital conditions required for solution of
      * the MDE for each propagator is known before it is processed. The
      * resulting sequence can be accessed by the propagatorId(int id)
      * member function. This function should be called in the 
      * readParameters member function.
      *
      * The algorithm implemented by this function will yield a valid
      * sequence for any any acyclic branched block polymer. There can,
      * however, exist two or more valid sequences of computations that 
      * satisfy the requirements imposed here. Different valid sequences
      * should yield to nearly identical computational performance in a 
      * strictly sequential algorithm. The algorithm used in this function 
      * was not designed to choose from among different possible valid 
      * sequences based on consideration of possible extension to a 
      * threaded implementation, or on any other performance consideration.
      */
      virtual void makePlan();

      /**
      * Create a matrix of vertex-to-vertex path signposts.
      *
      * This function constructs the data structure that is accessed by
      * the PolymerSpecies::path member function. This function should be
      * called in the readParameters member function.
      */
      void makePaths();

      /**
      * Check validity of the polymer graph.
      *
      * An Exception is thrown if the graph is not valid.
      */
      void isValid();

   private:

      /// Array of Vertex objects in this polymer.
      DArray<Vertex> vertices_;

      /// Propagator ids, indexed in order of computation.
      DArray< Pair<int> > propagatorIds_;

      /// Container for path-to-vertex signposts (see paths function).
      DArray< DArray< Pair<int> > > paths_;

      /// Number of blocks in this polymer
      int nBlock_;

      /// Number of vertices (ends and block junctions) in this polymer.
      int nVertex_;

      /// Number of propagators (two per block).
      int nPropagator_;

      /// Polymer type (Branched or Linear)
      PolymerType::Enum type_;

   };

   // Inline functions

   /*
   * Number of blocks in this polymer.
   */
   inline int PolymerSpecies::nBlock() const
   {  return nBlock_; }

   /*
   * Number of vertices (ends and/or junctions)
   */
   inline int PolymerSpecies::nVertex() const
   {  return nVertex_; }

   /*
   * Number of propagators in this polymer (2*nBlock).
   */
   inline int PolymerSpecies::nPropagator() const
   {  return nPropagator_; }

   /*
   * Get a specified Vertex by const reference.
   */
   inline
   Vertex const & PolymerSpecies::vertex(int id) const
   {  return vertices_[id]; }

   /*
   * Get a propagator id, indexed in order of computation.
   */
   inline
   Pair<int> const & PolymerSpecies::propagatorId(int id) const
   {
      UTIL_CHECK(id >= 0);
      UTIL_CHECK(id < nPropagator_);
      return propagatorIds_[id];
   }

   /*
   * Get a propagator id that leads from a source vertex towards a target.
   */
   inline
   Pair<int> const & PolymerSpecies::path(int is, int it) const
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
   inline
   PolymerType::Enum PolymerSpecies::type() const
   {  return type_; }

}
#endif

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
   * array of Edge objects that, together, provide a description of the
   * connectivity and properties of blocks with a polymer polymer, and of 
   * the associated graph.  A PolymerSpecies object is a descriptor for
   * for a polymer species, but does not contain functions or data 
   * structures required to solve the modified diffusion equation (MDE), 
   * and so is not a solver class.
   *
   * Each implementation level sub-namespace of Pscf (R1d, Rpc or
   * Rpg) contains a concrete class named Polymer that is a subclass of
   * Pscf::Polymer species, and that acts as both a descriptor and MDE 
   * solver for the associated species.  Each such implementation level 
   * sub-namespace also contains a class named Block that is a subclass
   * of Pscf::Edge, and that is a descriptor and MDE solver for a single 
   * block. The Polymer class in each such namespace is a subclass of 
   * a class template instantiation PolymerTmpl<Block> that is itself a
   * subclass of PolymerSpecies. 
   * 
   * A PolymerTmpl<Block> class has a private member variable that
   * is an array of Block objects.  The PolymerTmpl class template 
   * defines implementations of the pure virtual edge(int id) member 
   * functions declared by PolymerSpecies, which return a single element
   * of this array as a reference to an Edge. 
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
      * \param in input parameter stream
      */
      virtual void readParameters(std::istream& in);

      /// \name Accessors (objects, by reference)
      ///@{

      /**
      * Get a specified Edge (block descriptor) by non-const reference.
      *
      * This function is implemented by the class template 
      * Pscf::PolymerTmpl<class Block>, which is derived from Pscf::Edge.
      *
      * \param id block index, 0 <= id < nBlock
      */
      virtual Edge& edge(int id) = 0;

      /**
      * Get a specified Edge (block descriptor) by const reference.
      *
      * This function is implemented by the class template 
      * Pscf::PolymerTmpl<class Block>, which is derived from Pscf::Edge.
      *
      * \param id block index, 0 <= id < nBlock
      */
      virtual Edge const& edge(int id) const = 0;

      /**
      * Get a specified Vertex by const reference.
      *
      * Both chain ends and junctions are vertices.
      *
      * \param id vertex index, 0 <= id < nVertex
      */
      const Vertex& vertex(int id) const;

      /**
      * Get a propagator identifier, indexed by order of computation.
      *
      * The return value is a pair of integers that identifies a
      * directed edge, or a propagator. The first integer is a block
      * index between 0 and nBlock - 1, and the second is a propagator
      * direction id, which must be 0 or 1. By convention, direction 0
      * of edge propagates from vertex edge(i).vertexId(0) to vertex
      * edge(i).vertexId(1), while direction 1 propagates in the
      * reverse direction.
      *
      * \param id  propagator index, in order of computation plan
      */
      Pair<int> const & propagatorId(int id) const;

      /**
      * Get an id for a propagator from one vertex towards a target.
      *
      * For is != it, the return value is an identifier for an outgoing
      * propagator that begins at the source vertex (vertex index is)
      * and is the first edge of the directed path that leads from the
      * source vertex to the target vertex (vertex id it). The return
      * value is a pair of integers analogous to that returned by the
      * propagatorId(int) member function, for which the first element
      * is a block index and the second element is a direction id (0 or 1)
      * for the propagator direction for that block that is outgoing
      * from the source vertex.
      *
      * For the case is == it, the return value is a pair [-1, -1].
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
      * The algorithm creates a plan for computing propagators in an
      * that guarantees that the inital conditions required for each
      * propagator are known before it is processed. The algorithm is
      * works for any acyclic branched block copolymer. This function
      * is called in the default implementation of readParameters,
      * and must be called the readParameters method of any subclass.
      */
      virtual void makePlan();

      /**
      * Create a matrix of vertex-to-vertex path signposts.
      *
      * This function constructs the data structure that is accessed
      * by the paths member function.
      */
      void makePaths();

      /**
      * Check validity of the polymer graph.
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

      /// Number of vertices (ends and junctions) in this polymer
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
   inline int PolymerSpecies::nVertex() const
   {  return nVertex_; }

   /*
   * Number of blocks in this polymer.
   */
   inline int PolymerSpecies::nBlock() const
   {  return nBlock_; }

   /*
   * Number of propagators in this polymer.
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

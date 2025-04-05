#ifndef PSCF_POLYMER_SPECIES_H
#define PSCF_POLYMER_SPECIES_H

/*
* PSCF - Polymer Self-Consistent Field Theory
*
* Copyright 2016 - 2022, The Regents of the University of Minnesota
* Distributed under the terms of the GNU General Public License.
*/

#include <pscf/chem/Species.h>           // base class
#include <util/param/ParamComposite.h>   // base class

#include <pscf/chem/Vertex.h>            // member 
#include <pscf/chem/PolymerType.h>       // member
#include <pscf/chem/PolymerModel.h>      // member
#include <util/containers/Pair.h>        // member 
#include <util/containers/DArray.h>      // member 

#include <cmath>

namespace Pscf
{

   class Edge;
   using namespace Util;

   /**
   * Descriptor for an acyclic branched block polymer.
   *
   * A PolymerSpecies object has arrays of Block and Vertex objects.
   * Each Block has two propagator MDE solver objects.  The solve() 
   * member function solves the modified diffusion equation (MDE) for 
   * all propagators in molecule.
   *
   * \ingroup Pscf_Chem_Module
   */
   class PolymerSpecies : public Species, public ParamComposite
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
      * Read and initialize.
      *
      * \param in input parameter stream
      */
      virtual void readParameters(std::istream& in);

      /// \name Accessors (objects, by reference)
      ///@{

      /**
      * Get a specified Edge (block descriptor) by non-const reference.
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
      * Both chain ends and junctions are vertices.
      *
      * \param id vertex index, 0 <= id < nVertex
      */
      const Vertex& vertex(int id) const;

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
      * Allocate array of blocks.
      */
      virtual void allocateBlocks() = 0;

      /**
      * Read array of block data.
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
      * Construct the paths data structure.
      */
      void makePaths();

      /**
      * Check the validity of the polymer graph.
      */
      void isValid();

   private:

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
   inline int PolymerSpecies::nVertex() const
   {  return nVertex_; }

   /*
   * Number of blocks (or edges of the graph).
   */
   inline int PolymerSpecies::nBlock() const
   {  return nBlock_; }

   /*
   * Number of propagators.
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

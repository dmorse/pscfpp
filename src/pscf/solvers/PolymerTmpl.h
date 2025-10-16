#ifndef PSCF_POLYMER_TMPL_H
#define PSCF_POLYMER_TMPL_H

/*
* PSCF - Polymer Self-Consistent Field
*
* Copyright 2015 - 2025, The Regents of the University of Minnesota
* Distributed under the terms of the GNU General Public License.
*/

#include <pscf/chem/PolymerSpecies.h>    // base class
#include <util/param/ParamComposite.h>   // base class

#include <util/containers/Pair.h>        // member template
#include <util/containers/DArray.h>      // member template

namespace Pscf {

   class Edge;
   using namespace Util;

   /**
   * Template for an MDE solver and descriptor for a block polymer.
   *
   * Class template argument BT is an alias for a class that represents a 
   * block of a block polymer. By convention, this is a class named Block
   * defined in each implementation-level sub-namespace of Pscf.  Each 
   * such namespace contains a class named Block that is a subclass of 
   * Pscf::Edge, and a class named Polymer that is a subclass of the
   * template instantiation PolymerTmpl<Block>.
   * 
   * A PolymerTmpl<Block> object has an array of Block objects, as well as
   * an array of Vertex objects inherited from the PolymerSpecies base
   * class.  Each Block has two Propagator MDE solver objects associated
   * with the two directions along each block.
   *
   * The solve() member function solves the modified diffusion equation
   * (MDE) for all propagators in the molecule (i.e., all blocks, in both
   * directions) in a pre-defined order.
   *
   * \ingroup Pscf_Solver_Module
   */
   template <class BT, class PT>
   class PolymerTmpl : public PolymerSpecies
   {

   public:

      // Public typename aliases

      /**
      * Block of a block polymer.
      */
      using BlockT = BT;

      /**
      * Modified diffusion equation solver for one block, in one direction.
      */
      using PropagatorT = PT;

      // Public member functions

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
      * Solve modified diffusion equation for all propagators.
      *
      * Upon return, all propagators for all blocks, molecular partition
      * function q, phi and mu are all set.  The implementation of
      * PolymerTmpl::solve() the calls the solve() function for all
      * propagators in the molecule in a predeternined order, then
      * computes q, and finally computes mu from phi or phi from mu,
      * depending on the ensemble (open or closed) for this species.
      * This solve() function does not compute monomer concentration
      * field.
      *
      * Each implementation-level namespace defines a concrete subclass 
      * of PolymerTmpl<BT,PT> that is named Polymer by convention. Each such 
      * Polymer class defines a function named "compute" that takes an 
      * array of chemical fields (w-fields) as an argument, and that calls
      * PolymerTmpl<BT,PT>::solve internally.  Before calling the solve()
      * function declared here, the Polymer::compute() function must pass
      * pass the w-fields and any other required mutable data to all Block
      * objects in order to set up the MDE solver for each block. After
      * calling the solve() function, the Polymer::compute function must
      * also compute monomer concentrations for all blocks, each of which
      * is stored in a field container owned by the associated Block.
      *
      * The optional parameter phiTot is only relevant to problems 
      * involving a Mask that excludes material from part of the unit cell,
      * as done to define thin film problems.
      *
      * \param phiTot  fraction of unit cell volume occupied by material
      */
      virtual void solve(double phiTot = 1.0);

      /// \name Accessors (objects, by reference)
      ///@{

      /**
      * Get a specified Edge (block descriptor) by non-const reference.
      *
      * The edge member functions implements pure virtual functions
      * defined by the PolymerSpecies base class, and provide access to
      * each Block as a reference to an Edge (a block descriptor).
      *
      * \param id block index, 0 <= id < nBlock
      */
      Edge& edge(int id) final;

      /**
      * Get a specified Edge (block descriptor) by const reference.
      *
      * \param id block index, 0 <= id < nBlock
      */
      Edge const& edge(int id) const final;

      /**
      * Get a specified Block (solver and descriptor).
      *
      * \param id block index, 0 <= id < nBlock
      */
      BT& block(int id);

      /**
      * Get a specified Block (solver and descriptor) by const reference.
      *
      * \param id block index, 0 <= id < nBlock
      */
      BT const& block(int id) const ;

      /**
      * Get the propagator for a specific block and direction (non-const).
      *
      * For an edge that terminates at vertices with vertex indices given
      * by the return values of Edge::vertexId(0) and Edge::vertexId(1):
      *
      *    - direction 0 propagates from vertexId(0) to vertexId(1)
      *    - direction 1 propagates from vertexId(1) to vertexId(0)
      *
      * \param blockId integer index of associated block
      * \param directionId integer index for direction (0 or 1)
      */
      PT& propagator(int blockId, int directionId);

      /**
      * Get the propagator for a specific block and direction (const).
      *
      * \param blockId integer index of associated block
      * \param directionId integer index for direction (0 or 1)
      */
      PT const & propagator(int blockId, int directionId) const;

      /**
      * Get a propagator indexed in order of computation (non-const).
      *
      * The propagator index argument must satisfy 0 <= id < 2*nBlock.
      *
      * \param id  propagator index, in order of computation plan
      */
      PT& propagator(int id);

      ///@}

      // Inherited public members

      using PolymerSpecies::vertex;
      using PolymerSpecies::propagatorId;
      using PolymerSpecies::path;
      using PolymerSpecies::nBlock;
      using PolymerSpecies::nVertex;
      using PolymerSpecies::nPropagator;
      using PolymerSpecies::length;
      using PolymerSpecies::nBead;
      using PolymerSpecies::type;

      using Species::phi;
      using Species::mu;
      using Species::q;
      using Species::ensemble;

   protected:

      /**
      * Allocate array of Block objects.
      */
      void allocateBlocks() final;

      /**
      * Read array of data for blocks from parameter file
      *
      * \param in parameter input stream
      */
      void readBlocks(std::istream& in) final;

   private:

      /// Array of Block objects in this polymer.
      DArray<BT> blocks_;

      /**
      * Check validity of internal data structures set by readParameters.
      *
      * An Exception is thrown if any error is detected.
      */
      void isValid();

   };

   // Inline functions

   /*
   * Get a specified Edge (block descriptor) by non-const reference.
   */
   template <class BT, class PT>
   inline
   Edge& PolymerTmpl<BT,PT>::edge(int id)
   {  return blocks_[id]; }

   /*
   * Get a specified Edge (block descriptor) by const reference.
   */
   template <class BT, class PT>
   inline
   Edge const & PolymerTmpl<BT,PT>::edge(int id) const
   {  return blocks_[id]; }

   /*
   * Get a specified Block solver by non-const reference.
   */
   template <class BT, class PT>
   inline
   BT& PolymerTmpl<BT,PT>::block(int id)
   {  return blocks_[id]; }

   /*
   * Get a specified Block solver by const reference.
   */
   template <class BT, class PT>
   inline
   BT const & PolymerTmpl<BT,PT>::block(int id) const
   {  return blocks_[id]; }

}
#endif

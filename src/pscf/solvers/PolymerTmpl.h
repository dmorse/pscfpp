#ifndef PSCF_POLYMER_TMPL_H
#define PSCF_POLYMER_TMPL_H

/*
* PSCF - Polymer Self-Consistent Field
*
* Copyright 2015 - 2026, David Morse (morse012@um.edu) and co-authors
* Distributed under the terms of the GNU General Public License.
*/

#include <pscf/chem/PolymerSpecies.h>    // base class template
#include <util/containers/DArray.h>      // member

#include <functional>
#include <vector>

namespace Pscf {

   // Forward declaration
   class Edge;

   using namespace Util;

   /**
   * Descriptor and MDE solver for a block polymer.
   * 
   * Specializations of this template are used as base classes for classes
   * named "Polymer" that are defined in each program-level namespace. 
   * Each such program-level Polymer class is both a subclass of 
   * PolymerSpecies, which provides a description of polymer structure, 
   * and a modified diffusion equation solver. The PolymerTmpl template 
   * provides templates for several data structures and algorithms that 
   * are needed by all such Polymer classes. Specializations of this 
   * template should only be used within instances of derived classes,
   * The constructor and destructor are protected to document this intent.
   *
   * Class template argument BT is an alias for a class that represents 
   * a block of a block polymer. By convention, this is a class named 
   * Block defined in each program-level sub-namespace of Pscf.  
   *
   * Class template PT is an alias for the class that represents a
   * propagator, which holds the solution of the modified diffusion
   * equation for one block, in one direction. By convention, this is a
   * class named Propagator defined in each program-level sub-namespace
   * of Pscf.
   *
   * Class template WT is an alias for the data type of the value of 
   * a chemical potential field defined at a grid point, which can be
   * either a real (e.g., double) or complex data type. This type is 
   * also the type used to store several thermodynamic properties 
   * stored as members of the Species<WT> base class.
   * 
   * A PolymerTmpl<Block, Propagator, WT> object has an array of Block 
   * objects, as well as an array of Vertex objects inherited from the 
   * PolymerSpecies<WT> base class.  Each Block owns two Propagator 
   * objects associated with the two directions along each block.
   *
   * The solve() member function solves the modified diffusion equation
   * (MDE) for all propagators in the molecule (i.e., all blocks, in 
   * both directions) in a pre-defined order.
   *
   * \ingroup Pscf_Solver_Module
   */
   template <class BT, class PT, typename WT = double>
   class PolymerTmpl : public PolymerSpecies<WT>
   {

   public:

      // Public member functions

      /**
      * Read and initialize.
      *
      * \param in input parameter stream
      */
      void readParameters(std::istream& in) override;

      /**
      * Solve modified diffusion equation for all propagators.
      *
      * Upon return, solutions are stored for both propagators of all
      * blocks, and values are set for the molecular partition function
      * q, and phi or mu.  This function calls the solve() function for 
      * all propagators in the molecule in a predeternined order, then 
      * computes q, and finally computes mu from phi or phi from mu, 
      * depending on the ensemble (open or closed) for this polymer 
      * species. The order in which propgator solutions are computed is
      * determined by a plan that is defined by the PolymerSpecies base
      * class. This function does not compute monomer concentration 
      * fields.
      *
      * Each program-level namespace defines a class or class template
      * named Polymer that is derived from a specialization of PolymerTmpl.
      * Each such program-level Polymer class defines a member function
      * named "compute" that takes an array of chemical fields (w-fields) 
      * as an argument, and that calls PolymerTmpl<BT,PT,WT>::solve 
      * internally.  Before calling the PolymerTmpl::solve() function,
      * the Polymer::compute() function must also pass the w-fields and 
      * any other required mutable data to all Block objects in order to 
      * set up the MDE solver for each block. After calling this solve() 
      * function, the Polymer::compute function must also compute monomer
      * concentrations for all blocks, each of which is stored in a field 
      * container owned by the associated Block.
      *
      * The optional parameter phiTot is only relevant to problems 
      * involving a Mask that excludes material from part of the unit 
      * cell, as done to define thin film problems. In problems that do
      * not use such a mask, the default parameter value of phiTot = 1.0 
      * should be used.
      *
      * Each block solves the MDE in the field w - wShift of its monomer
      * type, where wShift = BT::wShift() is a constant chosen when the
      * solver is set up. The partition function of the propagators so
      * obtained, which is returned by qShifted(), then exceeds q by a
      * factor exp(S), where S is the sum over blocks of wShift times
      * the block length (thread model) or the number of beads in the
      * block (bead model). This function stores ln q = ln qShifted() - S,
      * from which mu or phi is computed, and never forms q itself.
      *
      * \param phiTot  fraction of unit cell volume occupied by material
      */
      virtual void solve(double phiTot = 1.0);

      /**
      * Get the partition function of the propagators as solved.
      *
      * This is the spatial average of the product of the propagators
      * at any point of the chain, computed with the shifted fields of
      * all blocks (see solve()). Products of propagators of this
      * polymer, e.g. in block concentrations and stresses, must be
      * normalized by this value rather than by q().
      */
      WT qShifted() const;

      /// \name Accessors (objects, by reference)
      ///@{

      /**
      * Get a specified Edge (block descriptor) by non-const reference.
      *
      * The edge member function implements a pure virtual function
      * defined by the PolymerSpecies base class, and provides access to
      * a specific Block as a reference to an Edge (a block descriptor),
      * which is a base class of the Block (BT) class.
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

      // Public typename aliases

      /**
      * Block of a block polymer.
      */
      using BlockT = BT;

      /**
      * Modified diffusion equation solution (propagator).
      */
      using PropagatorT = PT;

   protected:

      /**
      * Constructor.
      */
      PolymerTmpl();

      /**
      * Destructor.
      */
      ~PolymerTmpl() override = default;

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

      /**
      * Raise the field shifts of all blocks if qShifted() is too small.
      *
      * The field shifts chosen by the compute function of a Polymer
      * class keep the propagators bounded, so that qShifted() is at
      * most about one, but qShifted() can be very small, or zero, for
      * a long chain in a field that varies strongly over distances the
      * chain cannot follow. If qShifted() is below qMin, about 1e-294,
      * this function raises the shifts so that qShifted() becomes one,
      * without letting a propagator overflow. First each block is
      * calibrated alone: its shift is raised until the spatial average
      * of the tail of propagator(0) of the block, solved from a uniform
      * head, is one. Then the shifts of all blocks are raised by a common
      * constant, and the polymer is solved again, until qShifted() is
      * one. If that result is not finite or still below qMin (as for a
      * chain of many blocks confined to a small part of the cell, whose
      * blocks are calibrated alone with too large shifts), the shifts
      * given on entry are instead raised by a common constant only.
      *
      * Each raise of a shift by a constant c multiplies an average q of
      * the propagators concerned, of a total length (or number of
      * beads) l, by about exp(c*l). While q is below qMin, where it may
      * be zero or an inaccurate subnormal number, the shift is raised
      * in steps of c = (ln DBL_MAX - 32)/l, which cannot make q
      * overflow (at most 1000 steps). Then it is raised once by
      * c = -ln(q)/l.
      *
      * On entry, the propagators must be solved with the shifts given
      * in the array shifts. On return, they are solved with the shifts
      * that the array then contains, which the blocks also store.
      *
      * \param shifts  field shift of each block (input and output)
      * \param setupBlock  function(j, s) that sets up the MDE solver of
      *                    block j in the field w - s
      * \param solveBlock  function(j) that solves propagator(0) of block
      *                    j from a uniform head equal to one, and returns
      *                    the spatial average of its last slice
      * \param phiTot  fraction of unit cell volume occupied by material
      */
      void raiseShifts(std::vector<double> & shifts,
                       std::function<void(int, double)> const & setupBlock,
                       std::function<double(int)> const & solveBlock,
                       double phiTot = 1.0);

      // Protected typename aliases (base classes)

      /**
      * Indirect (grandparent) base class.
      */
      using SpeciesT = Species<WT>;

      /**
      * Direct (parent) base class.
      */
      using PolymerSpeciesT = PolymerSpecies<WT>;

   private:

      /// Array of Block objects in this polymer.
      DArray<BT> blocks_;

      /// Partition function of the propagators of the shifted fields.
      WT qShifted_;

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
   template <class BT, class PT, typename WT> inline
   Edge& PolymerTmpl<BT,PT,WT>::edge(int id)
   {  return blocks_[id]; }

   /*
   * Get a specified Edge (block descriptor) by const reference.
   */
   template <class BT, class PT, typename WT> inline
   Edge const & PolymerTmpl<BT,PT,WT>::edge(int id) const
   {  return blocks_[id]; }

   /*
   * Get the partition function of the propagators as solved.
   */
   template <class BT, class PT, typename WT> inline
   WT PolymerTmpl<BT,PT,WT>::qShifted() const
   {  return qShifted_; }

   /*
   * Get a specified Block solver by non-const reference.
   */
   template <class BT, class PT, typename WT> inline
   BT& PolymerTmpl<BT,PT,WT>::block(int id)
   {  return blocks_[id]; }

   /*
   * Get a specified Block solver by const reference.
   */
   template <class BT, class PT, typename WT> inline
   BT const & PolymerTmpl<BT,PT,WT>::block(int id) const
   {  return blocks_[id]; }

}
#endif

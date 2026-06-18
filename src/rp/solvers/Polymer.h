#ifndef RP_POLYMER_H
#define RP_POLYMER_H

/*
* PSCF - Polymer Self-Consistent Field
*
* Copyright 2015 - 2025, The Regents of the University of Minnesota
* Distributed under the terms of the GNU General Public License.
*/

#include <pscf/solvers/PolymerTmpl.h>    // base class template
#include <util/containers/FSArray.h>     // member template

// Forward declarations
namespace Util {
   template <typename T> class DArray;
}

namespace Pscf {
namespace Rp {

   using namespace Util;

   /**
   * Descriptor and MDE solver for one polymer species.
   *
   * \ref user_param_polymer_sec "Manual Page"
   *
   * This template adds functions required for periodic systems to an
   * appropriate specialization of the Pscf::PolymerTmpl base class 
   * template.
   *
   * Specializations of this class template are used as base classes for 
   * two closely analogous class templates, also named Polymer, that are
   * defined in the Rpc and Rpg program-level name spaces for use in the
   * pscf_rpc and pscf_rpg programs, respectively. 
   *
   * <b> Template parameters </b>:
   *
   *   - D  dimension of space
   *   - T  Types class (Rpc::Types<D> or Rpg::Types<D>)
   *
   * \ingroup Rp_Solver_Module
   */
   template <int D, class T>
   class Polymer 
     : public PolymerTmpl<typename T::Block, typename T::Propagator, double>
   {

   public:

      // Public type name aliases

      /// Block type, for a block within a block polymer.
      using BlockT = typename T::Block;

      /// Propagator type, for one direction within a block.
      using PropagatorT = typename T::Propagator;

      /// Direct base class, specialization of PolymerTmpl class template.
      using PolymerTmplT = PolymerTmpl<BlockT, PropagatorT, double>;

      /// Indirect base, specialization of PolymerSpecies (inherited).
      using typename PolymerTmplT::PolymerSpeciesT;

      /// Indirect base, specialization of Species (inherited).
      using typename PolymerTmplT::SpeciesT;

      // Public member functions

      /**
      * Constructor.
      */
      Polymer();

      /**
      * Destructor.
      */
      ~Polymer();

      /**
      * Set the number of unit cell parameters.
      *
      * \param nParam  the number of unit cell parameters
      */
      void setNParams(int nParam);

      /**
      * Clear all data that depends on unit cell parameters.
      *
      * This function should be called after each change in the unit cell.
      * It calls Block::clearUnitCellData() for all blocks in this
      * polymer.
      */
      void clearUnitCellData();

      /**
      * Compute MDE solutions and block concentrations.
      *
      * This function sets up w-fields in the MDE solvers for all blocks,
      * calls the base class PolymerTmpl solve function to solve the MDE
      * for all blocks, and then computes concentrations associated with
      * all blocks. On return, the associated Block objects all contain
      * propagator solutions and block volume fraction fields, while q and
      * phi or mu are set to new values.
      *
      * The parameter phiTot is only relevant to problems such as thin
      * films in which the material is excluded from part of the unit
      * cell by imposing an inhogeneous constraint on the sum of the
      * monomer concentrations (i.e., a "mask").
      *
      * \param wFields array of chemical potential fields.
      * \param phiTot  volume fraction of unit cell occupied by material
      */
      void compute(DArray<typename T::RField> const & wFields,
                   double phiTot = 1.0);

      /**
      * Compute SCFT stress contribution from this polymer species.
      *
      * This function computes contributions from this species to the
      * derivatives of SCFT free energy per monomer with respect to unit
      * cell parameters and stores the values. It requires that the MDE
      * has been solved for all blocks prior to entry, and so must be
      * called after the compute function.
      */
      void computeStress();

      /**
      * Get the precomputed contribution to stress from this species.
      *
      * Get the contribution from this polymer species to the derivative of
      * free energy per monomer with respect to unit cell parameter n, as
      * computed by the most recent call to computeStress().
      *
      * \param n  unit cell parameter index
      */
      double stress(int n) const;

      // Inherited non-dependent public member functions (for convenience)
      using PolymerTmplT::block;
      using PolymerSpeciesT::nBlock;
      using PolymerSpeciesT::length;
      using PolymerSpeciesT::nBead;
      using SpeciesT::phi;
      using SpeciesT::mu;
      using SpeciesT::q;
      using SpeciesT::ensemble;

   private:

      /// Stress contribution from this polymer species
      FSArray<double, 6> stress_;

      /// Number of unit cell parameters
      int nParam_;

      // Restrict access to inherited functions
      using PolymerTmplT::solve;
      using SpeciesT::setQ;

   };

   /// Get stress component n.
   template <int D, class T> inline
   double Polymer<D,T>::stress(int n) const
   {  return stress_[n]; }

} // namespace Rp
} // namespace Pscf
#endif

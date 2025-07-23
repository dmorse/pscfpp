#ifndef RPC_POLYMER_H
#define RPC_POLYMER_H

/*
* PSCF - Polymer Self-Consistent Field Theory
*
* Copyright 2015 - 2025, The Regents of the University of Minnesota
* Distributed under the terms of the GNU General Public License.
*/

#include "Block.h"                       // base class template parameter
#include <pscf/solvers/PolymerTmpl.h>    // base class template
#include <util/containers/FSArray.h>     // member template

// Forward declarations
namespace Util {
   template <typename T> class DArray;
}
namespace Pscf { 
   namespace Prdc { 
      namespace Cpu { 
         template <int D> class RField;
      }
   }
}

namespace Pscf { 
namespace Rpc { 

   using namespace Prdc;
   using namespace Pscf::Prdc::Cpu;

   /**
   * Descriptor and solver for one polymer species.
   *
   * The phi() and mu() accessor functions, which are inherited from
   * PolymerSpecies, return the value of phi (spatial average volume
   * fraction of a species) or mu (species chemical potential) computed
   * in the last call of the compute() function. If the ensemble for this
   * species is closed, phi is read from the parameter file and mu is 
   * computed. If the ensemble is open, mu is read from the parameter 
   * file and phi is computed.
   *
   * The block concentrations stored in the constituent Block<D> objects
   * contain the block concentrations (i.e., volume fractions) computed in
   * the most recent call of the compute function. These can be accessed
   * using the Block<D>::cField() function.
   *
   * \ref user_param_polymer_sec "Manual Page"
   *
   * \ingroup Rpc_Solver_Module
   */
   template <int D>
   class Polymer : public PolymerTmpl< Block<D> >
   {

   public:

      // Public type name aliases

      /// Base class, partial template specialization.
      using Base = PolymerTmpl< Block<D> >;

      /// Block type, for a block within a block polymer.
      using BlockT = Block<D>;

      /// Propagator type, for one direction within a block.
      using PropagatorT = typename BlockT::PropagatorT;

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
      * It calls Block<D>::clearUnitCellData() for all blocks in this 
      * polymer.
      */ 
      void clearUnitCellData();

      /**
      * Compute solution to MDE and block concentrations.
      * 
      * This function sets up w-fields in the MDE solvers for all blocks,
      * calls the base class PolymerTmpl solve function to solve the MDE
      * for all blocks, and then computes concentrations associated with 
      * all blocks. On return, the associated Block objects all contain
      * propagator solutions and block volume fraction fields, while q and
      * phi or mu are set to new values.
      *
      * The parameter phiTot is only relevant to problems such as thin
      * films in which the material is excluded from part of the unit cell
      * by imposing an inhogeneous constraint on the sum of the monomer
      * concentrations (i.e., a "Mask"). 
      *
      * \param wFields array of chemical potential fields.
      * \param phiTot  volume fraction of unit cell occupied by material
      */ 
      void compute(DArray< RField<D> > const & wFields, 
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
      * Get precomputed contribution to stress from this species.
      *  
      * Get the contribution from this polymer species to the deritative of
      * free energy per monomer with respect to unit cell parameter n, as
      * computed by the most recent call to computeStress().
      *  
      * \param n index of unit cell parameter
      */
      double stress(int n) const;

      // Inherited public member functions
      
      using Base::edge;
      using Base::block;
      using Base::propagator;
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
      using Species::setPhi;
      using Species::setMu;

   private: 

      /// Stress contribution from this polymer species
      FSArray<double, 6> stress_;

      /// Number of unit cell parameters.
      int nParam_;

      // Restricting access to inherited functions
      using Base::solve;
      using Species::setQ;

   };

   /// Get stress component n.
   template <int D>
   inline double Polymer<D>::stress(int n) const
   {  return stress_[n]; }
  
   #ifndef RPC_POLYMER_TPP
   // Supress implicit instantiation
   extern template class Polymer<1>;
   extern template class Polymer<2>;
   extern template class Polymer<3>;
   #endif

}
}
#endif

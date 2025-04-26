#ifndef RPC_SIM_STATE_H
#define RPC_SIM_STATE_H

/*
* PSCF - Polymer Self-Consistent Field Theory
*
* Copyright 2015 - 2025, The Regents of the University of Minnesota
* Distributed under the terms of the GNU General Public License.
*/


#include <prdc/cpu/RField.h>
#include <util/containers/DArray.h>

namespace Pscf {
namespace Rpc {

   using namespace Util;
   using namespace Prdc::Cpu;

   /**
   * SimState stores the state used by an FTS simulation.
   *
   * This class is used to restore the state of FTS simulation after an
   * attempted move or step that is rejected or fails to converge. It is
   * used in Monte Carlo (MC) simulations to restore the state after a 
   * rejected move. It is also used less frequently in Brownian dynamics 
   * (BD) simulations to restore the previous state after the compressor 
   * algorithm (the search for a partial saddle point) fails to converge 
   * after an attempted unconstrained BD step.
   *
   * \ingroup Rpc_Fts_Module
   */
   template <int D>
   struct SimState 
   {

   public:

      // Public member functions

      /**
      * Constructor.
      */
      SimState();

      /**
      * Destructor.
      */
      ~SimState();

      /**
      * Allocate memory for fields.
      *
      * \param nMonomer  number of monomer types
      * \param dimensions  dimensions of discretization grid
      */ 
      void allocate(int nMonomer, IntVec<D> const & dimensions);
 
      // Public data members

      /**
      * Chemical potential fields, r-grid format, indexed by monomer.
      *
      * Field w[i] is the chemical potential field for monomer type i,
      * for i = 0, ..., nMonomer - 1.
      */
      DArray< RField<D> > w;

      /**
      * Chemical potential fields, r-grid format, indexed by eigenvector.
      *
      * Field wc[i] is a pointwise projection of the w fields onto 
      * eigenvector number i of the projected chi matrix. for values
      * i = 0, ..., nMonomer - 1.
      */
      DArray< RField<D> > wc;
      
      /**
      * Eigenvector components of c fields on a real space grid.
      *
      * Field cc[i] is a point-wise projection of the c fields onto
      * eigenvector number i of the projected chi matrix , for values
      * i = 0, ..., nMonomer - 1.
      */
      DArray< RField<D> > cc;
      
      /**
      * Functional derivatives of the Hamiltonian on a real space grid.
      *
      * Field dc[i] is the functional derivative of H[W] with respect to
      * w-field component wc[i], indexed by eigenvector index i.
      */
      DArray< RField<D> > dc;
      
      /// Field theoretic Hamiltonian value (total).
      double hamiltonian;
      
      /// Ideal gas contribution to Hamiltonian.
      double idealHamiltonian;
      
      /// Quadratic field contribution to Hamiltonian value.
      double fieldHamiltonian;
            
      /// Perturbation to Hamiltonian value (if any).
      double perturbationHamiltonian;
            
      /// True iff cc fields need to be saved.
      bool needsCc;
      
      /// True iff dc fields need to be saved.
      bool needsDc;
      
      /// True iff Hamiltonian components need to be saved.
      bool needsHamiltonian;

      /// Does this object currently store data?
      bool hasData;
      
      /// Has memory been allocated for the fields?
      bool isAllocated;

   };

   #ifndef RPC_SIM_STATE_TPP
   // Suppress implicit instantiation
   extern template struct SimState<1>;
   extern template struct SimState<2>;
   extern template struct SimState<3>;
   #endif

}
}
#endif

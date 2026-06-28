#ifndef RP_SIM_STATE_H
#define RP_SIM_STATE_H

/*
* PSCF - Polymer Self-Consistent Field
*
* Copyright 2015 - 2025, The Regents of the University of Minnesota
* Distributed under the terms of the GNU General Public License.
*/

#include <pscf/math/IntVec.h>         // template with default parameter
#include <util/containers/DArray.h>   // member

namespace Pscf {
namespace Rp {

   using namespace Util;

   /**
   * SimState stores the state used by a FTS simulation.
   *
   * This class is used to restore the state of FTS simulation after an
   * attempted MC move or BD step is rejected or fails to converge. It
   * is used in Monte Carlo (MC) simulations to restore the state after a 
   * rejected move. It is also used less frequently in Brownian dynamics 
   * (BD) simulations to restore the previous state if the compressor 
   * algorithm (the search for a partial saddle point) fails to converge 
   * after an attempted unconstrained BD step.
   * 
   * Specializations of this class template are used as base classes for 
   * two closely analogous class templates, also named SimState, that
   * are defined in Rpc and Rpg namespaces and used in the pscf_rpc and
   * pscf_rpg programs, respectively.
   *
   * Template parameters:
   *   - D  : dimension of space
   *   - FT : field type, Rpc::RField<D> or Rpg::RField<D> 
   *
   * \ingroup Rp_Fts_Simulator_Module
   */
   template <int D, class T>
   class SimState 
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
      ~SimState() = default;

      /**
      * Allocate memory for stored fields.
      *
      * \param nMonomer  number of monomer types
      * \param dimensions  dimensions of discretization grid
      */ 
      void allocate(int nMonomer, IntVec<D> const & dimensions);

      // Public data members
      
      using FT = typename T::RField;

      /**
      * Chemical potential fields, r-grid format, indexed by monomer.
      *
      * Field w[i] is the chemical potential field for monomer type i,
      * for i = 0, ..., nMonomer - 1.
      */
      DArray<FT> w;

      /**
      * Chemical potential fields, r-grid format, indexed by eigenvector.
      *
      * Field wc[i] is a pointwise projection of the w fields onto 
      * eigenvector number i of the projected chi matrix. for values
      * i = 0, ..., nMonomer - 1.
      */
      DArray<FT> wc;
      
      /**
      * Eigenvector components of c fields on a real space grid.
      *
      * Field cc[i] is a point-wise projection of the c fields onto
      * eigenvector number i of the projected chi matrix , for values
      * i = 0, ..., nMonomer - 1.
      */
      DArray<FT> cc;
      
      /**
      * Functional derivatives of the Hamiltonian on a real space grid.
      *
      * Field dc[i] is the functional derivative of H[W] with respect to
      * w-field component wc[i], indexed by eigenvector index i.
      */
      DArray<FT> dc;
      
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

}
}
#endif

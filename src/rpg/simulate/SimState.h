#ifndef RPG_SIM_STATE_H
#define RPG_SIM_STATE_H

/*
* PSCF - Polymer Self-Consistent Field Theory
*
* Copyright 2016 - 2022, The Regents of the University of Minnesota
* Distributed under the terms of the GNU General Public License.
*/


#include <prdc/cuda/RField.h>              // memmber 
#include <util/containers/DArray.h>

namespace Pscf {
namespace Rpg {

   using namespace Util;
   using namespace Pscf::Prdc::Cuda;

   /**
   * SimState stores the state used by an fts simulation.
   *
   * \ingroup Rpg_Simulate_Module
   */
   template <int D>
   struct SimState 
   {
      public:

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
 
      /**
      * Chemical potential fields, r-grid format, indexed by monomer.
      */
      DArray< RField<D> > w;

      /**
      * Chemical potential fields, r-grid format, indexed by eigenvector.
      *
      * Each field is a component projected on pointwise onto a
      * eigenvector of the projected chi matrix, with indices that
      * correspond to eigenvector indices.
      */
      DArray< RField<D> > wc;
      
      /**
      * Eigenvector components of c fields on a real space grid.
      *
      * Each field component corresponds to a point-wise projection of c
      * onto an eigenvector of the projected chi matrix.
      */
      DArray< RField<D> > cc;
      
      /**
      * Components of functional derivatives of the Hamiltonian fields 
      * on a real space grid.
      *
      * Each field component is the functional derivative of H[W]
      * with respect to one eigenvector w-field component.
      */
      DArray< RField<D> > dc;
      
      /// Monte-Carlo Hamiltonian value.
      double hamiltonian;
      
      /// Monte-Carlo ideal gas contribution to Hamiltonian value.
      double idealHamiltonian;
      
      /// Monte-Carlo field part contribution to Hamiltonian value.
      double fieldHamiltonian;
            
      /// If cc fields needs to be saved.
      bool needsCc;
      
      /// If dc fields needs to be saved.
      bool needsDc;
      
      /// If hamiltonian needs to be saved.
      bool needsHamiltonian;

      /// Is this struct being used to store data?
      bool hasData;
      
      /// Has memory be allocated for the w field?
      bool isAllocated;

   };

   #ifndef RPG_SIM_STATE_TPP
   // Suppress implicit instantiation
   extern template struct SimState<1>;
   extern template struct SimState<2>;
   extern template struct SimState<3>;
   #endif

}
}
#endif

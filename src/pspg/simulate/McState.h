#ifndef PSPG_MC_STATE_H
#define PSPG_MC_STATE_H

/*
* PSCF - Polymer Self-Consistent Field Theory
*
* Copyright 2016 - 2022, The Regents of the University of Minnesota
* Distributed under the terms of the GNU General Public License.
*/


#include <pspg/field/RDField.h>
#include <pspg/field/DField.h>  
#include <util/containers/DArray.h>

namespace Pscf {
namespace Pspg
{

   using namespace Util;

   /**
   * McState stores the state used by an MC simulation.
   *
   * \ingroup Pspg_Simulate_Module
   */
   template <int D>
   struct McState 
   {
      public:

      /**
      * Constructor.
      */
      McState();

      /**
      * Allocate memory for w fields.
      *
      * \param nMonomer  number of monomer types
      * \param meshSize  number of grid points
      */ 
      void allocate(int nMonomer, int meshSize);
 
      /**
      * Chemical potential fields, r-grid format, indexed by monomer.
      */
      DArray< RDField<D> > w;

      /**
      * Chemical potential fields, r-grid format, indexed by eigenvector.
      *
      * Each field is a component projected on pointwise onto a
      * eigenvector of the projected chi matrix, with indices that
      * correspond to eigenvector indices.
      */
      DArray< RDField<D> > wc;

      /// Monte-Carlo Hamiltonian value.
      double mcHamiltonian;
      
      /// Monte-Carlo ideal gas contribution to Hamiltonian value.
      double mcIdealHamiltonian;
      
      /// Monte-Carlo field part contribution to Hamiltonian value.
      double mcFieldHamiltonian;

      /// Is this struct being used to store data?
      bool hasData;

      /// Has memory be allocated for the w field?
      bool isAllocated;

   };

   #ifndef PSPG_MC_STATE_TPP
   // Suppress implicit instantiation
   extern template struct McState<1>;
   extern template struct McState<2>;
   extern template struct McState<3>;
   #endif

}
}
#endif

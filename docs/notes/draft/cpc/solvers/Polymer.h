#ifndef CPC_POLYMER_H
#define CPC_POLYMER_H

/*
* PSCF - Polymer Self-Consistent Field
*
* Copyright 2015 - 2025, The Regents of the University of Minnesota
* Distributed under the terms of the GNU General Public License.
*/

#include <pscf/backends/CPT.h>            // backend class
#include <pscf/solvers/PolymerTmpl.h>  // base class template

// Forward declarations
namespace Util {
   template <typename T> class DArray;
}
namespace Pscf { 
   namespace Cpc {
      template <int D> class Block;
      template <int D> class Propagator;
   }
   namespace Prdc { 
      template <int D, class T> class CField;
   }
}

// Explicit instantiation declarations for base class
namespace Pscf { 
   extern template 
   class PolymerTmpl< Cpc::Block<1>, Cpc::Propagator<1>, std::complex<double> >;
   extern template 
   class PolymerTmpl< Cpc::Block<2>, Cpc::Propagator<2>, std::complex<double> >;
   extern template 
   class PolymerTmpl< Cpc::Block<3>, Cpc::Propagator<3>, std::complex<double> >;
}

namespace Pscf {
namespace Cpc {

   using namespace Prdc;

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
   * \ingroup Cpc_Solver_Module
   */
   template <int D>
   class Polymer 
    : public PolymerTmpl< Block<D>, Propagator<D>, std::complex<double> >
   {

   public:

      // Public type name aliases

      /// Base class, partial template specialization.
      using Base = PolymerTmpl< Block<D>, Propagator<D>, std::complex<double> >;

      /// PolymerSpecies template base class.
      using PolymerSpeciesT = PolymerSpecies< std::complex<double> >;

      /// Species template base class.
      using SpeciesT = Species< std::complex<double> >;

      /// Block type, for a block within a block polymer.
      using BlockT = Block<D>;

      /// Propagator type, for one direction within a block.
      using PropagatorT = Propagator<D>;

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
      * \param wFields array of chemical potential fields.
      */ 
      void compute(DArray< CField<D,CPT> > const & wFields);

      // Inherited public member functions
      
      using Base::edge;
      using Base::block;
      using Base::propagator;
      using PolymerSpeciesT::vertex;
      using PolymerSpeciesT::propagatorId;
      using PolymerSpeciesT::path;
      using PolymerSpeciesT::nBlock;
      using PolymerSpeciesT::nVertex;
      using PolymerSpeciesT::nPropagator;
      using PolymerSpeciesT::length;
      using PolymerSpeciesT::nBead;
      using PolymerSpeciesT::type;
      using SpeciesT::phi;
      using SpeciesT::mu;
      using SpeciesT::q;
      using SpeciesT::ensemble;
      using SpeciesT::setPhi;
      using SpeciesT::setMu;

   private: 

      // Restricting access to inherited functions
      using Base::solve;
      using SpeciesT::setQ;

   };

   // Explicit instantiation declarations
   extern template class Polymer<1>;
   extern template class Polymer<2>;
   extern template class Polymer<3>;

}
}
#endif

#ifndef RPC_BINARY_STRUCTURE_FACTOR_H
#define RPC_BINARY_STRUCTURE_FACTOR_H

/*
* PSCF - Polymer Self-Consistent Field
*
* Copyright 2015 - 2025, The Regents of the University of Minnesota
* Distributed under the terms of the GNU General Public License.
*/

#include <rp/fts/analyzer/BinaryStructureFactorBase.h> // base template
#include <rpc/system/Types.h>                          // base argument
#include <rpc/fts/analyzer/Analyzer.h>                 // indirect base
#include <prdc/field/cpu/RField.h>                     // base member
#include <prdc/field/cpu/RFieldDft.h>                  // base member

namespace Pscf {
namespace Rpc {

   using namespace Util;

   /**
   * Spherically averaged structure factor for a two-monomer system.
   *
   * Specializations of this template are derived from specializations of 
   * the base class template BinaryStructureFactorBase, and inherit most
   * of their source code from this base class.
   *
   * \see BinaryStructureFactorBase
   * \see \ref rp_BinaryStructureFactor_page "Manual Page"
   * \ingroup Rpc_Fts_Analyzer_Module
   */
   template <int D>
   class BinaryStructureFactor 
    : public Rp::BinaryStructureFactorBase<D, Types<D> >
   {

   public:

      /**
      * Constructor.
      *
      * \param simulator  parent Simulator object
      * \param system  parent System object
      */
      BinaryStructureFactor(Rp::Simulator<D, Rpc::Types<D> >& simulator, 
                            Rp::System<D, Rpc::Types<D> >& system);

      /**
      * Setup before the main loop.
      */
      void setup() override;

      /**
      * Compute structure factors and add to accumulators.
      *
      * \param iStep step counter
      */
      void sample(long iStep) override;

   protected:

      using BinaryStructureFactorBase 
         = typename Rp::BinaryStructureFactorBase<D, Types<D> >;
      using typename BinaryStructureFactorBase::AnalyzerT;

      using BinaryStructureFactorBase::wk_;
      using BinaryStructureFactorBase::allocate;
      using BinaryStructureFactorBase::findWaveBunches;
      using BinaryStructureFactorBase::computeW;
      using BinaryStructureFactorBase::computeS;

   };

}
}

// Explicit instantiation declarations
namespace Pscf {
   namespace Rp {
      extern template class BinaryStructureFactorBase<1, Rpc::Types<1> >;
      extern template class BinaryStructureFactorBase<2, Rpc::Types<2> >;
      extern template class BinaryStructureFactorBase<3, Rpc::Types<3> >;
   }
   namespace Rpc {
      extern template class BinaryStructureFactor<1>;
      extern template class BinaryStructureFactor<2>;
      extern template class BinaryStructureFactor<3>;
   }
}
#endif

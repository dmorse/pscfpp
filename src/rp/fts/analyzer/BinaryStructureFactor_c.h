#ifndef RPC_BINARY_STRUCTURE_FACTOR_H
#define RPC_BINARY_STRUCTURE_FACTOR_H

/*
* PSCF - Polymer Self-Consistent Field
*
* Copyright 2015 - 2025, The Regents of the University of Minnesota
* Distributed under the terms of the GNU General Public License.
*/

#include <rp/fts/analyzer/BinaryStructureFactorBase.h> // base template
#include <pscf/backends/CPT.h>                         // base argument

namespace Pscf {
namespace Rp {

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
   * \ingroup Rp_Fts_Analyzer_Module
   */
   template <int D>
   class BinaryStructureFactor<D,CPT>
    : public BinaryStructureFactorBase<D,CPT>
   {

   public:

      /**
      * Constructor.
      *
      * \param simulator  parent Simulator object
      * \param system  parent System object
      */
      BinaryStructureFactor(Simulator<D,CPT>& simulator, 
                            System<D,CPT>& system);

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

      using Base = BinaryStructureFactorBase<D,CPT>;
      using AnalyzerT = Analyzer<D,CPT>;

      using Base::wk_;
      using Base::allocate;
      using Base::findWaveBunches;
      using Base::computeW;
      using Base::computeS;

   };

}
}

// Explicit instantiation declarations
namespace Pscf {
   namespace Rp {
      extern template class BinaryStructureFactor<1,CPT>;
      extern template class BinaryStructureFactor<2,CPT>;
      extern template class BinaryStructureFactor<3,CPT>;
   }
}
#endif

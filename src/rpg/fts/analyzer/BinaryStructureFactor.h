#ifndef RPG_BINARY_STRUCTURE_FACTOR_H
#define RPG_BINARY_STRUCTURE_FACTOR_H

/*
* PSCF - Polymer Self-Consistent Field
*
* Copyright 2015 - 2025, The Regents of the University of Minnesota
* Distributed under the terms of the GNU General Public License.
*/

#include <rp/fts/analyzer/BinaryStructureFactorBase.h> // base template
#include <rpg/system/Types.h>                          // base argument
#include <rpg/fts/analyzer/Analyzer.h>                 // indirect base
#include <pscf/cuda/HostDArray.h>                      // member
#include <prdc/field/cuda/RField.h>                    // base class member
#include <prdc/field/cuda/RFieldDft.h>                 // base class member

namespace Pscf {
namespace Rpg {

   using namespace Util;

   /**
   * Spherically averaged structure factor for a two-monomer system.
   *
   * Specializations of this template are derived from specializations of 
   * the base class template Rp::BinaryStructureFactor, and inherit most
   * of their source code from this base class.
   *
   * \see Rp::BinaryStructureFactor
   * \see \ref rp_BinaryStructureFactor_page "Manual Page"
   * \ingroup Rpg_Fts_Analyzer_Module
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
      BinaryStructureFactor(
         Rp::Simulator<D, Rpg::Types<D> >& simulator, 
         Rp::System<D, Rpg::Types<D> >& system);

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

      using BinaryStructureFactorBase::allocate;
      using BinaryStructureFactorBase::findWaveBunches;
      using BinaryStructureFactorBase::computeW;
      using BinaryStructureFactorBase::computeS;

      using BinaryStructureFactorBase::wk_;

   private:

      // Copy of wk_ on host CPU
      HostDArray<cudaComplex> wkHost_;

   };

}
}

// Explicit instantiation declarations
namespace Pscf {
   namespace Rp {
      extern template class BinaryStructureFactorBase<1, Rpg::Types<1> >;
      extern template class BinaryStructureFactorBase<2, Rpg::Types<2> >;
      extern template class BinaryStructureFactorBase<3, Rpg::Types<3> >;
   }
   namespace Rpg {
      extern template class BinaryStructureFactor<1>;
      extern template class BinaryStructureFactor<2>;
      extern template class BinaryStructureFactor<3>;
   }
}
#endif

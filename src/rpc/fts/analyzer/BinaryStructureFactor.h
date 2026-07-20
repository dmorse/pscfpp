#ifndef RPC_BINARY_STRUCTURE_FACTOR_H
#define RPC_BINARY_STRUCTURE_FACTOR_H

/*
* PSCF - Polymer Self-Consistent Field
*
* Copyright 2015 - 2025, The Regents of the University of Minnesota
* Distributed under the terms of the GNU General Public License.
*/

#include <rp/fts/analyzer/BinaryStructureFactorBase.h> // base template
#include <pscf/cpu/CppTp.h>                          // base argument
#include <rpc/fts/analyzer/Analyzer.h>                 // indirect base
#include <prdc/field/cpu/RField.h>                     // base member
#include <prdc/field/cpu/RFieldDft.h>                  // base member

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
   class BinaryStructureFactor<D, CppTp<D> >
    : public BinaryStructureFactorBase<D, CppTp<D> >
   {

   public:

      /**
      * Constructor.
      *
      * \param simulator  parent Simulator object
      * \param system  parent System object
      */
      BinaryStructureFactor(Simulator<D, CppTp<D> >& simulator, 
                            System<D, CppTp<D> >& system);

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

      using Base = BinaryStructureFactorBase<D, CppTp<D> >;
      using AnalyzerT = Analyzer<D, CppTp<D> >;

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
      extern template class BinaryStructureFactorBase<1, CppTp<1> >;
      extern template class BinaryStructureFactorBase<2, CppTp<2> >;
      extern template class BinaryStructureFactorBase<3, CppTp<3> >;
      extern template class BinaryStructureFactor<1, CppTp<1> >;
      extern template class BinaryStructureFactor<2, CppTp<2> >;
      extern template class BinaryStructureFactor<3, CppTp<3> >;
   }
}
#endif

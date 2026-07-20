#ifndef RPG_BINARY_STRUCTURE_FACTOR_H
#define RPG_BINARY_STRUCTURE_FACTOR_H

/*
* PSCF - Polymer Self-Consistent Field
*
* Copyright 2015 - 2025, The Regents of the University of Minnesota
* Distributed under the terms of the GNU General Public License.
*/

#include <rp/fts/analyzer/BinaryStructureFactorBase.h> // base template
#include <pscf/cuda/CudaTp.h>                          // base argument
#include <rpg/fts/analyzer/Analyzer.h>                 // indirect base
#include <pscf/cuda/HostDArray.h>                      // member
#include <prdc/field/cuda/RField.h>                    // base class member
#include <prdc/field/cuda/RFieldDft.h>                 // base class member

namespace Pscf {
namespace Rp {

   using namespace Util;

   /**
   * Spherically averaged structure factor for a two-monomer system.
   *
   * Specializations of this template are derived from specializations of 
   * the base class template BinaryStructureFactor, and inherit most
   * of their source code from this base class.
   *
   * \see BinaryStructureFactorBase
   * \see \ref rp_BinaryStructureFactor_page "Manual Page"
   * \ingroup Rp_Fts_Analyzer_Module
   */
   template <int D>
   class BinaryStructureFactor<D, CudaTp<D> > 
    : public BinaryStructureFactorBase<D, CudaTp<D> >
   {

   public:

      /**
      * Constructor.
      *
      * \param simulator  parent Simulator object
      * \param system  parent System object
      */
      BinaryStructureFactor(
         Simulator<D, CudaTp<D> >& simulator, 
         System<D, CudaTp<D> >& system);

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

      using Base = BinaryStructureFactorBase<D, CudaTp<D> >;
      using AnalyzerT = Analyzer<D, CudaTp<D> > ;

      using Base::allocate;
      using Base::findWaveBunches;
      using Base::computeW;
      using Base::computeS;

      using Base::wk_;

   private:

      // Copy of wk_ on host CPU
      HostDArray<cudaComplex> wkHost_;

   };

}
}

// Explicit instantiation declarations
namespace Pscf {
   namespace Rp {
      extern template class BinaryStructureFactorBase<1, CudaTp<1> >;
      extern template class BinaryStructureFactorBase<2, CudaTp<2> >;
      extern template class BinaryStructureFactorBase<3, CudaTp<3> >;
      extern template class BinaryStructureFactor<1, CudaTp<1> >;
      extern template class BinaryStructureFactor<2, CudaTp<2> >;
      extern template class BinaryStructureFactor<3, CudaTp<3> >;
   }
}
#endif

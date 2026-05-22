#ifndef RPG_BINARY_STRUCTURE_FACTOR_H
#define RPG_BINARY_STRUCTURE_FACTOR_H

/*
* PSCF - Polymer Self-Consistent Field
*
* Copyright 2015 - 2025, The Regents of the University of Minnesota
* Distributed under the terms of the GNU General Public License.
*/

#include <rp/fts/analyzer/BinaryStructureFactor.h>  // base template
#include <rpg/system/Types.h>                       // template argument
#include <rpg/fts/analyzer/Analyzer.h>              // indirect base
#include <pscf/cuda/HostDArray.h>                   // member
#include <prdc/cuda/RField.h>                       // base class member
#include <prdc/cuda/RFieldDft.h>                    // base class member

namespace Pscf {
namespace Rpg {

   template <int D> class System;
   template <int D> class Simulator;

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
    : public Rp::BinaryStructureFactor<D, Types<D> >
   {

   public:

      /**
      * Constructor.
      *
      * \param simulator  parent Simulator object
      * \param system  parent System object
      */
      BinaryStructureFactor(Simulator<D>& simulator, System<D>& system);

      /*
      * Setup before main simulation loop
      */
      void setup() override;

      /**
      * Compute structure factors and add to accumulators.
      *
      * \param iStep step counter
      */
      void sample(long iStep) override;

   protected:

      using RpBinaryStructureFactor 
                      = typename Rp::BinaryStructureFactor<D, Types<D> >;
      using typename RpBinaryStructureFactor::AnalyzerT;

      using RpBinaryStructureFactor::wk_;
      using RpBinaryStructureFactor::computeW;
      using RpBinaryStructureFactor::computeS;

   private:

      // Copy of wk_ on host CPU
      HostDArray<cudaComplex> wkHost_

   };

}
}

// Explicit instantiation declarations
namespace Pscf {
   namespace Rp {
      extern template class BinaryStructureFactor<1, Rpg::Types<1> >;
      extern template class BinaryStructureFactor<2, Rpg::Types<2> >;
      extern template class BinaryStructureFactor<3, Rpg::Types<3> >;
   }
   namespace Rpg {
      extern template class BinaryStructureFactor<1>;
      extern template class BinaryStructureFactor<2>;
      extern template class BinaryStructureFactor<3>;
   }
}
#endif

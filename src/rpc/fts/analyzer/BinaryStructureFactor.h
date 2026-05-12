#ifndef RPC_BINARY_STRUCTURE_FACTOR_H
#define RPC_BINARY_STRUCTURE_FACTOR_H

/*
* PSCF - Polymer Self-Consistent Field
*
* Copyright 2015 - 2025, The Regents of the University of Minnesota
* Distributed under the terms of the GNU General Public License.
*/

#include "Analyzer.h"                             // base class

#include <prdc/cpu/RField.h>                      // member
#include <prdc/cpu/RFieldDft.h>                   // member
#include <util/accumulators/Average.h>            // member
#include <util/containers/DArray.h>               // member

#include <string>
#include <iostream>
#include <fstream>

namespace Pscf {
namespace Rpc {

   template <int D> class System;
   template <int D> class Simulator;

   using namespace Util;
   using namespace Pscf::Prdc::Cpu;

   /**
   * Spherically averaged structure factor for a two-monomer system.
   *
   * This class evaluates the structures factors for all wavevectors 
   * within the Fourier space grid, while grouping averages for 
   * wavevectors of equal norm.
   *
   * A structure factor for a wavevector k for a system with two types
   * monomer is given by an expectation value
   * \f[
   *     S(k) = 1/(v \chi )^2 <W_(k)W_(-k)>/V - 1/(2 \chi v)
   * \f]
   * where, V is system volume, v is monomer volume, and \f$W_(k)\f$ is 
   * a Fourier transform of the fluctuating exchange field 
   * \f$ W_{-} = (w_{A} - w_{B})/2 \f$.
   *
   * \see \ref rp_BinaryStructureFactor_page "Manual Page"
   * \ingroup Rpc_Fts_Analyzer_Module
   */
   template <int D>
   class BinaryStructureFactor : public Analyzer<D>
   {

   public:

      /**
      * Constructor.
      *
      * \param simulator  parent simulator object
      * \param system  parent system object
      */
      BinaryStructureFactor(Simulator<D>& simulator, System<D>& system);

      /**
      * Read parameters from file.
      *
      * Input format:
      *
      *   - int        interval          sampling interval
      *   - string     outputFileName    output file base name
      *   - int        nSamplePerBlock   output file base name
      *
      * \param in input parameter stream
      */
      void readParameters(std::istream& in);

      /**
      * Setup immediately before a simulation.
      */
      void setup();

      /**
      * Compute structure factors and add to accumulators.
      *
      * \param iStep step counter
      */
      void sample(long iStep);

      /**
      * Output results to predefined output file.
      */
      void output();

   protected:

      // Alias for base class.
      using AnalyzerT = Analyzer<D>;

      // Inherited protected member functions (selected).
      using AnalyzerT::system;
      using AnalyzerT::simulator;

   private:

      /// Exchange field W_(r):  wm = (wa-wb)/2 ].
      RField<D> wm_;

      /// Discrete Fourier transform (DFT) of wm_  .
      RFieldDft<D> wk_;

      /// Bunch ids for waves, indexed by wave id.
      DArray<int> bunchIds_;

      /// Weights of waves in bunch averages, indexed by wave id.
      DArray<double> weights_;

      /// Wavenumber values for wave bunches, indexed by bunch id.
      DArray<double> wavenumbers_;

      /// Storage for averaged S(q) values, indexed by bunch id.
      DArray<double> values_;

      /// Average accumulators for S(q), indexed by bunch id.
      DArray<Average> accumulators_;

      /// Constants used in S(q) calculation.
      double a_, b_;

      /// Dimensions of wavevector mesh for real-to-complex DFT.
      IntVec<D> kMeshDimensions_;

      /// Number of wavevectors in wavevector mesh.
      int kSize_;

      /// Number of wavevector bunches (sets of wavevectors of equal norm).
      int nBunch_;

      /// Number of samples per block average output.
      int nSamplePerBlock_;

      /// Has readParam been called?
      bool isInitialized_;

      /// File used for output.
      std::ofstream outputFile_;

   };

   // Explicit instantiation declarations
   extern template class BinaryStructureFactor<1>;
   extern template class BinaryStructureFactor<2>;
   extern template class BinaryStructureFactor<3>;

} // namespace Rpc
} // namespace Pscf
#endif

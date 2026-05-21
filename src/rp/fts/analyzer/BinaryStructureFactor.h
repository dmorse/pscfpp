#ifndef RP_BINARY_STRUCTURE_FACTOR_H
#define RP_BINARY_STRUCTURE_FACTOR_H

/*
* PSCF - Polymer Self-Consistent Field
*
* Copyright 2015 - 2025, The Regents of the University of Minnesota
* Distributed under the terms of the GNU General Public License.
*/

#include <util/accumulators/Average.h>            // member
#include <util/containers/DArray.h>               // member

#include <iostream>

namespace Pscf {
namespace Rp {

   using namespace Util;

   /**
   * Spherically averaged structure factor for a two-monomer system.
   *
   * This class evaluates the structures factors for all wavevectors 
   * within the Fourier space grid, while grouping averages for 
   * wavevectors of equal norm. Sets of wavevectors of equal norm are
   * referred to here as "bunches". 
   *
   * A structure factor for a wavevector k for a system with two types
   * monomer is given by an expectation value
   * \f[
   *     S(k) = 1/(v \chi )^2 <W_(k)W_(-k)>/V - 1/(2 \chi v)
   * \f]
   * where, V is system volume, v is monomer volume, and \f$W_(k)\f$ 
   * is a Fourier transform of the fluctuating exchange field 
   * \f$ W_{-} = (w_{A} - w_{B})/2 \f$. This analyzer outputs the 
   * average value of this quantity for each bunch (or wavenumber
   * value), averaged over waves in a "bunch".
   *
   * \see \ref rp_BinaryStructureFactor_page "Manual Page"
   * \ingroup Rp_Fts_Analyzer_Module
   */
   template <int D, class T>
   class BinaryStructureFactor : public T::Analyzer
   {

   public:

      /**
      * Constructor.
      *
      * \param simulator  parent simulator object
      * \param system  parent system object
      */
      BinaryStructureFactor(typename T::Simulator& simulator, 
                            typename T::System& system);

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
      * Setup immediately before the main loop.
      *
      * Allocates any private data structures that were not allocated
      * previously.
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

      // Alias for base class
      using AnalyzerT = typename T::Analyzer;

      // Inherited protected member functions (selected)
      using AnalyzerT::system;
      using AnalyzerT::simulator;

   protected:

      /**
      * Compute member variables wm_ and wk_.
      */
      void computeW();

      /**
      * Complete calculation of current structure factors.
      */
      void computeS(Array<typename T::Complex> const & wk);

   private:

      /// Exchange field W_(r):  wm = (wa-wb)/2  .
      typename T::RField wm_;

      /// Discrete Fourier transform (DFT) of wm_ . 
      typename T::RFieldDft wk_;

      /// Bunch ids, indexed by wave id (id of bunch containing a wave).
      DArray<int> waveBunchIds_;

      /// Weights of waves in bunch averages, indexed by wave id.
      DArray<double> waveWeights_;

      /// Wavenumber values for wave bunches, indexed by bunch id.
      DArray<double> bunchWavenumbers_;

      /// Storage for averaged S(q) values, indexed by bunch id.
      DArray<double> bunchValues_;

      /// Accumulators for bunch S(q) values, indexed by bunch id.
      DArray<Average> bunchAccumulators_;

      /// Accumulators for wavevector S(q) values, indexed by wave id.
      DArray<Average> waveAccumulators_;

      /// Constants used in S(q) calculation.
      double a_, b_;

      /// Dimensions of k-space mesh for DFT of real field.
      IntVec<D> kMeshDimensions_;

      /// Number of wavevectors in wavevector mesh.
      int nWave_;

      /// Number of wavevector bunches (sets of wavevectors of equal norm).
      int nBunch_;

      /// Should S(q) data be retained for individual waves?
      bool keepWaveData_;

      /// Has readParameters been called?
      bool isInitialized_;

   };

   // Explicit instantiation declarations
   extern template class BinaryStructureFactor<1>;
   extern template class BinaryStructureFactor<2>;
   extern template class BinaryStructureFactor<3>;

} // namespace Rp
} // namespace Pscf
#endif

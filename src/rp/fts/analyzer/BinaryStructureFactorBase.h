#ifndef RP_BINARY_STRUCTURE_FACTOR_BASE_H
#define RP_BINARY_STRUCTURE_FACTOR_BASE_H

/*
* PSCF - Polymer Self-Consistent Field
*
* Copyright 2015 - 2025, The Regents of the University of Minnesota
* Distributed under the terms of the GNU General Public License.
*/

#include <rp/fts/analyzer/Analyzer.h>         // base class
#include <pscf/math/IntVec.h>                 // member
#include <prdc/field/RField.h>                // member
#include <prdc/field/RFieldDft.h>             // member
#include <util/accumulators/Average.h>        // member
#include <util/containers/DArray.h>           // member

#include <pscf/backends/TmplDeclare.h>
#include <iostream>

// Forward declaration
namespace Util {
   template <typename T> class Array;
}
namespace Pscf {
   namespace Prdc {
      template <int D, class T> class RField;
      template <int D, class T> class RFieldDft;
      template <int D, class T> class FFT;
   }
}

namespace Pscf {
namespace Rp {

   using namespace Util;
   using namespace Prdc;

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
   * \f$ W_{-} = (w_{A} - w_{B})/2 \f$. This analyzer outputs the average
   * value of this quantity for each bunch (or wavenumber value), averaged
   * over waves in a "bunch".
   *
   * \see \ref rp_BinaryStructureFactor_page "Manual Page"
   * \ingroup Rp_Fts_Analyzer_Module
   */
   template <int D, class T>
   class BinaryStructureFactorBase : public Analyzer<D,T>
   {

   public:

      /**
      * Constructor.
      *
      * \param simulator  parent simulator object
      * \param system  parent system object
      */
      BinaryStructureFactorBase(Simulator<D,T>& simulator, 
                                System<D,T>& system);

      /**
      * Read parameters from file.
      *
      * For parameter file format see
      * \ref rp_BinaryStructureFactor_page "Manual Page"
      *
      * \param in input parameter stream
      */
      void readParameters(std::istream& in) override;

      // Virual setup and sample functions are implemented by subclasses

      /**
      * Output results to predefined output file.
      *
      * For output files and file format, see
      * \ref rp_BinaryStructureFactor_page "Manual Page"
      */
      void output() override;

   protected:

      /**
      * Allocate member arrays with dimensions that depend only on mesh.
      */
      void allocate();

      /**
      * Allocate and initialize data structures involving wave bunches.
      * 
      * \param kSq  values of square wavenumbers
      * \param implicit  bools indicating existence of implicit inverse
      */
      void findWaveBunches(
                 Array<double> const & kSq,
                 Array<bool> const & implicit);

      /**
      * Compute member variables wm_ and wk_.
      */
      void computeW();

      /**
      * Complete calculation of current structure factors.
      */
      void computeS(Array<typename T::Complex> const & wk);

      /// Discrete Fourier transform (DFT) of wm_ . 
      RFieldDft<D,T> wk_;

      /// Alias for base class
      using AnalyzerT = Analyzer<D,T>;

      // Inherited protected member functions (selected)
      using AnalyzerT::system;
      using AnalyzerT::simulator;

   private:

      /// Exchange field W_(r):  wm = (wa-wb)/2  .
      RField<D,T> wm_;

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

      /// Should S(q) data be outuput for all individual waves?
      bool writeWaveData_;

      /// Has readParameters been called?
      bool isInitialized_;

      // Alias for FFT wrapper class
      using FFTT = FFT<D,T>;

   };

   // Explicit instantiation declarations
   PSCF_TMPL_DECLARE(BinaryStructureFactorBase)

   // Declaration of primary template for subclasses
   template <int D, class T> class BinaryStructureFactor;

} // namespace Rp
} // namespace Pscf
#endif

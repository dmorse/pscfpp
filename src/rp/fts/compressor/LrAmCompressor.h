#ifndef RP_LR_AM_COMPRESSOR_H
#define RP_LR_AM_COMPRESSOR_H

/*
* PSCF - Polymer Self-Consistent Field
*
* Copyright 2015 - 2025, The Regents of the University of Minnesota
* Distributed under the terms of the GNU General Public License.
*/

#include <pscf/iterator/AmIteratorTmpl.h>  // base class template
#include <rp/fts/compressor/Compressor.h>  // base template argument
#include <pscf/math/IntVec.h>              // member
#include <util/containers/DArray.h>        // member

// Forward declarations
namespace Pscf {
   namespace Rp {
      template <int D, class T> class System;
      template <int D, class T> class Simulator;
      template <int D, class T> class Compressor;
      template <int D, class T> class IntraCorrelation;
   }
}

namespace Pscf {
namespace Rp {

   using namespace Util;

   /**
   * Anderson mixing compressor with linear-response correction step.
   *
   * Class LrAmCompressor implements an Anderson mixing algorithm in
   * which the second "correction" step is treated as a quasi-Newton
   * step, while the Jacobian is approximated by the linear response
   * of a hypothetical homogenous liquid. The residual is an r-grid
   * vector in which each element represents the deviation of the
   * sum of volume fractions from unity.
   *
   * Template parameters:
   *
   *    - D : dimension
   *    - T : Types class (e.g., CppTp<D> or CudaTp<D>)
   *
   * \see \ref rp_LrAmCompressor_page "Manual Page"
   * \ingroup Rp_Fts_Compressor_Module
   */
   template <int D, class T>
   class LrAmCompressor
    : public AmIteratorTmpl<Compressor<D,T>, typename T::RDevArray>
   {

   public:

      /**
      * Constructor.
      *
      * \param system  parent System object
      */
      LrAmCompressor(System<D,T>& system);

      /**
      * Destructor.
      */
      ~LrAmCompressor() = default;

      /**
      * Read body of parameter file block and initialize.
      *
      * \param in  input parameter file stream
      */
      void readParameters(std::istream& in) override;

      /**
      * Initialize just before entry to iterative loop.
      *
      * This function is called by the solve function before entering the
      * loop over iterations. It stores the initial values of the fields
      * prior to iteration.
      *
      * \param isContinuation  true iff continuation within a sweep
      */
      void setup(bool isContinuation) override;

      /**
      * Compress to obtain partial saddle point field.
      *
      * \return 0 for convergence, 1 for failure
      */
      int compress() override;

      /**
      * Return compressor time contributions.
      *
      * \param out  output stream
      */
      void outputTimers(std::ostream& out) const override;

      /**
      * Clear all timers and MDE solution counter.
      */
      void clearTimers() override;

   protected:

      // Inherited member function.
      using Compressor<D,T>::system;

   private:

      /// Alias for base class.
      using CompressorT = Compressor<D,T>;

      /// Type of field and residual vectors.
      using VectorT = typename T::RDevArray;

      /// Typename alias for base class.
      using AmTmpl = AmIteratorTmpl< Compressor<D,T>, VectorT >;

      /// Typename alias for FFT class
      using FFTT = typename T::FFT;

      /**
      * Initial values of all w fields.
      */
      DArray< typename T::RField > w0_;

      /**
      * Temporary array of w fields used in update function.
      */
      DArray< typename T::RField > wFieldTmp_;

      /**
      * Residual in real space.
      */
      typename T::RField resid_;

      /**
      * Residual in Fourier space.
      */
      typename T::RFieldDft residK_;

      /**
      * Intramolecular correlation function in Fourier space.
      */
      typename T::RField intraCorrelationK_;

      /**
      * IntraCorrelation object, used to compute intraCorrelationK_.
      */
      IntraCorrelation<D,T> intra_;

      /**
      * Dimensions of wavevector mesh for real-to-complex transform.
      */
      IntVec<D> kMeshDimensions_;

      /**
      * Number of points in k-space wavevector mesh.
      */
      int kSize_;

      /**
      * Has intraCorrelationK_ been calculated?
      */
      bool isIntraCalculated_;

      /**
      * Has required memory been allocated?
      */
      bool isAllocated_;

      // Private AM algorithm operations

      /**
      * Compute and return the number of elements in a field vector.
      *
      * Called during allocation and then stored.
      */
      int nElements() override;

      /**
      * Does the system have an initial guess for the field?
      */
      bool hasInitialGuess() override;

      /**
      * Gets the current field vector from the system.
      *
      * \param curr  current field vector
      */
      void getCurrent(VectorT& curr) override;

      /**
      * Have the system perform a computation using new field.
      *
      * Solves the modified diffusion equations, computes concentrations,
      * and optionally computes stress components.
      */
      void evaluate() override;

      /**
      * Compute the residual vector.
      *
      * \param resid  current residual vector value
      */
      void getResidual(VectorT& resid) override;

      /**
      * Add a correction based on the predicted residual.
      *
      * \param fieldTrial  trial field (in/out)
      * \param resTrial  predicted residual for current trial (in)
      */
      void addCorrection(VectorT& fieldTrial,
                         VectorT const & resTrial) override;

      /**
      * Update the system field with the new trial field.
      *
      * \param newGuess  trial field vector
      */
      void update(VectorT& newGuess) override;

      /**
      * Outputs relevant system details to the iteration log.
      */
      void outputToLog() override;

   };

} // namespace Rp
} // namespace Pscf
#endif

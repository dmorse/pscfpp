#ifndef RP_LR_COMPRESSOR_H
#define RP_LR_COMPRESSOR_H

/*
* PSCF - Polymer Self-Consistent Field
*
* Copyright 2015 - 2025, The Regents of the University of Minnesota
* Distributed under the terms of the GNU General Public License.
*/

#include <rp/fts/compressor/Compressor.h>    // base class
#include <pscf/math/IntVec.h>                // member
#include <util/containers/DArray.h>          // member
#include <util/misc/Timer.h>                 // member

// Forward declarations
namespace Pscf {
   namespace Rp {
      template <int D, class T> class System;
      template <int D, class T> class Simulator;
      template <int D, class T> class IntraCorrelation;
   }
}

namespace Pscf {
namespace Rp {

   using namespace Util;

   /**
   * Linear response compressor.
   *
   * This class implements a compressor that is an approximate Newton's
   * method in which the Jacobian is approximated by the analytically
   * calculated linear response of a homogeneous liquid with the same
   * composition as the system of interest.

   * Specializations of this class template are used as base classes for 
   * two closely analogous class templates, also both named LrCompressor,
   * that are defined in Rpc and Rpg namespaces for use in the pscf_rpc
   * and pscf_rpg programs, respectively.
   *
   * Template parameters:
   *
   *   - D : dimension of space (D=1, 2, or 3)
   *   - T : Types class (CppTp<D> or CudaTp<D>)
   *
   * \see \ref rp_LrCompressor_page
   * \ingroup Rp_Fts_Compressor_Module
   */
   template <int D, class T>
   class LrCompressor : public Compressor<D,T>
   {

   public:

      /**
      * Constructor.
      *
      * \param system  parent System object
      */
      LrCompressor(System<D,T>& system);

      /**
      * Destructor.
      */
      ~LrCompressor() = default;

      /**
      * Read all parameters and initialize.
      *
      * \param in  input parameter file stream
      */
      void readParameters(std::istream& in);

      /**
      * Initialize just before entry to iterative loop.
      *
      * This function is called by the solve function before entering the
      * loop over iterations. Store the current values of the fields at the
      * beginning of iteration
      */
      void setup();

      /**
      * Iterate pressure field to obtain partial saddle point.
      *
      * \return 0 for convergence, 1 for failure
      */
      int compress();

      /**
      * Return compressor times contributions.
      *
      * \param out  output stream
      */
      void outputTimers(std::ostream& out) const;

      /**
      * Clear all timers.
      */
      void clearTimers();

   protected:

      // Inherited protected members
      using CompressorT = Compressor<D,T>;
      using CompressorT::mdeCounter_;
      using CompressorT::system;

   private:

      using RFieldT = typename T::RField;
      using RFieldDftT = typename T::RFieldDft;
      using FFTT = typename T::FFT;

      // IntraCorrelation object
      IntraCorrelation<D,T> intra_;

      // Template w Field used in update function
      DArray< RFieldT > wFieldTmp_;

      // Residual in real space
      RFieldT resid_;

      // Residual in Fourier space
      RFieldDftT residK_;

      // Intramolecular correlation in Fourier space
      RFieldT intraCorrelationK_;

      // Dimensions of wavevector mesh in real-to-complex transform
      IntVec<D> kMeshDimensions_;

      // Timers for analyzing performance
      Timer timerTotal_;
      Timer timerMDE_;

      // Type of error criterion used to test convergence
      std::string errorType_;

      // Error tolerance
      double epsilon_;

      // Current iteration counter
      int itr_;

      // Maximum number of iterations
      int maxItr_;

      // Total iteration counter
      int totalItr_;

      // Verbosity level
      int verbose_;

      // Has required memory been allocated?
      bool isAllocated_;

      // Has the IntraCorrelation been calculated?
      bool isIntraCalculated_;

      // Private member functions

      /**
      * Compute the residual vector.
      */
      void computeResidual();

      /**
      * Update system w fields.
      */
      void updateWFields();

      /**
      * Compute and return error used to test for convergence.
      *
      * \param verbose  verbosity level of output report
      * \return error  measure used to test for convergence.
      */
      double computeError(int verbose);

      /**
      * Outputs relevant system details to the iteration log.
      */
      void outputToLog();

   };

} // namespace Rp
} // namespace Pscf
#endif

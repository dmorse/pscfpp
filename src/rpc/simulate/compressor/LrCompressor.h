#ifndef RPC_LR_COMPRESSOR_H
#define RPC_LR_COMPRESSOR_H

/*
* PSCF - Polymer Self-Consistent Field Theory
*
* Copyright 2016 - 2022, The Regents of the University of Minnesota
* Distributed under the terms of the GNU General Public License.
*/

#include "Compressor.h"
#include <prdc/cpu/RField.h>
#include <prdc/cpu/RFieldDft.h>
#include <util/containers/DArray.h>
#include <util/containers/DMatrix.h>
#include <util/misc/Timer.h>

namespace Pscf {
namespace Rpc
{

   template <int D> class System;

   using namespace Util;
   using namespace Pscf::Prdc::Cpu;

   /**
   * Linear response compressor.
   *
   * This class implements a compressor that is an approximate Newton's
   * method in which the Jacobian is approximated by the analytically
   * calculated linear response of a homogeneous liquid with the same
   * composition as the system of interest.
   *
   * \ingroup Rpc_Compressor_Module
   */
   template <int D>
   class LrCompressor : public Compressor<D>
   {

   public:

      /**
      * Constructor.
      *
      * \param system System object associated with this compressor.
      */
      LrCompressor(System<D>& system);

      /**
      * Destructor.
      */
      ~LrCompressor();

      /**
      * Read all parameters and initialize.
      *
      * \param in input stream for parameter file
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

      double subspacePercent(){return 0;};

      double correctionPercent(){return 0;};

      /**
      * Return compressor times contributions.
      */
      void outputTimers(std::ostream& out);

      void clearTimers();

   protected:

      // Inherited protected members
      using Compressor<D>::mdeCounter_;
      using ParamComposite::read;
      using ParamComposite::readOptional;
      using ParamComposite::setClassName;

   private:

      // Error tolerance.
      double epsilon_;

      // Current iteration counter.
      int itr_;

      // Maximum number of iterations.
      int maxItr_;

      // Timers for analyzing performance.
      Timer timerTotal_;
      Timer timerMDE_;

      // Type of error criterion used to test convergence
      std::string errorType_;

      // Verbosity level.
      int verbose_;

      // Has required memory been allocated.
      bool isAllocated_;

      /**
      * Dimensions of wavevector mesh in real-to-complex transform
      */
      IntVec<D> kMeshDimensions_;

      /**
      * IntraCorrelation.
      */
      RField<D> intraCorrelation_;

      /**
      * Residual in real space used for linear response anderson mixing.
      */
      RField<D> resid_;

      /**
      * Residual in Fourier space used for linear response anderson mixing.
      */
      RFieldDft<D> residK_;

      /**
      * Current values of the fields
      */
      DArray< RField<D> > w0_;

      /**
      * Template w Field used in update function
      */
      DArray< RField<D> > wFieldTmp_;

      /**
      * Compute the residual vector.
      *
      * \param resid current residual vector value
      */
      void getResidual();

      /**
      * Use homogenous linear respose analytic approximation to update wFields
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
      * Find the maximum magnitude element of a residual field.
      */
      double maxAbs(RField<D> const & a);

      /**
      * Compute the inner product of two fields.
      *
      * \param a first field
      * \param b second field
      */
      double dotProduct(RField<D> const & a, RField<D> const & b);

      /**
      * Compute L2 norm
      */
      double norm(RField<D> const & a);

      /**
      * Outputs relevant system details to the iteration log.
      */
      void outputToLog();
      
      // Private inherited members
      using Compressor<D>::system;

   };

   #ifndef RPC_LR_COMPRESSOR_TPP
   // Suppress implicit instantiation
   extern template class LrCompressor<1>;
   extern template class LrCompressor<2>;
   extern template class LrCompressor<3>;
   #endif

} // namespace Rpc
} // namespace Pscf
#endif

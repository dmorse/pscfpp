#ifndef PSPC_LR_AM_COMPRESSOR_H
#define PSPC_LR_AM_COMPRESSOR_H

/*
* PSCF - Polymer Self-Consistent Field Theory
*
* Copyright 2016 - 2022, The Regents of the University of Minnesota
* Distributed under the terms of the GNU General Public License.
*/

#include "Compressor.h"
#include <prdc/cpu/RField.h>
#include <prdc/cpu/RFieldDft.h>
#include <pscf/iterator/AmIteratorTmpl.h>                 

namespace Pscf {
namespace Pspc
{

   template <int D> class System;

   using namespace Util;
   using namespace Pscf::Prdc::Cpu;

   /**
   * Anderson Mixing compressor with linear-response preconditioning.
   *
   * Class LrAmCompressor implements an Anderson mixing algorithm in
   * which the residual is defined using a preconditioning scheme that
   * would yield a Jacobian of unity if applied to a homogeneous system.
   * The residual in the unpreconditioned form of Anderson mixing is a
   * vector in which each that represents a deviations in the sum of 
   * volume fractions from unity. In this preconditioned algorithm, each
   * Fourier component of this deviation is multiplied by the inverse of
   * Fourier representation of the linear response of total concentration 
   * to changes in pressure in a homogeneous system of the same chemical
   * composition as the system of interest.
   *
   * \ingroup Pspc_Compressor_Module
   */
   template <int D>
   class LrAmCompressor 
         : public AmIteratorTmpl<Compressor<D>, DArray<double> >
   {

   public:

      /**
      * Constructor.
      * 
      * \param system System object associated with this compressor.
      */
      LrAmCompressor(System<D>& system);

      /**
      * Destructor.
      */
      ~LrAmCompressor();

      /**
      * Read all parameters and initialize.
      *
      * \param in input filestream
      */
      void readParameters(std::istream& in);

      /**
      * Initialize just before entry to iterative loop.
      *
      * This function is called by the solve function before entering the
      * loop over iterations. Store the current values of the fields at the 
      * beginning of iteration
      *
      * \param isContinuation true iff continuation within a sweep
      */ 
      void setup(bool isContinuation);      
      
      
      /**
      * Compress to obtain partial saddle point w+
      *
      * \return 0 for convergence, 1 for failure
      */
      int compress();    
      
      /**
      * Write a report of time contributions used by this algorithm.
      * 
      * \param out  output stream to which to write
      */
      void outputTimers(std::ostream& out);

      /**
      * Reset / clear all timers.
      */
      void clearTimers();

      /**
      * Compute and return the scalar error.
      *
      * \param verbose  verbosity level (higher is more verbose)
      */
      double computeError(int verbose);
      
      
      // Inherited public member functions

      using AmIteratorTmpl<Compressor<D>, DArray<double> >::setClassName;
      
   protected:
  
      // Inherited protected members 
      using ParamComposite::readOptional;
      using Compressor<D>::mdeCounter_;

   private:

      /**
      * Type of error criterion used to test convergence 
      */ 
      std::string errorType_;
      
      /**
      * How many times MDE has been solved for each mc move.
      */
      int itr_;
  
      /**
      * Current values of the fields
      */
      DArray< RField<D> > w0_;  
      
      DArray<double> error_;
      
      /**
      * Residual in real space used for linear response anderson mixing.
      */
      RField<D> resid_;
      
      /**
      * Residual in Fourier space used for linear response anderson mixing.
      */
      RFieldDft<D> residK_;
     
      /**
      * IntraCorrelation.
      */
      RField<D> intraCorrelation_;
      
      /**
      * Dimensions of wavevector mesh in real-to-complex transform
      */ 
      IntVec<D> kMeshDimensions_;
      
      /**
      * Has the variable been allocated?
      */
      bool isAllocated_;
      
      /**
      * Template w Field used in update function.
      */
      
      DArray< RField<D> > wFieldTmp_;
      
      /**
      * New Basis variable used in updateBasis function. 
      */
      DArray<double> newBasis_;

      /**
      * Assign one field to another.
      * 
      * \param a the field to be set (lhs of assignment)
      * \param b the field for it to be set to (rhs of assigment)
      */
      void setEqual(DArray<double>& a, DArray<double> const & b);

      /**
      * Compute the inner product of two vectors
      */
      double dotProduct(DArray<double> const & a, DArray<double> const & b);

      /**
      * Find the maximum magnitude element of a residual vector.
      */
      double maxAbs(DArray<double> const & hist);

      /**
      * Update the basis for residual or field vectors.
      * 
      * \param basis RingBuffer of residual or field basis vectors
      * \param hists RingBuffer of past residual or field vectors
      */
      void updateBasis(RingBuffer<DArray<double> > & basis, 
                       RingBuffer<DArray<double> > const & hists);

      /**
      * Add linear combination of basis vectors to trial field.
      * 
      * \param trial trial vector (input-output)
      * \param basis RingBuffer of basis vectors
      * \param coeffs array of coefficients of basis vectors
      * \param nHist number of histories stored at this iteration
      */
      void addHistories(DArray<double>& trial, 
                        RingBuffer<DArray<double> > const & basis, 
                        DArray<double> coeffs, 
                        int nHist);

      /**
      * Add predicted error to field trial.
      * 
      * \param fieldTrial trial field (in-out)
      * \param resTrial predicted error for current trial
      * \param lambda Anderson-Mixing mixing 
      */
      void addPredictedError(DArray<double>& fieldTrial, 
                             DArray<double> const & resTrial, 
                             double lambda);

      /**
      * Does the system has an initial guess for the field?
      */
      bool hasInitialGuess();
     
      /** 
      * Compute and returns the number of elements in field vector.
      *
      * Called during allocation and then stored.
      */
      int nElements();

      /**
      * Gets the current field vector from the system.
      * 
      * \param curr current field vector
      */ 
      void getCurrent(DArray<double>& curr);

      /**
      * Have the system perform a computation using new field.
      *
      * Solves the modified diffusion equations, computes concentrations,
      * and optionally computes stress components.
      */
      void evaluate();

      /**
      * Compute the residual vector.
      *
      * \param resid current residual vector value
      */
      void getResidual(DArray<double>& resid);

      /**
      * Updates the system field with the new trial field.
      *
      * \param newGuess trial field vector
      */
      void update(DArray<double>& newGuess);

      /**
      * Outputs relevant system details to the iteration log.
      */
      void outputToLog();
      
      /**
      * Compute Debye function
      */
      double computeDebye(double x);
      
      /**
      * Compute intramolecular correlation at specific sqSquare
      */
      double computeIntraCorrelation(double qSquare);
      
      /**
      * Compute intramolecular correlation  
      */
      void computeIntraCorrelation();
      
      using Compressor<D>::system;

   };

} // namespace Pspc
} // namespace Pscf
#endif

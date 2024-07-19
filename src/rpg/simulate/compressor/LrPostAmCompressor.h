#ifndef RPG_LR_POST_AM_COMPRESSOR_H
#define RPG_LR_POST_AM_COMPRESSOR_H

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
#include <rpg/simulate/compressor/intra/IntraCorrelation.h> 

namespace Pscf {
namespace Rpg
{

   template <int D> class System;
   template <int D> class IntraCorrelation;

   using namespace Util;
   using namespace Pscf::Prdc;
   using namespace Pscf::Prdc::Cuda;

   /**
   * Anderson Mixing compressor with linear-response mixing step.
   *
   * Class LrPostAmCompressor implements an Anderson mixing algorithm 
   * which modifies the second mixing step, estimating Jacobian by linear 
   * response of homogenous liquid instead of unity. The residual is a 
   * vector in which each that represents a deviations 
   * in the sum of volume fractions from unity.
   *
   * \ingroup Rpg_Simulate_Compressor_Module
   */
   template <int D>
   class LrPostAmCompressor 
         : public AmIteratorTmpl<Compressor<D>, Field<cudaReal> >
   {

   public:

      /**
      * Constructor.
      * 
      * \param system System object associated with this compressor.
      */
      LrPostAmCompressor(System<D>& system);

      /**
      * Destructor.
      */
      ~LrPostAmCompressor();

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
      * Return compressor times contributions.
      */
      void outputTimers(std::ostream& out);

      /**
      * Clear all timers (reset accumulated time to zero).
      */
      void clearTimers();
      
      // Inherited public member functions
      using AmIteratorTmpl<Compressor<D>, Field<cudaReal> >::setClassName;
      
   protected:
  
      // Inherited protected members 
      using ParamComposite::readOptional;
      using Compressor<D>::mdeCounter_;

   private:
   
      /**
      * How many times MDE has been solved for each mc move 
      */
      int itr_;
      
      /**
      * Current values of the fields
      */
      DArray< RField<D> > w0_;  

      /**
      * Has the variable been allocated?
      */
      bool isAllocated_;
      
      /**
      * Template w Field used in update function
      */
      DArray< RField<D> > wFieldTmp_;
      
      /**
      * New Basis variable used in updateBasis function 
      */
      Field<cudaReal> newBasis_;
      
      /**
      * Residual in real space used for linear response anderson mixing.
      */
      RField<D> resid_;
      
      /**
      * Residual in Fourier space used for linear response anderson mixing.
      */
      RFieldDft<D> residK_;
     
      /**
      * IntraCorrelation in fourier space calculated by IntraCorrlation class
      */
      RField<D> intraCorrelationK_;
      
      /**
      * Dimensions of wavevector mesh in real-to-complex transform
      */ 
      IntVec<D> kMeshDimensions_;
      
      /**
      * Number of points in k-space grid
      */
      int kSize_;

      /**
      * Assign one field to another.
      * 
      * \param a the field to be set (lhs of assignment)
      * \param b the field for it to be set to (rhs of assigment)
      */
      void setEqual(Field<cudaReal>& a, Field<cudaReal> const & b);

      /**
      * Compute the inner product of two vectors
      */
      double dotProduct(Field<cudaReal> const & a, Field<cudaReal> const & b);

      /**
      * Find the maximum magnitude element of a residual vector.
      */
      double maxAbs(Field<cudaReal> const & hist);

      /**
      * Update the basis for residual or field vectors.
      * 
      * \param basis RingBuffer of residual or field basis vectors
      * \param hists RingBuffer of past residual or field vectors
      */
      void updateBasis(RingBuffer<Field<cudaReal> > & basis, 
                       RingBuffer<Field<cudaReal> > const & hists);

      /**
      * Add linear combination of basis vectors to trial field.
      * 
      * \param trial trial vector (input-output)
      * \param basis RingBuffer of basis vectors
      * \param coeffs array of coefficients of basis vectors
      * \param nHist number of histories stored at this iteration
      */
      void addHistories(Field<cudaReal>& trial, 
                        RingBuffer<Field<cudaReal> > const & basis, 
                        DArray<double> coeffs, 
                        int nHist);

      /**
      * Add predicted error to field trial.
      * 
      * \param fieldTrial trial field (in-out)
      * \param resTrial predicted error for current trial
      * \param lambda Anderson-Mixing mixing 
      */
      void addPredictedError(Field<cudaReal>& fieldTrial, 
                             Field<cudaReal> const & resTrial, 
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
      void getCurrent(Field<cudaReal>& curr);

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
      void getResidual(Field<cudaReal>& resid);

      /**
      * Updates the system field with the new trial field.
      *
      * \param newGuess trial field vector
      */
      void update(Field<cudaReal>& newGuess);

      /**
      * Outputs relevant system details to the iteration log.
      */
      void outputToLog();
      
      /**
      * Set mixing parameter lambda
      */
      double setLambda();
      
      /**
      * IntraCorrelation (homopolymer) object
      */
      IntraCorrelation<D> intraCorrelation_;
    
      // Inherited private members 
      using Compressor<D>::system;

   };
   
   #ifndef RPG_LR_POST_AM_COMPRESSOR_TPP
   // Suppress implicit instantiation
   extern template class LrPostAmCompressor<1>;
   extern template class LrPostAmCompressor<2>;
   extern template class LrPostAmCompressor<3>;
   #endif

} // namespace Rpg
} // namespace Pscf
#endif

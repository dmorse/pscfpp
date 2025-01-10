#ifndef RPG_AM_COMPRESSOR_H
#define RPG_AM_COMPRESSOR_H

/*
* PSCF - Polymer Self-Consistent Field Theory
*
* Copyright 2016 - 2022, The Regents of the University of Minnesota
* Distributed under the terms of the GNU General Public License.
*/

#include "Compressor.h"
#include <prdc/cuda/RField.h>         
#include <pscf/cuda/DeviceArray.h>
#include <pscf/iterator/AmIteratorTmpl.h>     

#include <util/containers/DArray.h>                 

namespace Pscf {
namespace Rpg
{

   template <int D> class System;

   using namespace Util;
   using namespace Pscf::Prdc;
   using namespace Pscf::Prdc::Cuda;

   /**
   * Rpg implementation of the Anderson Mixing compressor.
   *
   * \ingroup Rpg_Fts_Compressor_Module
   */
   template <int D>
   class AmCompressor : public AmIteratorTmpl<Compressor<D>, 
                               DeviceArray<cudaReal> >
   {

   public:

      /**
      * Constructor.
      * 
      * \param system System object associated with this compressor.
      */
      AmCompressor(System<D>& system);

      /**
      * Destructor.
      */
      ~AmCompressor();

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
      * Compute mixing parameter lambda
      */
      double computeLambda(double r);  
      
      /**
      * Return how many times MDE has been solved.
      */
      int mdeCounter(); 
      
      /**
      * Return compressor times contributions.
      */
      void outputTimers(std::ostream& out);
      
      /**
      * Clear all timers (reset accumulated time to zero).
      */
      void clearTimers();
      
      // Inherited public member functions
      using AmIteratorTmpl<Compressor<D>, 
                           DeviceArray<cudaReal> >::setClassName;

   protected:
  
      // Inherited protected members 
      using ParamComposite::readOptional;
      using Compressor<D>::system;
      using Compressor<D>::mdeCounter_;

   private:
      /**
      * Count how many times MDE has been solved.
      */
      int counter_;
      
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
      * Assign one field to another.
      * 
      * \param a the field to be set (lhs of assignment)
      * \param b the field for it to be set to (rhs of assigment)
      */
      void setEqual(DeviceArray<cudaReal>& a, 
                    DeviceArray<cudaReal> const & b);

      /**
      * Compute the inner product of two vectors
      */
      double dotProduct(DeviceArray<cudaReal> const & a, 
                        DeviceArray<cudaReal> const & b);

      /**
      * Find the maximum magnitude element of a residual vector.
      */
      double maxAbs(DeviceArray<cudaReal> const & hist);

      /**
      * Update the basis for residual or field vectors.
      * 
      * \param basis RingBuffer of residual or field basis vectors
      * \param hists RingBuffer of past residual or field vectors
      */
      void updateBasis(RingBuffer<DeviceArray<cudaReal> > & basis, 
                       RingBuffer<DeviceArray<cudaReal> > const & hists);

      /**
      * Add linear combination of basis vectors to trial field.
      * 
      * \param trial trial vector (input-output)
      * \param basis RingBuffer of basis vectors
      * \param coeffs array of coefficients of basis vectors
      * \param nHist number of histories stored at this iteration
      */
      void addHistories(DeviceArray<cudaReal>& trial, 
                        RingBuffer<DeviceArray<cudaReal> > const & basis, 
                        DArray<double> coeffs, 
                        int nHist);

      /**
      * Add predicted error to field trial.
      * 
      * \param fieldTrial trial field (in-out)
      * \param resTrial predicted error for current trial
      * \param lambda Anderson-Mixing mixing 
      */
      void addPredictedError(DeviceArray<cudaReal>& fieldTrial, 
                             DeviceArray<cudaReal> const & resTrial, 
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
      void getCurrent(DeviceArray<cudaReal>& curr);

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
      void getResidual(DeviceArray<cudaReal>& resid);

      /**
      * Updates the system field with the new trial field.
      *
      * \param newGuess trial field vector
      */
      void update(DeviceArray<cudaReal>& newGuess);

      /**
      * Outputs relevant system details to the iteration log.
      */
      void outputToLog();
      

   };
   
} // namespace Rpg
} // namespace Pscf
#endif

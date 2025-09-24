#ifndef RPG_LR_AM_COMPRESSOR_H
#define RPG_LR_AM_COMPRESSOR_H

/*
* PSCF - Polymer Self-Consistent Field
*
* Copyright 2015 - 2025, The Regents of the University of Minnesota
* Distributed under the terms of the GNU General Public License.
*/

#include "Compressor.h"                          // base class argument
#include <pscf/cuda/DeviceArray.h>               // base class argument
#include <pscf/iterator/AmIteratorTmpl.h>        // base class template

#include <rpg/fts/compressor/IntraCorrelation.h> // member
#include <prdc/cpu/RField.h>                     // member
#include <prdc/cpu/RFieldDft.h>                  // member
#include <util/containers/DArray.h>              // member

#include <iostream>                 

namespace Pscf {
namespace Rpg {

   // Forward declaration
   template <int D> class System;

   using namespace Util;
   using namespace Pscf::Prdc;
   using namespace Pscf::Prdc::Cuda;

   /**
   * Anderson Mixing compressor with linear-response mixing step.
   *
   * Class LrAmCompressor implements an Anderson mixing algorithm 
   * which modifies the second mixing step, estimating Jacobian by linear 
   * response of homogenous liquid instead of unity. The residual is a 
   * vector in which each that represents a deviations 
   * in the sum of volume fractions from unity.
   *
   * \ingroup Rpg_Fts_Compressor_Module
   */
   template <int D>
   class LrAmCompressor 
         : public AmIteratorTmpl<Compressor<D>, DeviceArray<cudaReal> >
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
      * Return compressor times contributions.
      */
      void outputTimers(std::ostream& out) const;

      /**
      * Clear all timers (reset accumulated time to zero).
      */
      void clearTimers();
      
      /// Typename alias for base class.
      using Base = AmIteratorTmpl< Compressor<D>, DeviceArray<cudaReal> >;

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
      * Template w Field used in update function
      */
      DArray< RField<D> > wFieldTmp_;
      
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
      * Has the IntraCorrelation been calculated?
      */
      bool isIntraCalculated_;
      
      /**
      * IntraCorrelation object
      */
      IntraCorrelation<D> intra_;
      
      /**
      * Has the variable been allocated?
      */
      bool isAllocated_;
   
      // AM algorithm operations
    
      /**
      * Add predicted error to field trial.
      * 
      * \param fieldTrial trial field (in-out)
      * \param resTrial predicted error for current trial
      * \param lambda Anderson-Mixing mixing 
      */
      void addPredictedError(DeviceArray<cudaReal>& fieldTrial, 
                             DeviceArray<cudaReal> const & resTrial, 
                             double lambda) override;

      #if 0
      /**
      * Compute mixing parameter lambda
      *
      * \param r  ratio in ramp 
      */
      double computeLambda(double r);
      #endif
      
      // Private virtual that interact with the parent System
      
      /** 
      * Compute and returns the number of elements in field vector.
      *
      * Called during allocation and then stored.
      */
      int nElements() override;

      /**
      * Does the system has an initial guess for the field?
      */
      bool hasInitialGuess() override;
     
      /**
      * Gets the current field vector from the system.
      * 
      * \param curr current field vector
      */ 
      void getCurrent(DeviceArray<cudaReal>& curr) override;

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
      * \param resid current residual vector value
      */
      void getResidual(DeviceArray<cudaReal>& resid) override;

      /**
      * Updates the system field with the new trial field.
      *
      * \param newGuess trial field vector
      */
      void update(DeviceArray<cudaReal>& newGuess) override;

      /**
      * Outputs relevant system details to the iteration log.
      */
      void outputToLog() override;
 
      // Private virtual functions for vector math
      
      /**
      * Assign one field to another.
      * 
      * \param a the field to be set (lhs of assignment)
      * \param b the field for it to be set to (rhs of assigment)
      */
      void setEqual(DeviceArray<cudaReal>& a, 
                    DeviceArray<cudaReal> const & b) override;

      /**
      * Compute the inner product of two vectors
      *
      * \param a  first input vector
      * \param a  second input vector
      */
      double dotProduct(DeviceArray<cudaReal> const & a, 
                        DeviceArray<cudaReal> const & b) override;

      /**
      * Find the maximum magnitude element of a residual vector.
      */
      double maxAbs(DeviceArray<cudaReal> const & hist) override;

      /**
      * Compute the difference a = b - c for vectors a, b and c.
      *
      * \param a result vector (LHS)
      * \param b first vector (RHS)
      * \param c second vector (RHS)
      */
      void subVV(DeviceArray<cudaReal>& a, 
                 DeviceArray<cudaReal> const & b, 
		 DeviceArray<cudaReal> const & c) override;

      /**
      * Compute a += c*b for vectors a and b and scalar c.
      *
      * \param a result vector (LHS)
      * \param b input vector (RHS)
      * \param c scalar coefficient (RHS)
      */
      void addEqVc(DeviceArray<cudaReal>& a, 
		   DeviceArray<cudaReal> const & b, 
                   double c) override;

      // Inherited private member
      using Compressor<D>::system;

   };
   
   // Explicit instantiation declarations
   extern template class LrAmCompressor<1>;
   extern template class LrAmCompressor<2>;
   extern template class LrAmCompressor<3>;

} // namespace Rpg
} // namespace Pscf
#endif

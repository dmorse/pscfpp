#ifndef PSPC_AM_COMPRESSOR_H
#define PSPC_AM_COMPRESSOR_H

/*
* PSCF - Polymer Self-Consistent Field Theory
*
* Copyright 2016 - 2022, The Regents of the University of Minnesota
* Distributed under the terms of the GNU General Public License.
*/

#include "Compressor.h"
#include <pscf/iterator/AmIteratorTmpl.h>                 

namespace Pscf {
namespace Pspc
{

   template <int D>
   class System;

   using namespace Util;

   /**
   * Pspc implementation of the Anderson Mixing compressor.
   *
   * \ingroup Pspc_Compressor_Module
   */
   template <int D>
   class AmCompressor 
         : public AmIteratorTmpl<Compressor<D>, DArray<double> >
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
      * compress to obtain partial saddle point w+
      *
      * \return 0 for convergence, 1 for failure
      */
      int compress();    
      
      
      // Inherited public member functions
      using AmIteratorTmpl<Compressor<D>, DArray<double> >::setClassName;

   protected:
  
      // Inherited protected members 
      using ParamComposite::readOptional;
      using Compressor<D>::system;

   private:
       /**
       * Current values of the fields
       */
      DArray< DArray<double> > w0;  

      /**
      * Find L2 norm of a residual vector.
      */
      double findNorm(DArray<double> const & hist);

      /**
      * Find the maximum magnitude element of a residual vector.
      */
      double findMaxAbs(DArray<double> const & hist);

      /**
      * Update the basis for residual or field vectors.
      * 
      * \param basis RingBuffer of residual or field basis vectors
      * \param hists RingBuffer of past residual or field vectors
      */
      void updateBasis(RingBuffer<DArray<double> > & basis, 
                       RingBuffer<DArray<double> > const & hists);

      /**
      * Compute the dot product for one element of the U matrix.
      * 
      * \param resBasis RingBuffer of residual basis vectors.
      * \param m row of the U matrix
      * \param n column of the U matrix
      */
      double computeUDotProd(RingBuffer<DArray<double> > const & resBasis, 
                             int m, int n);

      /**
      * Compute the dot product for one element of the v vector.
      * 
      * \param resCurrent current residual vector 
      * \param resBasis RingBuffer of residual basis vectors
      * \param m row index of the v vector
      */
      double computeVDotProd(DArray<double> const & resCurrent, 
                             RingBuffer<DArray<double> > const & resBasis, 
                             int m);

      /**
      * Update the U matrix.
      * 
      * \param U U matrix
      * \param resBasis RingBuffer of residual basis vectors.
      * \param nHist number of past states
      */
      void updateU(DMatrix<double> & U, 
                   RingBuffer<DArray<double> > const & resBasis, 
                   int nHist);

      /**
      * Update the v vector.
      * 
      * \param v v vector
      * \param resCurrent current residual vector 
      * \param resBasis RingBuffer of residual basis vectors.
      * \param nHist number of past states 
      */
      void updateV(DArray<double> & v, 
                   DArray<double> const & resCurrent, 
                   RingBuffer<DArray<double> > const & resBasis, 
                   int nHist);

      /**
      * Assign one field to another.
      * 
      * \param a the field to be set (lhs of assignment)
      * \param b the field for it to be set to (rhs of assigment)
      */
      void setEqual(DArray<double>& a, DArray<double> const & b);

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

   };
   #ifndef PSPC_AM_COMPRESSOR_TPP
   // Suppress implicit instantiation
   extern template class AmCompressor<1>;
   extern template class AmCompressor<2>;
   extern template class AmCompressor<3>;
   #endif
} // namespace Pspc
} // namespace Pscf
#endif

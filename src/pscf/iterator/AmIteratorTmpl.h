#ifndef PSCF_AM_ITERATOR_TMPL_H
#define PSCF_AM_ITERATOR_TMPL_H

/*
* PSCF - Polymer Self-Consistent Field Theory
*
* Copyright 2016 - 2019, The Regents of the University of Minnesota
* Distributed under the terms of the GNU General Public License.
*/

#include <pscf/math/LuSolver.h>

#include <util/containers/DArray.h>
#include <util/containers/FArray.h>
#include <util/containers/FSArray.h>
#include <util/containers/DMatrix.h>
#include <util/containers/RingBuffer.h>

namespace Pscf {

   using namespace Util;
   
   /**
   * Anderson mixing iterator for the pseudo spectral method.
   *
   * \ingroup Pscf_Iterator_Module
   */
   template <typename Iterator, typename T>
   class AmIteratorTmpl : virtual public Iterator
   {
   public:

      /**
      * Constructor
      *
      */
      AmIteratorTmpl();

      /**
      * Destructor
      */
      ~AmIteratorTmpl();

      /**
      * Read all parameters and initialize.
      *
      * \param in input filestream
      */
      void readParameters(std::istream& in);

      /**
      * Setup and allocate required memory.
      */
      void setup();

      /**
      * Iterate to a solution
      */
      int solve();

   private:

      /// Error tolerance
      double epsilon_;

      /// Type of error checked for convergence.
      /// Either maxResid or normResid.
      std::string errorType_;

      /// Free parameter for minimization
      double lambda_;

      /// Number of previous steps to use to compute next state. [0,maxHist_]
      int nHist_;

      /// Number of histories to retain.
      int maxHist_;

      /// Maximum number of iterations to attempt.
      int maxItr_;

      /// Number of elements in the field.
      int nElem_; 

      /// Field histories
      RingBuffer< T > fieldHists_;

      /// Basis vectors of field histories.
      RingBuffer< T > fieldBasis_;

      /// Residual histories
      RingBuffer< T > resHists_;

      /// Basis vectors of residual histories
      RingBuffer< T > resBasis_;

      /// Matrix containing the dot products of residual history
      /// basis vectors in resBasis_.
      DMatrix<double> U_;

      /// Coefficients for mixing previous states.
      DArray<double> coeffs_;

      /// Vector of dot products of current residual with residual history
      /// basis vectors in resBasis_.
      DArray<double> v_;

      /// New trial field (big W in Arora et al. 2017)
      T fieldTrial_;

      /// Predicted field residual for trial state (big D)
      T resTrial_;

      /// Workspace for calculations.
      T temp_;

      /**
      * Compute the deviation of wFields from a mean field solution
      */
      void computeResidual();

      /**
      * Check if solution is converge within specified tolerance.
      *
      * \return true for error < epsilon and false for error >= epsilon
      */
      bool isConverged();

      /**
      * Determine the coefficients that would minimize U_
      */
      void findResidCoeff();

      /**
      * Rebuild wFields for the next iteration from minimized coefficients
      */
      void updateGuess();

      /**
      * Clean up after a call to solve(), enabling future calls to solve.
      */
      void cleanUp();

      // ---- Methods for doing AM iterator math ---- //

      /**
      * Find norm of a residual vector.
      */
      virtual double findNorm(T const & hist) = 0;

      /**
      * Find the maximum magnitude element of a residual vector.
      */
      virtual double findMaxAbs(T const & hist) = 0;

      /**
      * Update the series of residual vectors.
      * 
      * \param basis RingBuffer object storing the list of residual or field basis vectors.
      * \param hists RingBuffer object storing the histories of residual or field vectors.
      */
      virtual void updateBasis(RingBuffer<T> & basis, RingBuffer<T> const & hists) = 0;

      /**
      * Compute the dot product for an element of the U matrix.
      * 
      * \param resBasis RingBuffer object storing the list of residual basis vectors.
      * \param m row of the U matrix
      * \param n column of the U matrix
      */
      virtual double computeUDotProd(RingBuffer<T> const & resBasis, int m, int n) = 0;
      
      /**
      * Compute the dot product for an element of the v vector.
      * 
      * \param resCurrent the residual vector calculated at the present iteration step
      * \param resBasis RingBuffer object storing the list of residual basis vectors.
      * \param m row of the v vector
      */
      virtual double computeVDotProd(T const & resCurrent, RingBuffer<T> const & resBasis, int m) = 0;
      
      /**
      * Compute the series of necessary dot products and update the U matrix.
      * 
      * \param U U matrix
      * \param resBasis RingBuffer object storing the list of residual basis vectors.
      * \param nHist number of histories stored at this iteration
      */      
      virtual void updateU(DMatrix<double> & U, RingBuffer<T> const & resBasis, int nHist) = 0;

      /**
      * Compute the series of necessary dot products and update the v vector.
      * 
      * \param v v vector
      * \param resCurrent the residual vector calculated at the present iteration step
      * \param resBasis RingBuffer object storing the list of residual basis vectors.
      * \param nHist number of histories stored at this iteration
      */
      virtual void updateV(DArray<double> & v, T const & resCurrent, RingBuffer<T> const & resBasis, int nHist) = 0;
      
      /**
      * Set a field equal to another. Essentially a = b, but potentially more complex
      * in certain implementations of the AmIterator.
      * 
      * \param a the field to be set
      * \param b the field for it to be set to
      */
      virtual void setEqual(T& a, T const & b) = 0;

      /**
      * Mix histories, scaled by their respective coefficients, into the trial field.
      * 
      * \param trial object for calculation results to be stored in.
      * \param basis list of history basis vectors.
      * \param coeffs list of coefficients for each history.
      * \param nHist number of histories stored at this iteration
      */
      virtual void addHistories(T& trial, RingBuffer<T> const & basis, DArray<double> coeffs, int nHist) = 0;

      /**
      * Add predicted error into the field trial guess to attempt to correct for it.
      * 
      * \param fieldTrial field for calculation results to be stored in.
      * \param resTrial predicted error for current mixing of histories.
      * \param lambda Anderson-Mixing parameter for mixing in histories
      */
      virtual void addPredictedError(T& fieldTrial, T const & resTrial, double lambda) = 0;

      // ---- Methods for getting data from and sending data to the system ---- //
      
      /// Checks if the system has an initial guess
      virtual bool hasInitialGuess() = 0;
      
      /// Calculates and returns the number of elements in the
      /// array to be iterated
      virtual int nElements() = 0;

      /// Gets a reference to the current state of the system
      virtual void getCurrent(T& curr) = 0;

      /// Runs calculation to evaluate function for fixed point.
      virtual void evaluate() = 0;

      /// Gets residual values from system
      virtual void getResidual(T& resid) = 0;

      /// Updates the system with a passed in state of the iterator.
      virtual void update(T& newGuess) = 0;

      /// Outputs relevant system details to the iteration log
      virtual void outputToLog() = 0;

      // Members of parent classes with non-dependent names
      using Iterator::setClassName;
      using ParamComposite::read;
      using ParamComposite::readOptional;

   };

}
#include "AmIteratorTmpl.cpp"
#endif

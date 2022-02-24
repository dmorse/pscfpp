#ifndef PSCF_AM_ITERATOR_H
#define PSCF_AM_ITERATOR_H

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
   template <typename T>
   class AmIteratorTmpl : public ParamComposite 
   {
   public:

      /**
      * Constructor
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

   protected:
      
      using ParamComposite::read;
      using ParamComposite::readOptional;
      using ParamComposite::setClassName;

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
      * Check if ensemble is canonical. Returns false if grand-canonical 
      * or mixed ensemble.
      */
      bool isCanonical();

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

      // Pure virtual functions to be implemented by subclasses

      // System/iterator mediation classes
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

      // AmIterator mathematical operations
      /// Find the norm of the residual vector.
      virtual double findNorm(T const & hist) const = 0;

      /// Find the element of the residual vector with the maximum magnitude.
      virtual double findMaxAbs(T const & hist) const = 0;

      /// Update the list of residual basis vectors used for combining histories.
      virtual void updateBasis(RingBuffer<T> & basis, RingBuffer<T> const & hists) const = 0;

      /// Compute the dot product for constructing the U matrix. 
      virtual double computeUDotProd(RingBuffer<T> const & resBasis, int m, int n) const = 0;
      
      /// Compute the dot product for constructing the v vector. 
      virtual double computeVDotProd(T const & resCurrent, RingBuffer<T> const & resBasis, int m) const = 0;
      
      /// Update the U matrix containing dot products of residual histories basis vectors.
      virtual void updateU(DMatrix<double> & U, RingBuffer<T> const & resBasis, int nHist) const = 0;

      /// Update the v vector containing dot products of current residuals with residual basis
      /// vectors. 
      virtual void updateV(DArray<double> & v, T const & resCurrent, RingBuffer<T> const & resBasis, int nHist) const = 0;
      
      /// Set two things equal to each other.
      virtual void setEqual(T& a, T const & b) const = 0;

      /// Mix histories, scaled by their respective coefficients, into the trial field.
      virtual void addHistories(T& trial, RingBuffer<T> const & basis, DArray<double> coeffs, int nHist) const = 0;

      /// Add predicted error into the trial guess to approximately correct for it.
      virtual void addPredictedError(T& fieldTrial, T const & resTrial, double lambda) const = 0;

   };

}
#include "AmIteratorTmpl.cpp"
#endif

#ifndef PSCF_AM_ITERATOR_TMPL_H
#define PSCF_AM_ITERATOR_TMPL_H

/*
* PSCF - Polymer Self-Consistent Field Theory
*
* Copyright 2016 - 2022, The Regents of the University of Minnesota
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
   * Template for Anderson mixing iterator algorithm.
   *
   * Anderson mixing is an algorithm for solving a system of N nonlinear
   * equations of the form r_{i}(x) = 0 for i = 0, ..., N-1, where x
   * denotes a vector or array of unknown coordinate values. A vector 
   * of array of unknowns is referred here a "field" vector, while a 
   * vector or array of values of the errors r_{0},..., r_{N-1} is 
   * referred to as a residual vector.
   *
   * The type T is the type of the data structure used to store both 
   * field and residual vectors.
   *
   * \ingroup Pscf_Iterator_Module
   */
   template <typename Iterator, typename T>
   class AmIteratorTmpl : virtual public Iterator
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
      * Iterate to a solution
      *
      * \param isContinuation true iff continuation within a sweep
      * \return 0 for convergence, 1 for failure
      */
      int solve(bool isContinuation = false);

   protected:

      /// Type of error criterion used to test convergence 
      std::string errorType_;

      // Members of parent classes with non-dependent names
      using Iterator::setClassName;
      using ParamComposite::read;
      using ParamComposite::readOptional;

      /**
      * Allocate memory required by AM algorithm, if necessary.
      *
      * If the required memory has been allocated previously, this 
      * function does nothing and returns.
      */
      void allocateAM();

   private:

      // Private member variables

      /// Error tolerance.
      double epsilon_;

      /// Free parameter for minimization.
      double lambda_;

      /// Number of previous steps to use to compute next state. 
      int nHist_;

      /// Maximum number of previous states to retain.
      int maxHist_;

      /// Maximum number of iterations to attempt.
      int maxItr_;

      /// Number of elements in field or residual vectors.
      int nElem_; 

      /// Boolean indicating whether the calculation has diverged to NaN.
      bool diverged_;

      /// Has the allocate function been called.
      bool isAllocated_;

      /// History of previous field vectors.
      RingBuffer< T > fieldHists_;

      /// Basis vectors of field histories (differences of fields).
      RingBuffer< T > fieldBasis_;

      /// History of residual vectors.
      RingBuffer< T > resHists_;

      /// Basis vectors of residual histories (differences of residuals)
      RingBuffer< T > resBasis_;

      /// Matrix containing the dot products of vectors in resBasis_
      DMatrix<double> U_;

      /// Coefficients for mixing previous states.
      DArray<double> coeffs_;

      /// Dot products of current residual with residual basis vectors.
      DArray<double> v_;

      /// New trial field (big W in Arora et al. 2017)
      T fieldTrial_;

      /// Predicted field residual for trial state (big D)
      T resTrial_;

      /// Workspace for calculations
      T temp_;

      // --- Non-virtual private functions (implemented here) ---- //

      #if 0
      /**
      * Compute a vector of residuals, add to history.
      */
      void computeResidual();
      #endif

      /**
      * Compute the coefficients that would minimize residual L2 norm.
      */
      void computeResidCoeff();

      /**
      * Compute the trial field and set in the system.
      *
      * This implements a two-step process:
      *   - Correct the current state and predicted residual by adding 
      *     linear combination of state and residual basis vectors. 
      *   - Call addPredictedError to attempt to correct predicted 
      *     residual
      */
      void updateGuess();

      // --- Private virtual functions with default implementations --- //

      /**
      * Find the L2 norm of a vector (calls dotProduct internally).
      *
      * \param hist residual vector
      */
      virtual double norm(T const & hist);

      /**
      * Initialize just before entry to iterative loop.
      *
      * This function is called by the solve method just before entering
      * the loop over iterations. It must call the protected allocateAM()
      * to allocate memory required by the AM algorithm. The default
      * implementation just calls allocateAM(). Re-implementations by
      * subclasses may add additional operations.
      */ 
      virtual void setup();
     
      /**
      * Update the U matrix.
      * 
      * \param U U matrix
      * \param resBasis RingBuffer of residual basis vectors.
      * \param nHist number of histories stored at this iteration
      */      
      virtual void updateU(DMatrix<double> & U, 
                           RingBuffer<T> const & resBasis, int nHist);

      /**
      * Update the v vector.
      * 
      * \param v v vector
      * \param resCurrent  current residual vector 
      * \param resBasis RingBuffer of residual basis vectors.
      * \param nHist number of histories stored at this iteration
      */
      virtual void updateV(DArray<double> & v, T const & resCurrent, 
                           RingBuffer<T> const & resBasis, int nHist);
      
      /**
      * Check if solution is converge within specified tolerance.
      *
      * \return true if error < epsilon and false if error >= epsilon
      */
      virtual bool isConverged();

      /**
      * Clean up after a call to solve(), enabling future calls to solve.
      */
      virtual void cleanUp();

      // --- Pure virtual methods for doing AM iterator math --- //

      /**
      * Set one field equal to another.
      *
      * \param a the field to be set (lhs of assignment)
      * \param b the field value to assign (rhs of assignment)
      */
      virtual void setEqual(T& a, T const & b) = 0;

      /**
      * Compute the inner product of two vectors.
      *
      * \param a first vector
      * \param b second vector
      */
      virtual double dotProduct(T const & a, T const & b) = 0;

      /**
      * Find the maximum magnitude element of a residual vector.
      *
      * \param hist residual vector
      */
      virtual double maxAbs(T const & hist) = 0;

      /**
      * Update a basis that spans differences of past vectors.
      *
      * This function is applied to update bases for both residual
      * vectors and field vectors.
      * 
      * \param basis RingBuffer of residual basis vectors
      * \param hists RingBuffer of history of residual vectors
      */
      virtual 
      void updateBasis(RingBuffer<T> & basis, 
                       RingBuffer<T> const & hists) = 0;

      /**
      * Add linear combination of field basis vectors to the trial field.
      * 
      * \param trial object for calculation results to be stored in
      * \param basis list of history basis vectors
      * \param coeffs list of coefficients of basis vectors
      * \param nHist number of histories stored at this iteration
      */
      virtual 
      void addHistories(T& trial, 
                        RingBuffer<T> const & basis, 
                        DArray<double> coeffs, 
                        int nHist) = 0;

      /**
      * Remove predicted error from trial in attempt to correct for it.
      * 
      * \param fieldTrial field for calculation results to be stored in
      * \param resTrial predicted error for current mixing of histories
      * \param lambda Anderson-Mixing parameter for mixing in histories
      */
      virtual 
      void addPredictedError(T& fieldTrial, T const & resTrial, 
                             double lambda) = 0;

      // -- Pure virtual functions to exchange data with parent system -- //
     
      /** 
      * Does the system have an initial guess?
      */
      virtual bool hasInitialGuess() = 0;
     
      /** 
      * Compute and retur the number of residual or field vector elements.
      */
      virtual int nElements() = 0;

      /**
      * Get the current field vector from the system.
      *
      * \param curr current field vector (output)
      */
      virtual void getCurrent(T& curr) = 0;

      /**
      * Run a calculation to update the system state.
      */
      virtual void evaluate() = 0;

      /**
      * Compute residual vector.
      *
      * \param resid current residual vector.
      */
      virtual void getResidual(T& resid) = 0;

      /**
      * Update the system with a passed in trial field.
      *
      * \param newGuess new field vector (input)
      */
      virtual void update(T& newGuess) = 0;

      /**
      * Output relevant system details to the iteration log.
      */
      virtual void outputToLog() = 0;

   };
   
}
#include "AmIteratorTmpl.tpp"
#endif

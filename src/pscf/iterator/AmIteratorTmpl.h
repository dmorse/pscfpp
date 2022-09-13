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
   * Template for Anderson mixing iterator algorithm.
   *
   * Anderson mixing is an algorithm for solving a system of N nonlinear
   * equations of the form g_{i}(x) = 0 for i = 0, ..., N-1, where x
   * denotes a vector or array of unknown coordinate values. A vector or
   * array of unknowns is referred here a "field", while a vector or array 
   * of values of the errors g_{0},..., g_{N-1} is referred to as a residual.
   *
   * The type T is the type of the data structure used to store both field
   * residual vectors.
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
      * Setup and allocate required memory.
      */
      void setup();

      /**
      * Iterate to a solution
      *
      * \param isContinuation true iff continuation within a sweep
      * \return 0 for convergence, 1 for failure
      */
      int solve(bool isContinuation = false);

   protected:

      // Members of parent classes with non-dependent names
      using Iterator::setClassName;
      using ParamComposite::read;
      using ParamComposite::readOptional;

   private:

      /// Error tolerance.
      double epsilon_;

      /// Type of error checked for convergence (maxResid or normResid).
      std::string errorType_;

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

      /**
      * Compute a vector of residuals, add to history.
      */
      void computeResidual();

      /**
      * Check if solution is converge within specified tolerance.
      *
      * \return true if error < epsilon and false if error >= epsilon
      */
      bool isConverged();

      /**
      * Compute the coefficients that would minimize residual L2 norm.
      */
      void findResidCoeff();

      /**
      * Rebuild predicted wFields from minimized coefficients
      */
      void updateGuess();

      /**
      * Clean up after a call to solve(), enabling future calls to solve.
      */
      void cleanUp();

      // ---- Methods for doing AM iterator math ---- //

      /**
      * Find the L2 norm of a residual vector.
      *
      * \param hist residual vector
      */
      virtual double findNorm(T const & hist) = 0;

      /**
      * Find the maximum magnitude element of a residual vector.
      *
      * \param hist residual vector
      */
      virtual double findMaxAbs(T const & hist) = 0;

      /**
      * Update the basis of residual vectors.
      * 
      * \param basis RingBuffer of residual basis vectors
      * \param hists RingBuffer of history of residual vectors
      */
      virtual 
      void updateBasis(RingBuffer<T> & basis, 
                       RingBuffer<T> const & hists) = 0;

      /**
      * Update the U matrix.
      * 
      * \param U U matrix
      * \param resBasis RingBuffer of residual basis vectors.
      * \param nHist number of histories stored at this iteration
      */      
      virtual void updateU(DMatrix<double> & U, 
                           RingBuffer<T> const & resBasis, int nHist) = 0;

      /**
      * Update the v vector.
      * 
      * \param v v vector
      * \param resCurrent  current residual vector 
      * \param resBasis RingBuffer of residual basis vectors.
      * \param nHist number of histories stored at this iteration
      */
      virtual void updateV(DArray<double> & v, T const & resCurrent, 
                           RingBuffer<T> const & resBasis, int nHist) = 0;
      
      /**
      * Set one field equal to another.
      *
      * \param a the field to be set (lhs of assignment)
      * \param b the field value to assign (rhs of assignment)
      */
      virtual void setEqual(T& a, T const & b) = 0;

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

      // -- Functions to exchange data with the parent system - //
     
      /** 
      * Does the system has an initial guess?
      */
      virtual bool hasInitialGuess() = 0;
     
      /** 
      * Return the number of elements in a field or residual vector.
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

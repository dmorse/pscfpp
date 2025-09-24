#ifndef PSCF_AM_ITERATOR_TMPL_H
#define PSCF_AM_ITERATOR_TMPL_H

/*
* PSCF - Polymer Self-Consistent Field
*
* Copyright 2015 - 2025, The Regents of the University of Minnesota
* Distributed under the terms of the GNU General Public License.
*/

#include <util/containers/DArray.h>       // member template
#include <util/containers/DMatrix.h>      // member template
#include <util/containers/RingBuffer.h>   // member template
#include <util/accumulators/Average.h>    // member template
#include <util/misc/Timer.h>              // member

// Uncomment to test details of Anderson-Mixing algorithm performance
//#define PSCF_AM_TEST

namespace Pscf {

   using namespace Util;

   /**
   * Template for Anderson mixing iterator algorithm.
   *
   * Anderson mixing is an algorithm for solving a system of N nonlinear
   * equations of the form R{i}(X) = 0 for i = 0, ..., N-1, where X
   * denotes a vector or array of N unknown coordinate values. A vector
   * of array of unknowns is referred here as the state vector, while 
   * a vector or array of values of the errors R{0} ,..., R{N-1} is
   * referred to as the residual vector.
   *
   * The template type parameter Iterator is a base class that must 
   * be derived from Util::ParamComposite. In applications to SCFT 
   * iterators, the Iterator class may declare a virtual solve function
   * that is overridden by the solve() function defined by this template.
   *
   * The template type parameter T is the type of data structure that is 
   * used to represent state and residual vectors.
   *
   * \sa \ref pscf_AmIteratorTmpl_page
   *
   * \ingroup Pscf_Iterator_Module
   */
   template <typename Iterator, typename T>
   class AmIteratorTmpl : virtual public Iterator
   {
   public:

      /**
      * Constructor.
      */
      AmIteratorTmpl();

      /**
      * Destructor.
      */
      ~AmIteratorTmpl();

      /**
      * Read all parameters and initialize.
      *
      * \param in input filestream
      */
      virtual void readParameters(std::istream& in);

      /**
      * Iterate to a solution.
      *
      * \param isContinuation true iff continuation within a sweep
      * \return 0 for convergence, 1 for failure
      */
      int solve(bool isContinuation = false);

      /**
      * Log output timing results.
      * 
      * \param out  output stream
      */
      void outputTimers(std::ostream& out) const;

      /**
      * Clear timers
      */
      void clearTimers();

   protected:

      /**
      * Error tolerance.
      */
      double epsilon_;

      /**
      * Maximum number of iterations to attempt.
      */
      int maxItr_;

      /**
      * Maximum number of basis vectors in AM algorithm.
      */
      int maxHist_;

      /**
      * Verbosity level.
      */
      int verbose_;

      /**
      * Type of error criterion used to test convergence.
      */
      std::string errorType_;

      // Parameter initialization

      /**
      * Read and validate the optional errorType string parameter.
      *
      * \param in  input filestream
      */
      void readErrorType(std::istream& in);

      /**
      * Checks if a string is a valid error type.
      *
      * Virtual to allow extension of allowed error type string values.
      *
      * \return true  if type is valid, false otherwise.
      */
      virtual bool isValidErrorType();

      /**
      * Read optional parameters used in default correction algorithm.
      *
      * Sets a default value for useLambdaRamp, and then optionally
      * reads lambda, useLambdaRamp, and (if useLambdaRamp == true)
      * the ratio r used in the ramp.
      *
      * \param in  input filestream
      * \param useLambdaRamp  default value for useLambdaRamp
      */
      void readMixingParameters(std::istream& in, 
                                bool useLambdaRamp = true);

      // Protected AM mixing operations

      /**
      * Allocate memory required by AM algorithm, if necessary.
      *
      * If the required memory has been allocated previously, this
      * function does nothing and returns.
      */
      void allocateAM();

      /**
      * Have data structures required by the AM algorithm been allocated?
      */
      bool isAllocatedAM() const;

      /**
      * Initialize just before entry to iterative loop.
      *
      * This function is called by the solve function before entering the
      * loop over iterations. The default implementation calls allocateAM()
      * if isAllocatedAM() is false, and clears state and residual 
      * histories. If isContinuation is false, it also clears state and
      * residual basis vector lists. 
      *
      * \param isContinuation  true iff continuation within a sweep
      */
      virtual void setup(bool isContinuation);

      /**
      * Clear information about history.
      *
      * This function clears the the history and basis vector ring
      * buffer containers.
      */
      virtual void clear();

      /**
      * Compute and return error used to test for convergence.
      *
      * \param residTrial  current residual vector
      * \param stateTrial  current state vector
      * \param errorType  type of error
      * \param verbose  verbosity level of output report.
      * \return error  measure used to test for convergence.
      */
      virtual 
      double computeError(T& residTrial, T& stateTrial,
                          std::string errorType, int verbose);

      /**
      * Compute and return error used to test for convergence.
      *
      * \param verbose  verbosity level of output report
      * \return error  error value used to test for convergence.
      */
      double computeError(int verbose);

      /**
      * Compute ramped prefactor of mixing parameter lambda.
      *
      * \return prefactor 1 - r^{nBasis} multiplying lambda.
      */
      double lambdaRampFactor();

      /**
      * Compute mixing parameter for correction step of Anderson mixing.
      *
      * If useLambdaRamp_ == true and nBasis_ < maxHist_, then return the
      * ramped parameter lambda_ * ( 1 - r_^(nBasis) ). Otherwise, return
      * the parameter lambda_.
      *
      * \return lambda mixing parameter
      */
      virtual double computeLambda();

      // Protected accessors for member variables

      /**
      * Verbosity level, allowed values 0, 1, or 2.
      */
      int verbose() const;

      /**
      * Get error type string.
      */
      std::string errorType() const;

      /**
      * Get the current residual vector by const reference.
      */
      T const & residual() const;

      /**
      * Return the current state vector by const reference.
      */
      T const & state() const;

      /**
      * Return the total number of iterations needed to converge.
      */
      int totalItr();

      // --- Timer value accessors --- //

      /**
      * Get total time.
      */
      double timerTotal();

      /**
      * Get time spent solving the modified diffusion equation (MDE).
      */
      double timerMDE();

      /**
      * Get total time for AM algorithm, excluding MDE solution.
      */
      double timerAM();

      /**
      * Get time spent computing residual.
      */
      double timerResid();

      /**
      * Get time evaluating scalar error.
      */
      double timerError();

      /**
      * Get time spent evaluating Anderson mixing coefficients.
      */
      double timerCoeff();

      /**
      * Get time spent updating w states.
      */
      double timerOmega();

      // Inherited members of parent classes with non-dependent names
      using ParamComposite::read;
      using ParamComposite::readOptional;

   private:

      // Private member variables

      /// Current scalar error.
      double error_;

      /// Number of basis vectors defined as differences.
      int nBasis_;

      /// Current iteration counter.
      int itr_;

      /// Total iteration counter.
      int totalItr_;

      /// Number of elements in state or residual vectors.
      int nElem_;

      /// Mixing coefficient for standard AM mixing / correction step.
      double lambda_;

      /// Ramp parameter for ramped AM mixing / correction step.
      double r_;

      /// Should the AM mixing coefficient lambda be ramped up ?
      bool useLambdaRamp_;

      /// Has the allocateAM function been called?
      bool isAllocatedAM_;

      /// History of previous state vectors.
      RingBuffer<T> stateHistory_;

      /// Basis vectors for state (differences of state vectors).
      RingBuffer<T> stateBasis_;

      /// History of previous residual vectors.
      RingBuffer<T> residualHistory_;

      /// Basis vectors for residuals (differences of residual vectors).
      RingBuffer<T> residualBasis_;

      /// Matrix containing the dot products of vectors in residualBasis_
      DMatrix<double> U_;

      /// Dot products of current residual with residual basis vectors.
      DArray<double> v_;

      /// Coefficients of basis vectors that minimize predicted residual.
      DArray<double> coeffs_;

      /// New trial state vector (first created by updateTrial)
      T stateTrial_;

      /// Predicted residual for trial state (created by updateTrial)
      T residualTrial_;

      /// Workspace for calculations
      T temp_;

      // Timers for analyzing performance
      Timer timerMDE_;
      Timer timerAM_;
      Timer timerResid_;
      Timer timerError_;
      Timer timerCoeff_;
      Timer timerOmega_;
      Timer timerTotal_;

      #ifdef PSCF_AM_TEST
      bool hasAmTest_{false};
      double preError_{0};
      double projectionError_{0};
      double correctionError_{0};
      double projectionRatio_{0};
      double correctionRatio_{0};
      int testCounter{0};
      #endif

      // Private non-virtual functions used in AM algorithm

      /**
      * Update a basis that spans sequential differences of past vectors.
      *
      * This function is applied to update bases for both state vectors
      * and residual vectors.
      *
      * \param basis  RingBuffer of basis vectors
      * \param history  RingBuffer of history of prior vectors
      */
      void updateBasis(RingBuffer<T> & basis,
                       RingBuffer<T> const & history);

      /**
      * Update the U matrix.
      *
      * \param U U matrix
      * \param resBasis RingBuffer of residual basis vectors.
      */
      void updateU(DMatrix<double> & U,
                   RingBuffer<T> const & resBasis);

      /**
      * Compute the v vector.
      *
      * \param v  v vector
      * \param resCurrent  current residual vector
      * \param resBasis  RingBuffer of residual basis vectors.
      */
      void computeV(DArray<double> & v, 
                   T const & resCurrent,
                   RingBuffer<T> const & resBasis);

      /**
      * Compute coefficients for the trial state vector.
      *
      * Computes coefficients of basis vectors chosen to minimize the l2
      * norm of the predicted residual vector. Resulting coefficients are
      * stored in the private member variable coeffs_
      */
      void computeTrialCoeff();

      /**
      * Compute trial state and predicted trial residual vectors.
      * 
      * Adds linear combinations of basis vectors with coefficients 
      * computed by the computeTrialCoeff() functions to the current 
      * state vector and residual vector to obtain a trial state vector 
      * and a corresponding predicted trial residual vector. Resulting
      * vectors are stored in private member variables stateTrial_ and 
      * residualTrial_.
      */
      void updateTrial();

      /**
      * Add linear combination of basis vectors to a vector.
      *
      * This function is used within updateTrial to update both state 
      * and residual vectors, using different bases but the same list 
      * of coefficients.
      *
      * \param v  vector to be modified (input / output)
      * \param basis  list of basis vectors (input)
      * \param coeffs  list of coefficients of basis vectors (input)
      */
      void addEqVectors(T& v,
                        RingBuffer<T> const & basis,
                        DArray<double> coeffs);

      // Private virtual functions with a default implementation 

      /**
      * Add correction based on the residual vector. 
      *
      * This is the second "correction" stage of an Anderson mixing
      * algorithm. The default implementation simply adds a correction 
      * proportional to the predicted residual vector, residualTrial, 
      * multiplied by the coefficient given by the computeLambda() 
      * function.
      *
      * \param stateTrial  trial state vector (in/out)
      * \param residualTrial  predicted residual for trial state (in)
      */
      virtual 
      void addCorrection(T& stateTrial, T const & residualTrial);

      // Private pure virtual functions that interact with parent system

      /**
      * Compute and return the number of residual or state vector elements.
      *
      * The private variable nElem_ is assigned the return value of this
      * function on entry to allocateAM to set the number of elements of
      * the residual and state vectors.
      */
      virtual int nElements() = 0;

      /**
      * Does the system have an initial guess for the state vector?
      */
      virtual bool hasInitialGuess() = 0;

      /**
      * Get the current state vector from the system.
      *
      * \param curr  current state vector (output)
      */
      virtual void getCurrent(T& curr) = 0;

      /**
      * Run a calculation to update the system state.
      *
      * This function normally solves the modified diffusion equations,
      * calculates monomer concentration states, and computes stresses
      * if appropriate.
      */
      virtual void evaluate() = 0;

      /**
      * Compute the residual vector from the current system state.
      *
      * \param resid current residual vector.
      */
      virtual void getResidual(T& resid) = 0;

      /**
      * Update the system with a passed in trial state vector.
      *
      * \param newGuess new state vector (input)
      */
      virtual void update(T& newGuess) = 0;

      /**
      * Output relevant system details to the iteration log.
      */
      virtual void outputToLog() = 0;

      // Private functions for vector math

      /**
      * Assignment for vectors of type T.
      *
      * This function must perform an assignment a = b.
      *
      * \param a  vector to be set (lhs of assignment)
      * \param b  vector value to assign (rhs of assignment)
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
      * Find the L2 norm of a vector.
      *
      * This calls dotProduct internally, returning the square root of 
      * dotProduct(a, a).
      *
      * \param a residual vector
      */
      double norm(T const & a);

      /**
      * Return the maximum magnitude element of a vector.
      *
      * \param a  input vector
      */
      virtual double maxAbs(T const & a) = 0;

      /**
      * Compute the difference a = b - c for vectors a, b and c.
      *
      * \param a result vector (LHS)
      * \param b first vector (RHS)
      * \param c second vector (RHS)
      */
      virtual void subVV(T& a, T const & b, T const & c) = 0;

      /**
      * Compute a += c*b for vectors a and b and scalar c.
      *
      * \param a result vector (LHS)
      * \param b input vector (RHS)
      * \param c scalar coefficient (RHS)
      */
      virtual void addEqVc(T& a, T const & b, double c) = 0;

   };

   /*
   * Return integer level for verbosity of the log output (0-2).
   */
   template <typename Iterator, typename T>
   int AmIteratorTmpl<Iterator,T>::verbose() const
   {  return verbose_; }

   /*
   * Return error type string.
   */
   template <typename Iterator, typename T>
   std::string AmIteratorTmpl<Iterator,T>::errorType() const
   {  return errorType_; }

   /*
   * Return the current state vector by const reference.
   */
   template <typename Iterator, typename T>
   T const & AmIteratorTmpl<Iterator,T>::state() const
   {  return stateHistory_[0]; }

   /*
   * Return the current residual vector by const reference.
   */
   template <typename Iterator, typename T>
   T const & AmIteratorTmpl<Iterator,T>::residual() const
   {  return residualHistory_[0]; }

   /*
   * Has memory required by the AM algorithm been allocated?
   */
   template <typename Iterator, typename T>
   bool AmIteratorTmpl<Iterator,T>::isAllocatedAM() const
   {  return isAllocatedAM_; }

   /*
   * Return total iteration counter.
   */
   template <typename Iterator, typename T>
   int AmIteratorTmpl<Iterator,T>::totalItr()
   {  return totalItr_; }

   /*
   * Return computing MDE time cost.
   */
   template <typename Iterator, typename T>
   double AmIteratorTmpl<Iterator,T>::timerMDE()
   {  return timerMDE_.time(); }

   /*
   * Return computing AM time cost.
   */
   template <typename Iterator, typename T>
   double AmIteratorTmpl<Iterator,T>::timerAM()
   {  return timerAM_.time(); }

   /*
   * Return computing Resid time cost.
   */
   template <typename Iterator, typename T>
   double AmIteratorTmpl<Iterator,T>::timerResid()
   {  return timerResid_.time(); }

   /*
   * Return computing Error time cost.
   */
   template <typename Iterator, typename T>
   double AmIteratorTmpl<Iterator,T>::timerError()
   {  return timerError_.time(); }

   /*
   * Return computing Coeff time cost.
   */
   template <typename Iterator, typename T>
   double AmIteratorTmpl<Iterator,T>::timerCoeff()
   {  return timerCoeff_.time(); }

   /*
   * Return computing Omega time cost.
   */
   template <typename Iterator, typename T>
   double AmIteratorTmpl<Iterator,T>::timerOmega()
   {  return timerOmega_.time(); }

   /*
   * Return total time cost.
   */
   template <typename Iterator, typename T>
   double AmIteratorTmpl<Iterator,T>::timerTotal()
   {  return timerTotal_.time(); }

}
#include "AmIteratorTmpl.tpp"
#endif

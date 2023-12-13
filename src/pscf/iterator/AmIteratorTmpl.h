#ifndef PSCF_AM_ITERATOR_TMPL_H
#define PSCF_AM_ITERATOR_TMPL_H
/*
* Uncomment to test the contirbution of Anderson-Mixing for error reduction
* from linear mixing step and correction step 
*/ 
//#define PSCF_AM_TEST

/*
* PSCF - Polymer Self-Consistent Field Theory
*
* Copyright 2016 - 2022, The Regents of the University of Minnesota
* Distributed under the terms of the GNU General Public License.
*/

#include <util/containers/DArray.h>     // member template
#include <util/containers/DMatrix.h>    // member template
#include <util/containers/RingBuffer.h> // member template
#include <util/misc/Timer.h>            // member template
#include <util/accumulators/Average.h>  // member template

namespace Pscf {

   using namespace Util;
   
   /**
   * Template for Anderson mixing iterator algorithm.
   *
   * Anderson mixing is an algorithm for solving a system of N nonlinear
   * equations of the form r_{i}(x) = 0 for i = 0, ..., N-1, where x
   * denotes a vector or array of unknown coordinate values. A vector 
   * of array of unknowns is referred here the "field" vector, while a 
   * vector or array of values of the errors r_{0},..., r_{N-1} is 
   * referred to as the residual vector.
   *
   * The template parameter Iterator is a base class that must be derived
   * from Util::ParamComposite, and must declare a virtual solve(bool) 
   * function with the same interface as that declared here.
   *
   * The template type parameter T is the type of the data structure used 
   * to store both field and residual vectors. 
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
      
      /**
      * Log output timing results 
      */
      void outputTimers(std::ostream& out);
      
      /**
      * Clear timers 
      */
      void clearTimers();
      
      /**
      * Obtain error type
      */
      std::string errorType();

   protected:

      /// Type of error criterion used to test convergence 
      std::string errorType_;

      // Members of parent classes with non-dependent names
      using Iterator::setClassName;
      using ParamComposite::read;
      using ParamComposite::readOptional;

      /**
      * Set value of maxItr.
      *
      * Provided to allow subclasses to set a modified default value 
      * before calling readParameters, in which maxItr is optional.
      * Global default, set in constructor, is maxItr = 200.
      *
      * \param maxItr  maximum number of iterations attempted
      */
      void setMaxItr(int maxItr);

      /**
      * Set value of maxHist (number of retained previous states)
      *
      * Provided to allow subclasses to set a modified default value 
      * before calling readParameters, in which maxItr is optional.
      * Global default, set in constructor, is maxHist = 50.
      *
      * \param maxHist  maximum number of retained previous states
      */
      void setMaxHist(int maxHist);

      /**
      * Set and validate value of errorType string.
      *
      * Provided to allow subclasses to set a modified default value 
      * before calling readParameters, in which errorType is optional.
      * Global default, set in constructor, is relNormResid = 50.
      *
      * \param errorType error type string
      */
      void setErrorType(std::string errorType);

      /**
      * Read and validate the optional errorType string parameter.
      *
      * \param in input filestream
      */
      void readErrorType(std::istream& in);

      /**
      * Checks if a string is a valid error type.
      *
      * Virtual to allow extension of allowed error type string values.
      *
      * \return true if type is valid, false otherwise.
      */
      virtual bool isValidErrorType();

      /**
      * Find the L2 norm of a vector. 
      *
      * The default implementation calls dotProduct internally. 
      * Virtual to allow more optimized versions.
      *
      * \param hist residual vector
      */
      virtual double norm(T const & hist);

      /**
      * Allocate memory required by AM algorithm, if necessary.
      *
      * If the required memory has been allocated previously, this 
      * function does nothing and returns.
      */
      void allocateAM();

      /**
      * Clear information about history.
      *
      * This function clears the the history and basis vector ring 
      * buffer containers.
      */
      virtual void clear();

      /**
      * Initialize just before entry to iterative loop.
      *
      * This function is called by the solve function before entering the
      * loop over iterations. The default functions calls allocateAM() 
      * if isAllocatedAM() is false, and otherwise calls clear() if 
      * isContinuation is false.
      *
      * \param isContinuation true iff continuation within a sweep
      */ 
      virtual void setup(bool isContinuation);
     
      /**
      * Compute and return error used to test for convergence.
      *
      * \param verbose  verbosity level of output report
      * \return error  measure used to test for convergence.
      */
      virtual double computeError(int verbose);
      
      /**
      * Set mixing parameter for second step of an Anderson Mixing Algorithm
      *
      * \return lambda mixing parameter
      */
      virtual double setLambda();
      
      #ifdef PSCF_AM_TEST
      double computeError(T a);
      #endif

      /**
      * Return the current residual vector by const reference.
      */
      T const & residual() const;

      /**
      * Return the current field or state vector by const reference.
      */
      T const & field() const;

      /**
      * Verbosity level, allowed values 0, 1, or 2.
      */
      int verbose() const;
      
      /**
      * Return the total number of iterations needed to converge.
      */
      int totalItr();
      
      /// Get total time.
      double timerTotal();
      
      /// Get time solving modified diffusion equation (MDE).
      double timerMDE();

      /// Get total time for AM algorithm, excluding MDE solution.
      double timerAM();

      /// Get time computing residual.
      double timerResid();

      /// Get time evaluating scalar error.
      double timerError();

      /// Get time evaluating Anderson mixing coefficients.
      double timerCoeff();

      /// Get time updating w fields.
      double timerOmega();

      /**
      * Have data structures required by the AM algorithm been allocated?
      */
      bool isAllocatedAM() const;
      
   private:
      
      #ifdef PSCF_AM_TEST
      double preError_{0};
      double mixingError_{0};
      double correctionError_{0};
      double mixingRatio_{0};
      double correctionRatio_{0};
      int testCounter{0};
      #endif
      
      // Private member variables
      /// Error
      double error_;
      /// Error tolerance.
      double epsilon_;

      /// Free parameter for minimization.
      double lambda_;

      /// Maximum number of iterations to attempt.
      int maxItr_;

      /// Maximum number of basis vectors
      int maxHist_;

      /// Number of basis vectors defined as differences.
      int nBasis_;

      /// Current iteration counter.
      int itr_;
      
      /// Total iteration counter.
      int totalItr_;

      /// Number of elements in field or residual vectors.
      int nElem_; 

      /// Verbosity level.
      int verbose_;

      /// Has the allocateAM function been called.
      bool isAllocatedAM_;

      /// History of previous field vectors.
      RingBuffer<T> fieldHists_;

      /// Basis vectors of field histories (differences of fields).
      RingBuffer<T> fieldBasis_;

      /// History of residual vectors.
      RingBuffer<T> resHists_;

      /// Basis vectors of residual histories (differences of residuals)
      RingBuffer<T> resBasis_;

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
      
      // Timers for analyzing performance
      Timer timerMDE_;
      Timer timerAM_;
      Timer timerResid_;
      Timer timerError_;
      Timer timerCoeff_;
      Timer timerOmega_;
      Timer timerTotal_;

      // --- Non-virtual private functions (implemented here) ---- //

      /**
      * Compute optimal coefficients of residual basis vectors.
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
      
      // --- Pure virtual methods for doing AM iterator math --- //

      /**
      * Set one vector equal to another.
      *
      * \param a the vector to be set (lhs of assignment)
      * \param b the vector value to assign (rhs of assignment)
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
      * Compute and return the number of residual or field vector elements.
      *
      * The private variable nElem_ is assigned the return value of this 
      * function on entry to allocateAM to set the number of elements of
      * the residual and field vectors.
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
      *
      * This function normally solves the modified diffusion equations,
      * calculates monomer concentration fields, and computes stresses 
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
      * Update the system with a passed in trial field vector.
      *
      * \param newGuess new field vector (input)
      */
      virtual void update(T& newGuess) = 0;

      /**
      * Output relevant system details to the iteration log.
      */
      virtual void outputToLog() = 0;

   };
  
   /*
   * Return the current residual vector by const reference.
   */
   template <typename Iterator, typename T>
   T const & AmIteratorTmpl<Iterator,T>::residual() const
   {  return resHists_[0]; }

   /*
   * Return the current field/state vector by const reference.
   */
   template <typename Iterator, typename T>
   T const & AmIteratorTmpl<Iterator,T>::field() const
   {  return fieldHists_[0]; }

   /*
   * Integer level for verbosity of the log output (0-2).
   */
   template <typename Iterator, typename T>
   int AmIteratorTmpl<Iterator,T>::verbose() const
   {  return verbose_; }

   /*
   * Has memory required by AM algorithm been allocated?
   */
   template <typename Iterator, typename T>
   bool AmIteratorTmpl<Iterator,T>::isAllocatedAM() const
   {  return isAllocatedAM_; }
   
   /*
   * Return error type
   */ 
   template <typename Iterator, typename T>
   std::string AmIteratorTmpl<Iterator,T>::errorType() 
   {  return errorType_; }
   
   /*
   * Return total iteration counter
   */
   template <typename Iterator, typename T>
   int AmIteratorTmpl<Iterator,T>::totalItr() 
   {  return totalItr_; }
   
   /*
   * Return computing MDE time cost
   */
   template <typename Iterator, typename T>
   double AmIteratorTmpl<Iterator,T>::timerMDE() 
   {  return timerMDE_.time(); }
  
   /*
   * Return computing AM time cost
   */
   template <typename Iterator, typename T>
   double AmIteratorTmpl<Iterator,T>::timerAM() 
   {  return timerAM_.time(); }
   
   /*
   * Return computing Resid time cost
   */
   template <typename Iterator, typename T>
   double AmIteratorTmpl<Iterator,T>::timerResid() 
   {  return timerResid_.time(); }
   
   /*
   * Return computing Error time cost
   */
   template <typename Iterator, typename T>
   double AmIteratorTmpl<Iterator,T>::timerError() 
   {  return timerError_.time(); }
   
   /*
   * Return computing Coeff time cost
   */
   template <typename Iterator, typename T>
   double AmIteratorTmpl<Iterator,T>::timerCoeff() 
   {  return timerCoeff_.time(); }
   
   /*
   * Return computing Omega time cost
   */
   template <typename Iterator, typename T>
   double AmIteratorTmpl<Iterator,T>::timerOmega() 
   {  return timerOmega_.time(); }
   
   /*
   * Return total time cost
   */
   template <typename Iterator, typename T>
   double AmIteratorTmpl<Iterator,T>::timerTotal() 
   {  return timerTotal_.time(); }
   
}
#include "AmIteratorTmpl.tpp"
#endif

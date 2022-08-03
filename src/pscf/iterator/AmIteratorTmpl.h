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
      */
      int solve();

      /**
       * Accepts a field array representing the concentration profile
       * of the mask that will be imposed upon the unit cell, and stores 
       * it for use during iteration.
       * 
       * \param field mask concentration profile, in basis format
       */
      void setMask(T const & field);

      /**
       * Accepts an array of fields representing the external potential
       * field felt by each monomer species, and stores it for use 
       * during iteration.
       * 
       * \param fields mask concentration profile, in basis format
       */
      void setExternalFields(DArray<T> const & fields);

      /**
      * Returns the field that represents the mask imposed upon the 
      * unit cell, or returns an unallocated field if the iterator 
      * does not have a mask. Field is in symmetry-adapted basis format.
      */
      T const & maskField() const;

      /**
      * Returns the array of fields that represents the external potential
      * field felt by each monomer species, or returns an unallocated field
      * if the iterator does not have an external field. Fields are in 
      * symmetry-adapted basis format.
      */
      DArray<T> const & externalFields() const;

      /**
      * Returns the field that represents the external field imposed
      * on the monomer species monomerId, or returns an unallocated
      * field if the iterator does not have an external field. Field
      * is in symmetry-adapted basis format.
      * 
      * \param monomerId integer monomer type index
      */
      T const & externalField(int monomerId) const;

      /**
      * Return true if this iterator imposes a mask within the unit cell,
      * false if no mask is imposed.
      */
      bool const hasMask() const;

      /**
      * Return true if this iterator imposes an external field on any
      * monomer species, returns false if no external field is present.
      */
      bool const hasExternalField() const;

   protected:

      // Members of parent classes with non-dependent names
      using Iterator::setClassName;
      using ParamComposite::read;
      using ParamComposite::readOptional;

   private:

      /// Error tolerance
      double epsilon_;

      /// Type of error checked for convergence (maxResid or normResid)
      std::string errorType_;

      /// Free parameter for minimization
      double lambda_;

      /// Number of previous steps to use to compute next state. [0,maxHist_]
      int nHist_;

      /// Number of histories to retain
      int maxHist_;

      /// Maximum number of iterations to attempt
      int maxItr_;

      /// Number of elements in the field
      int nElem_; 

      /// Field histories
      RingBuffer< T > fieldHists_;

      /// Basis vectors of field histories
      RingBuffer< T > fieldBasis_;

      /// Residual histories
      RingBuffer< T > resHists_;

      /// Basis vectors of residual histories (differences of residuals)
      RingBuffer< T > resBasis_;

      /// Matrix containing the dot products of vectors in resBasis_
      DMatrix<double> U_;

      /// Coefficients for mixing previous states
      DArray<double> coeffs_;

      /// Dot products of current residual with residual basis vectors
      DArray<double> v_;

      /// New trial field (big W in Arora et al. 2017)
      T fieldTrial_;

      /// Predicted field residual for trial state (big D)
      T resTrial_;

      /// Workspace for calculations
      T temp_;

      /// Field representing the mask imposed on the unit cell
      DArray<double> maskField_;

      /// Array containing the external fields felt by each species
      DArray<DArray<double> > externalFields_;

      /// Does this iterator have a mask?
      bool hasMask_;

      /// Does this iterator have an external field?
      bool hasExternalField_;

      /**
      * Compute the deviation of wFields from a mean field solution
      */
      void computeResidual();

      /**
      * Check if solution is converge within specified tolerance.
      *
      * \return true if error < epsilon and false if error >= epsilon
      */
      bool isConverged();

      /**
      * Compute the coefficients that would minimize residual norm.
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
      * Find norm of a residual vector.
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
      * Update the series of residual vectors.
      * 
      * \param basis RingBuffer storing the list of residual basis vectors
      * \param hists RingBuffer storing the histories of residual vectors
      */
      virtual 
      void updateBasis(RingBuffer<T> & basis, 
                       RingBuffer<T> const & hists) = 0;

      /**
      * Compute the dot product for an element of the U matrix.
      * 
      * \param resBasis RingBuffer storing residual basis vectors.
      * \param m row index for U matrix
      * \param n column index of the U matrix
      */
      virtual 
      double computeUDotProd(RingBuffer<T> const & resBasis, 
                             int m, int n) = 0;
      
      /**
      * Compute the dot product for an element of the v vector.
      * 
      * \param resCurrent  current residual vector 
      * \param resBasis RingBuffer of past residual basis vectors.
      * \param m row index of element of v vector
      */
      virtual 
      double computeVDotProd(T const & resCurrent, 
                             RingBuffer<T> const & resBasis, int m) = 0;
      
      /**
      * Compute required dot products and update the U matrix.
      * 
      * \param U U matrix
      * \param resBasis RingBuffer of residual basis vectors.
      * \param nHist number of histories stored at this iteration
      */      
      virtual void updateU(DMatrix<double> & U, 
                           RingBuffer<T> const & resBasis, int nHist) = 0;

      /**
      * Compute required dot products and update the v vector.
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
      * Essentially a = b, but potentially more complex in some 
      * implementations of the AmIterator.
      * 
      * \param a the field to be set (value)
      * \param b the field value to assign (rhs)
      */
      virtual void setEqual(T& a, T const & b) = 0;

      /**
      * Mix histories, scaled by their respective coefficients, into the trial field.
      * 
      * \param trial object for calculation results to be stored in
      * \param basis list of history basis vectors
      * \param coeffs list of coefficients for each history
      * \param nHist number of histories stored at this iteration
      */
      virtual 
      void addHistories(T& trial, RingBuffer<T> const & basis, 
                        DArray<double> coeffs, int nHist) = 0;

      /**
      * Remove predicted from trial in attempt to correct for it.
      * 
      * \param fieldTrial field for calculation results to be stored in
      * \param resTrial predicted error for current mixing of histories
      * \param lambda Anderson-Mixing parameter for mixing in histories
      */
      virtual 
      void addPredictedError(T& fieldTrial, T const & resTrial, 
                             double lambda) = 0;

      // -- Functions to exchange data with the parent system - //
      
      /// Checks if the system has an initial guess
      virtual bool hasInitialGuess() = 0;
      
      /// Return the number of elements in the array or residual
      virtual int nElements() = 0;

      /// Get a reference to the current state of the system
      virtual void getCurrent(T& curr) = 0;

      /// Run a calculation to evaluate function for current state.
      virtual void evaluate() = 0;

      /// Get a residual values from system
      virtual void getResidual(T& resid) = 0;

      /// Update the system with a passed in state of the iterator.
      virtual void update(T& newGuess) = 0;

      /// Output relevant system details to the iteration log
      virtual void outputToLog() = 0;

   };

   // Inline member functions

   // Store the concentration field for the mask imposed on the unit cell
   template <typename Iterator, typename T>
   inline void AmIteratorTmpl<Iterator,T>::setMask(T const & field)
   {  
      maskField_ = field; // copy field into maskField_
      hasMask_ = true;
   }

   // Store the external fields imposed on each monomer species
   template <typename Iterator, typename T>
   inline 
   void AmIteratorTmpl<Iterator,T>::setExternalFields(DArray<T> const & fields)
   {  
      externalFields_ = fields; // copy fields into externalFields_
      hasExternalField_ = true;
   }

   // Get the concentration field for the mask imposed on the unit cell
   template <typename Iterator, typename T>
   inline T const & AmIteratorTmpl<Iterator,T>::maskField() const
   {  return maskField_; }

   // Get array of all external fields felt by monomer species
   template <typename Iterator, typename T>
   inline DArray<T> const & AmIteratorTmpl<Iterator,T>::externalFields() 
   const
   {  return externalFields_; }

   // Get the external field felt by one monomer species
   template <typename Iterator, typename T>
   inline 
   T const & AmIteratorTmpl<Iterator,T>::externalField(int monomerId) 
   const
   {  return externalFields_[monomerId]; }

   // Does this iterator have a mask?
   template <typename Iterator, typename T>
   inline bool const AmIteratorTmpl<Iterator,T>::hasMask() const
   {  return hasMask_; }

   // Does this iterator have any external field?
   template <typename Iterator, typename T>
   inline bool const AmIteratorTmpl<Iterator,T>::hasExternalField() const
   {  return hasExternalField_; }

}
#include "AmIteratorTmpl.tpp"
#endif

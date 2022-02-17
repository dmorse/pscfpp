#ifndef PSCF_AM_ITERATOR_H
#define PSCF_AM_ITERATOR_H

/*
* PSCF - Polymer Self-Consistent Field Theory
*
* Copyright 2016 - 2019, The Regents of the University of Minnesota
* Distributed under the terms of the GNU General Public License.
*/

#include "Iterator.h" // base class

#include <pscf/math/LuSolver.h>
#include "AmStrategy.h"
#include "IteratorMediator.h"

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
   class AmIterator : public Iterator<T>
   {
   public:

      /**
      * Constructor
      *
      * \param iterMed pointer to an iterator mediator
      */
      AmIterator(IteratorMediator<T>& iterMed, AmStrategy<T>& strategy);

      /**
      * Destructor
      */
      ~AmIterator();

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
      * Return the strategy object to manage math execution.
      */ 
      AmStrategy<T>& strategy();

   private:

      /// Strategy object for math execution.
      AmStrategy<T>* strategy_;

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

      /// History of field residuals.
      /// 1st index = history, 2nd index = monomer, 3rd index = basis func.
      RingBuffer< T > resHists_;

      /// History of previous fields
      RingBuffer< T > fieldHists_;

      /// Matrix, the dot products of residual differences.
      DMatrix<double> U_;

      /// Cn, coefficients for mixing previous states.
      DArray<double> coeffs_;

      /// Vector, dot products of residuals with differences from histories
      DArray<double> v_;

      /// New trial field (big W in Arora et al. 2017)
      T fieldTrial_;

      /// Predicted field residual for trial state (big D)
      T resTrial_;

      /// Workspace for extracting field.
      T fieldTemp_;

      /// Workspace for residual calculation.
      T resTemp_;

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

      // Members of parent classes with non-dependent names
      using Iterator<T>::setClassName;
      using Iterator<T>::iterMed;
      using ParamComposite::read;
      using ParamComposite::readOptional;

   };

   template <typename T>
   inline AmStrategy<T>& AmIterator<T>::strategy() 
   {  return *strategy_; }

}
#endif

#ifndef PSPC_AM_ITERATOR_H
#define PSPC_AM_ITERATOR_H

/*
* PSCF - Polymer Self-Consistent Field Theory
*
* Copyright 2016 - 2019, The Regents of the University of Minnesota
* Distributed under the terms of the GNU General Public License.
*/

#include <pspc/iterator/Iterator.h> // base class
#include <pspc/solvers/Mixture.h>
#include <pscf/math/LuSolver.h>
#include <util/containers/DArray.h>
#include <util/containers/FArray.h>
#include <util/containers/FSArray.h>
#include <util/containers/DMatrix.h>
#include <util/containers/RingBuffer.h>
#include <pspc/field/RField.h>

namespace Pscf {
namespace Pspc
{

   using namespace Util;
   
   /**
   * Anderson mixing iterator for the pseudo spectral method.
   *
   * \ingroup Pspc_Iterator_Module
   */
   template <int D>
   class AmIterator : public Iterator<D>
   {
   public:
      
      typedef RField<D> WField;
      typedef RField<D> CField;

      /**
      * Constructor
      *
      * \param system pointer to a parent System object
      */
      AmIterator(System<D>& system);

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
      * Should the unit cell be adjusted to find a zero stress state?
      */
      bool isFlexible();

      /**
      * Check if ensemble is canonical. Returns false if grand-canonical 
      * or mixed ensemble.
      */
      bool isCanonical();

      /**
      * Get epsilon (error threshhold).
      */
      double epsilon();

      /**
      * Get the maximum number of field histories retained.
      */
      int maxHist();

      /**
      * Get the maximum number of iteration before convergence.
      */
      int maxItr();

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
      void minimizeCoeff(int itr);

      /**
      * Rebuild wFields for the next iteration from minimized coefficients
      */
      void buildOmega(int itr);

      /**
      * Clean up after a call to solve(), enabling future calls to solve.
      */
      void cleanUp();

   private:

      /// Error tolerance
      double epsilon_;

      /// Flexible cell computation (1) or rigid (0), default value = 0
      bool isFlexible_;

      /// Free parameter for minimization
      double lambda_;

      /// Number of previous steps to use to compute next state. [0,maxHist_]
      int nHist_;

      // Number of histories to retain.
      int maxHist_;

      /// Maximum number of iterations to attempt.
      int maxItr_;

      // Ensemble shift factor. 1 if canonical, 0 if mixed or grand-canonical
      // Affects the number of basis functions treated by the residual. Ignore
      // the spatially homogeneous basis function if canonical.
      int shift_;  

      // Work Array for iterating on parameters 
      FSArray<double, 6> parameters_;

      /// holds histories of deviation for each monomer
      /// 1st index = history, 2nd index = monomer, 3rd index = ngrid
      // The ringbuffer used is now slightly modified to return by reference
      RingBuffer< DArray < DArray<double> > > resHists_;

      RingBuffer< DArray < DArray<double> > > wHists_;

      /// History of deviation for each cell parameter
      /// 1st index = history, 2nd index = cell parameter
      // The ringbuffer used is now slightly modified to return by reference
      RingBuffer< FArray <double, 6> > stressHists_;

      RingBuffer< FSArray<double, 6> > cellParamHists_;

      /// Matrix, the dot products of residual differences.
      DMatrix<double> U_;

      /// Cn, coefficients for mixing prevous histories
      DArray<double> coeffs_;

      /// Vector, dot products of residuals with differences from histories
      DArray<double> v_;

      /// bigW, blended omega fields
      DArray<DArray <double> > wArrays_;

      /// bigD, blened deviation fields. new wFields = bigW + lambda * bigD
      DArray<DArray <double> > dArrays_;

      /// bigWcP, blended parameter
      FArray <double, 6> wCpArrays_;

      /// bigDCp, blened deviation parameter. new wParameter = bigWCp + lambda * bigDCp
      FArray <double, 6> dCpArrays_;

      // workspace for residual calculation
      DArray< DArray<double> > resArrays_;

      using Iterator<D>::setClassName;
      using Iterator<D>::system;
      using ParamComposite::read;
      using ParamComposite::readOptional;

   //friend:


   };

   template<int D>
   inline double AmIterator<D>::epsilon()
   { return epsilon_; }

   template<int D>
   inline int AmIterator<D>::maxHist()
   { return maxHist_; }

   template<int D>
   inline int AmIterator<D>::maxItr()
   { return maxItr_; }

   #ifndef PSPC_AM_ITERATOR_TPP
   // Suppress implicit instantiation
   extern template class AmIterator<1>;
   extern template class AmIterator<2>;
   extern template class AmIterator<3>;
   #endif

}
}
#endif

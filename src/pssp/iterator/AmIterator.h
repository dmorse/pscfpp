#ifndef PSSP_AM_ITERATOR_H
#define PSSP_AM_ITERATOR_H

/*
* PSCF - Polymer Self-Consistent Field Theory
*
* Copyright 2016, The Regents of the University of Minnesota
* Distributed under the terms of the GNU General Public License.
*/

#include <pssp/iterator/Iterator.h> // base class
#include <pssp/solvers/Mixture.h>
#include <pscf/math/LuSolver.h>
#include <util/containers/DArray.h>
#include <util/containers/FArray.h>
#include <util/containers/FSArray.h>
#include <util/containers/DMatrix.h>
#include <util/containers/RingBuffer.h>
//#include <pssp/iterator/RingBuffer.h>
#include <pssp/field/RField.h>


namespace Pscf {
namespace Pssp
{

   using namespace Util;
   
   /**
   * Anderson mixing iterator for the pseudo spectral method.
   *
   * \ingroup Pssp_Iterator_Module
   */
   template <int D>
   class AmIterator : public Iterator<D>
   {
   public:
      
      typedef RField<D> WField;
      typedef RField<D> CField;

      #if 0
      /**
      * Default constructor
      */
      AmIterator();
      #endif

      /**
      * Constructor
      *
      * \param system pointer to a system object
      */
      AmIterator(System<D>* system);

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
      * Allocate all arrays
      */
      void allocate();

      /**
      * Iterate to a solution
      */
      int solve();

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
      void computeDeviation();

      /**
      * Check if solution is converge within specified tolerance.
      *
      * \return true for error < epsilon and false for error >= epsilon
      */
      bool isConverged();

      /**
      * Determine the coefficients that would minimize invertMatrix_ Umn
      */
      void minimizeCoeff(int itr);

      /**
      * Rebuild wFields for the next iteration from minimized coefficients
      */
      void buildOmega(int itr);

   private:

      /// Error tolerance
      double epsilon_;

      /// Variable cell computation (1) or fixed cell computation (0), default value = 0
      bool cell_;

      /// Free parameter for minimization
      double lambda_;

      /// Number of previous steps to use to compute next state. [0,maxHist_]
      int nHist_;

      // Number of histories to retain.
      int maxHist_;

      /// Maximum number of iterations to attempt.
      int maxItr_;

      // Work Array for iterating on parameters 
      FSArray<double, 6> parameters;

      /// holds histories of deviation for each monomer
      /// 1st index = history, 2nd index = monomer, 3rd index = ngrid
      // The ringbuffer used is now slightly modified to return by reference
      RingBuffer< DArray < DArray<double> > > devHists_;

      RingBuffer< DArray < DArray<double> > > omHists_;

      /// History of deviation for each cell parameter
      /// 1st index = history, 2nd index = cell parameter
      // The ringbuffer used is now slightly modified to return by reference
      RingBuffer< FArray <double, 6> > devCpHists_;

      RingBuffer< FSArray<double, 6> > CpHists_;

      /// Umn, matrix to be minimized
      DMatrix<double> invertMatrix_;

      /// Cn, coefficient to convolute previous histories with
      DArray<double> coeffs_;

      DArray<double> vM_;

      /// bigW, blended omega fields
      DArray<DArray <double> > wArrays_;

      /// bigD, blened deviation fields. new wFields = bigW + lambda * bigD
      DArray<DArray <double> > dArrays_;

      /// bigWcP, blended parameter
      FArray <double, 6> wCpArrays_;

      /// bigDCp, blened deviation parameter. new wParameter = bigWCp + lambda * bigDCp
      FArray <double, 6> dCpArrays_;

      DArray< DArray<double> > tempDev;

      using Iterator<D>::setClassName;
      using Iterator<D>::systemPtr_;
      using ParamComposite::read;
      using ParamComposite::readOptional;

   //friend:
   //for testing purposes


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

   #ifndef PSSP_AM_ITERATOR_TPP
   // Suppress implicit instantiation
   extern template class AmIterator<1>;
   extern template class AmIterator<2>;
   extern template class AmIterator<3>;
   #endif

}
}
#endif

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

      /**
      * Ensemble shift factor: 1 if isCanonical, 0 otherwise.
      *
      * A mixture ensemble is canonical iff all polymer and solvent species
      * have closed ensembles. If ensemble isCanonical, ignore coefficients
      * of the n=0 (spatially homogeneous) basis function for some purposes.
      */
      int shift_;  

      /// Work Array for iterating on parameters 
      FSArray<double, 6> parameters_;

      /**
      * History of field residuals.
      * 1st index = history, 2nd index = monomer, 3rd index = basis func.
      */
      RingBuffer< DArray < DArray<double> > > resHists_;

      /// History of previous w-fields
      RingBuffer< DArray < DArray<double> > > wHists_;

      /// History of previous stress values.
      RingBuffer< FArray <double, 6> > stressHists_;

      /// History of unit cell parameter values.
      RingBuffer< FSArray<double, 6> > cellParamHists_;

      /// Matrix, the dot products of residual differences.
      DMatrix<double> U_;

      /// Cn, coefficients for mixing previous states.
      DArray<double> coeffs_;

      /// Vector, dot products of residuals with differences from histories
      DArray<double> v_;

      /// New trial w field (big W in Arora et al. 2017)
      DArray<DArray <double> > wArrays_;

      /// Predicted field residual for trial state (big D)
      DArray<DArray <double> > dArrays_;

      /// New trial vector of cell parameters.
      FArray<double, 6> wCpArrays_;

      /// Predicted stress residual.
      FArray<double, 6> dCpArrays_;

      /// Workspace for residual calculation.
      DArray< DArray<double> > resArrays_;

      /**
      * Check if ensemble is canonical. Returns false if grand-canonical 
      * or mixed ensemble.
      */
      bool isCanonical();

      /**
      * Return the number of components for a given residual. This is 
      * either the number of spectral basis functions if the residual
      * is an SCF residual, or 1 if the residual is a stress residual. 
      */ 
      int nElem(int i);

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

      // Members of parent classes with non-dependent names
      using Iterator<D>::setClassName;
      using Iterator<D>::system;
      using ParamComposite::read;
      using ParamComposite::readOptional;

   //friend:

   };

   // Inline member functions

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

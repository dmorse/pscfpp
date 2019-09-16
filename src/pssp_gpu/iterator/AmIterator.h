#ifndef PSSP_GPU_AM_ITERATOR_H
#define PSSP_GPU_AM_ITERATOR_H

/*
* PSCF - Polymer Self-Consistent Field Theory
*
* Copyright 2016, The Regents of the University of Minnesota
* Distributed under the terms of the GNU General Public License.
*/

#include <pssp_gpu/iterator/Iterator.h> // base class
#include <pssp_gpu/solvers/Mixture.h>
#include <pscf/math/LuSolver.h>
#include <util/containers/DArray.h>
#include <util/containers/DMatrix.h>
#include <pssp_gpu/iterator/HistMat.h>
#include <pssp_gpu/iterator/RingBuffer.h>
#include <pssp_gpu/field/RDField.h>

//determines if Anderson mixing is done with gpu or cpu
#define GPU_OUTER

namespace Pscf {
namespace Pssp_gpu
{
   using namespace Util;
   
   /**
   * Anderson mixing iterator for the pseudo spectral method
   *
   * \ingroup Pssp_gpu_Iterator_Module
   */
   template <int D>
   class AmIterator : public Iterator<D>
   {
   public:
      
      typedef RDField<D> WField;
      typedef RDField<D> CField;
      /**
      * Default constructor
      */
      AmIterator();

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
      *
      */
      void allocate();

      /**
      * Iterate to a solution
      *
      */
      int solve();

      /**
      * Getter for epsilon
      */
      double epsilon();

      /**
      * Getter for the maximum number of field histories to 
      * convolute into a new field
      */
      int maxHist();

      /**
      * Getter for the maximum number of iteration before convergence
      */
      int maxItr();

      /**
      * Compute the deviation of wFields from a mean field solution
      */
      void computeDeviation();

      /**
      * Compute the error from deviations of wFields and compare with epsilon_
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

      ///error tolerance
      double epsilon_;

      ///free parameter for minimization
      double lambda_;

      ///number of histories to convolute into a new solution [0,maxHist_]
      int nHist_;

      //maximum number of histories to convolute into a new solution
      //AKA size of matrix
      int maxHist_;

      ///number of maximum iteration to perform
      int maxItr_;

      /// holds histories of deviation for each monomer
      /// 1st index = history, 2nd index = monomer, 3rd index = ngrid
      // The ringbuffer used is now slightly modified to return by reference
      #ifndef GPU_OUTER
      RingBuffer< DArray < DArray<double> > > devHists_;
      RingBuffer< DArray < DArray<double> > > omHists_;
      #else
      RingBuffer< DArray < RDField<D> > > devHists_;
      RingBuffer< DArray < RDField<D> > > omHists_;
      #endif


      /// Umn, matrix to be minimized
      DMatrix<double> invertMatrix_;

      /// Cn, coefficient to convolute previous histories with
      DArray<double> coeffs_;

      DArray<double> vM_;

      #ifndef GPU_OUTER
      /// bigW, blended omega fields
      DArray<DArray <double> > wArrays_;

      /// bigD, blened deviation fields. new wFields = bigW + lambda * bigD
      DArray<DArray <double> > dArrays_;

      DArray< DArray<double> > tempDev;
      #else
      /// bigW, blended omega fields
      DArray<RDField<D> > wArrays_;

      /// bigD, blened deviation fields. new wFields = bigW + lambda * bigD
      DArray<RDField<D> > dArrays_;

      DArray<RDField<D> > tempDev;
      #endif

#ifndef GPU_OUTER
	  HistMat <double> histMat_;
#else
	  HistMat <cufftReal> histMat_;
#endif

      using Iterator<D>::setClassName;
      using Iterator<D>::systemPtr_;
      using ParamComposite::read;

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

}
}
#include "AmIterator.tpp"
#endif
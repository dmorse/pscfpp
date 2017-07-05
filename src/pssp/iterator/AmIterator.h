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
#include <util/containers/DMatrix.h>
#include <util/containers/RingBuffer.h>


namespace Pscf {
namespace Pssp
{
	using namespace Util;
   
   /**
   * Anderson mixing iterator for the pseudo spectral method
   *
   * \ingroup Pssp_Iterator_Module
   */
   template <int D>
	class AmIterator : public Iterator<D>
	{
   public:
      
      typedef typename Iterator<D>::WField WField;
      typedef typename Iterator<D>::CField CField;
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

      int maxHist();

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
      int maxHist_;

   	///number of maximum iteration to perform
   	int maxItr_;

   	/// holds histories of deviation for each monomer
   	/// 1st index = history, 2nd index = monomer, 3rd index = ngrid
   	// data member must be integral type or pointer otherwise incur large
   	// data transfer penalty. Homebrewed RingPBuffer might provide better speed
   	RingBuffer< DArray < DArray<double> > > devHists_;

      RingBuffer< DArray < DArray<double> > > omHists_;

   	/// Umn, matrix to be minimized
   	DMatrix<double> invertMatrix_;//prefer a GMatrix here. Not available

   	/// Cn, coefficient to convolute previous histories with
   	DArray<double> coeffs_;

      DArray<double> vM_;

   	/// bigW, blended omega fields
   	DArray<DArray <double> > wArrays_;

   	/// bigD, blened deviation fields. new wFields = bigW + lambda * bigD
   	DArray<DArray <double> > dArrays_;

      using Iterator<D>::setClassName;
      using Iterator<D>::systemPtr_;
      using ParamComposite::read;


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
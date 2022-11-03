#ifndef PSPC_MC_SIMULATOR_H
#define PSPC_MC_SIMULATOR_H

/*
* Simpatico - Simulation Package for Polymeric and Molecular Liquids
*
* Copyright 2010 - 2017, The Regents of the University of Minnesota
* Distributed under the terms of the GNU General Public License.
*/

#include "McMove.h"                      // base class template parameter
#include "McState.h"                     // member
#include <util/param/Manager.h>          // base class template
#include <util/containers/DArray.h>      // member template

namespace Util { class Random; }

namespace Pscf {
namespace Pspc {

   using namespace Util;

   template <int D> class System;

   /**
   * Manager for a set of McMove objects.
   *
   * \sa \ref mcMd_mcMove_McSimulator_page "parameter file format"
   *
   * \ingroup Pspc_Manager_Module
   * \ingroup Pspc_McMove_Module
   */
   template <int D>
   class McSimulator : public Manager< McMove<D> >
   {

   public:

      /**
      * Constructor.
      *
      * \param system parent System
      */
      McSimulator(System<D>& system);

      /**
      * Destructor.
      */
      ~McSimulator();

      /**
      * Read instructions for creating McMove objects.
      *
      * \param in input parameter stream
      */
      virtual void readParameters(std::istream &in);

      /**
      * Perform a field theoretic Monte-Carlo simulation.
      *
      * Perform a field theoretic Monte-Carlo simulation using the 
      * partial saddle-point approximation. 
      * 
      * \param nStep  number of Monte-Carlo steps
      */
      void simulate(int nStep);

      /**
      * Compute the Hamiltonian used in Monte-Carlo simulations.
      */
      void computeMcHamiltonian();

      /**
      * Get the Hamiltonian used in Monte-Carlo simulations.
      */
      double mcHamiltonian() const;

      /**
      * Initialize at beginning of system run.
      *
      * This method calls the initialize method for every McMove.
      */
      void setup();

      /**
      * Choose an McMove at random, using specified probabilities.
      *
      * \return chosen McMove
      */
      McMove<D>& chooseMove();

      /**
      * Save a copy of the Monte-Carlo state.
      *
      * This function and restoreMcState() are intended for use in 
      * the implementation of field theoretic Monte Carlo moves. This
      * function stores the current w fields and the corresponding 
      * Hamiltonian value.  This is normally the first step of a MC 
      * move, prior to an attempted modification of the fields stored
      * in the system w field container.
      */
      void saveMcState();

      /**
      * Restore the saved copy of the Monte-Carlo state.
      *
      * This function  and saveMcState() are intended to be used 
      * together in the implementation of Monte-Carlo moves. If an 
      * attempted move is rejected, restoreMcState() is called to
      * restore the fields ahd Hamiltonian value that were saved 
      * by a previous call to the function saveMcState().
      */
      void restoreMcState();

      /**
      * Output statistics for all moves.
      */
      void output();

      /**
      * Return probability of move i.
      *
      * \param i index for McMove
      * \return probability of McMove number i
      */
      double probability(int i) const;

      double chiPvalue(int i) const
      {  return chiPvalues_[i]; }

      DArray<double> const & chiPvector(int i) const
      {  return chiPvectors_[i]; }

      /**
      * Get parent system by reference.
      */
      System<D>& system();

      /**
      * Get random number generator by reference.
      */
      Random& random();

      using Manager< McMove<D> >::size;

   protected:

      using Manager< McMove<D> >::setClassName;
       
   private:

      // Private data members

      /**
      * Array of McMove probabilities.
      */
      DArray<double>  probabilities_;

      /**
      * Projected chi matrix
      *
      * Projected matrix chiP_ = P*chi*P, where P is projection matrix that 
      * projects onto the subspace orthogonal to the vector e = [1, ... , 1].
      */
      DMatrix<double> chiP_;

      /**
      * Eigenvalues of the projected chi matrix.
      */
      DArray<double>  chiPvalues_;

      /**
      * Eigenvectors of the projected chi matrix.
      */
      DArray< DArray<double> >  chiPvectors_;

      /**
      * State saved during MC simulation.
      */
      mutable McState<D> mcState_;

      /**
      * Monte-Carlo System Hamiltonian (extensive value).
      */
      double mcHamiltonian_;

      /**
      * Has the MC Hamiltonian been computed for the current w and c fields?
      */ 
      bool hasMcHamiltonian_;

      /**
      * Pointer to parent System.
      */
      System<D>* systemPtr_;

      /**
      * Pointer to random number generator.
      */
      Random* randomPtr_;

      // Private member functions

      /**
      * Return pointer to a new McMoveFactory.
      */
      virtual Factory< McMove<D> >* newDefaultFactory() const;

      /**
      * Perform eigenvalue analysis of projected chi matrix.
      */
      void analyzeChi();

   };

   // Inline functions

   /*
   * Return probability of move number i
   */
   template <int D>
   inline double McSimulator<D>::probability(int i) const
   {
      assert(i >= 0);  
      assert(i < size());  
      return probabilities_[i];
   }

   // Get the precomputed MC Hamiltonian
   template <int D>
   inline double McSimulator<D>::mcHamiltonian() const
   {  
      UTIL_CHECK(hasMcHamiltonian_);
      return mcHamiltonian_; 
   }

   // Get the precomputed MC Hamiltonian
   template <int D>
   inline System<D>& McSimulator<D>::system()
   {  
      UTIL_CHECK(systemPtr_);
      return *systemPtr_; 
   }

   // Get the precomputed MC Hamiltonian
   template <int D>
   inline Random& McSimulator<D>::random()
   {  
      UTIL_CHECK(randomPtr_);
      return *randomPtr_; 
   }

}
}
#endif

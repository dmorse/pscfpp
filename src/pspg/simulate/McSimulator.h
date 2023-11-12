#ifndef PSPG_MC_SIMULATOR_H
#define PSPG_MC_SIMULATOR_H

/*
* PSCF - Polymer Self-Consistent Field Theory
*
* Copyright 2016 - 2022, The Regents of the University of Minnesota
* Distributed under the terms of the GNU General Public License.
*/

#include <util/param/ParamComposite.h>       // base class

#include <prdc/cuda/RField.h>                // member (template arg.)
#include <util/random/Random.h>              // member
#include <util/containers/DArray.h>          // member (template)
#include <util/containers/DMatrix.h>         // member (template)

#include "McState.h"                                 // member
#include <util/param/Manager.h>                      // member template
#include <pspg/simulate/mcmove/McMoveManager.h>      // member
#include <pspg/simulate/analyzer/AnalyzerManager.h>  // member

namespace Pscf {
namespace Pspg {

   template <int D> class System;

   template <int D> class McMove;
   template <int D> class TrajectoryReader;
   template <int D> class TrajectoryReaderFactory;

   using namespace Util;
   using namespace Pscf::Prdc;
   using namespace Pscf::Prdc::Cuda;

   /**
   * Field theoretic Monte-Carlo simulator.
   */
   template <int D>
   class McSimulator : public ParamComposite
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
      * Read parameters for a MC simulation.
      *
      * \param in input parameter stream
      */
      virtual void readParameters(std::istream &in);

      /// \name Primary Actions: Simulation and Analysis
      ///@{

      /**
      * Perform a field theoretic simulation.
      *
      * Perform a field theoretic simulation using the partial
      * saddle-point approximation.
      *
      * \param nStep  number of Monte-Carlo steps
      */
      void simulate(int nStep);

      /**
      * Read and analyze a trajectory file.
      *
      * This function uses an instance of the TrajectoryReader class
      * specified by the "classname" argument to read a trajectory
      * file.
      *
      * \param min  start at this frame number
      * \param max  end at this frame number
      * \param classname  name of the TrajectoryReader class to use
      * \param filename  name of the trajectory file
      */
      void analyzeTrajectory(int min, int max,
                             std::string classname,
                             std::string filename);

      /**
      * Output timing results.
      *
      * \param out  output stream
      */
      void outputTimers(std::ostream& out);

      /**
      * Clear all timers.
      */
      void clearTimers();

      /**
      * Get the current simulation step index
      */
      long iStep();

      ///@}
      /// \name Hamiltonian Computation
      ///@{

      /**
      * Compute the field theoretic Hamiltonian and its components.
      */
      void computeHamiltonian();

      /**
      * Get the pre-computed field theoretic Hamiltonian.
      */
      double hamiltonian() const;

      /**
      * Get the ideal gas Hamiltonian contribution.
      */
      double idealHamiltonian() const;

      /**
      * Get the quadratic field Hamiltonian contribution (H_W).
      */
      double fieldHamiltonian() const;

      /**
      * Has the Hamiltonian been computed for the current w and c fields?
      */
      bool hasHamiltonian() const;

      ///@}
      /// \name Projected Chi matrix
      ///@{

      /**
      * Perform eigenvalue analysis of projected chi matrix.
      */
      void analyzeChi();

      /**
      * Get an array of the eigenvalues of the projected chi matrix.
      *
      * The projected chi matrix is given by the matrix product P*chi*P,
      * where P is the symmetric projection matrix that projects onto the
      * subspace orthogonal to the vector e = (1,1,...,1). The projected
      * chi matrix is singular, and always has one zero eigenvalue, with
      * associated eigenvector e. By convention, this zero eigenvalue and
      * eigenvector e are listed last, with index nMonomer -  1.
      */
      DArray<double> const & chiEvals() const
      {  return chiEvals_; }

      /**
      * Get an single eigenvalue of the projected chi matrix.
      *
      * \param i index of eigenvalue (0, ... , nMonomer - 1)
      */
      double chiEval(int i ) const
      {  return chiEvals_[i]; }

      /**
      * Get a matrix of eigenvectors of the projected chi matrix.
      *
      * The first index of the matrix indexes the eigenvector, while
      * the second index indexes vector components. All eigenvectors
      * have a Euclidean norm equal to nMonomr. The sign of each vector
      * is chosen so as to make the first (0) component of each vector
      * positive.
      */
      DMatrix<double> const & chiEvecs() const
      {  return chiEvecs_; }

      ///@}
      /// \name WField Components
      ///@{

      /**
      *  Compute components of W fields in basis of chiP eigenvectors.
      *
      *  Compute and store the components of the values of the w fields
      *  on nodes of a real-space grid (r-grid) in a basis of the
      *  eigenvectors of the projected chi matrix.
      */
      void computeWC();

      /**
      * Get an eigenvector component of the w fields.
      *
      * Each component is a point-wise projection of the w fields onto a
      * corresponding eigenvector of the projected chi matrix. The last
      * index, i = nMonomer - 1, corresponds to the Lagrange multiplier
      * pressure component, with associated eigenvector e = (1,1,...,1).
      *
      * \param i index for eigenvector / eigenvalue pair
      */
      RField<D> const & wc(int i) const
      {   return wc_[i]; }

      /**
      * Are eigen-components of the current w fields valid?
      */
      bool hasWC() const;

      /**
      * Clear w field eigen-components and Hamiltonian components.
      */
      void clearData();

      ///@}
      /// \name Utilities for MC simulation
      ///@{

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
      * Clear the saved copy of the Monte-Carlo state.
      *
      * This function, restoreMcState(), and saveMcState() are intended to be used
      * together in the implementation of Monte-Carlo moves. If an
      * attempted move is accepted, clearMcState() is called to
      * clear mcState_.hasData
      */
      void clearMcState();

      ///@}
      /// \name Miscellaneous
      ///@{

      /**
      * Get parent system by reference.
      */
      System<D>& system();

      /**
      * Get random number generator by reference.
      */
      Random& random();

      /**
      * Get AnalyzerManger
      */
      AnalyzerManager<D>& analyzerManager();

      /**
      * Get McMoveManger
      */
      McMoveManager<D>& mcMoveManager();

      /**
      * Get the trajectory reader factory by reference.
      */
      Factory<TrajectoryReader<D>>& trajectoryReaderFactory();

   protected:

      using Util::ParamComposite::setClassName;

      /**
      * Random number generator
      */
      Random random_;

      /**
      * Eigenvector components of w on a real space grid.
      *
      * Each field component corresponds to a point-wise projection of w
      * onto an eigenvector of the projected chi matrix.
      */
      DArray< RField<D> > wc_;

      /**
      * Monte-Carlo System Hamiltonian (extensive value).
      */
      double hamiltonian_;

      /**
      * Ideal gas contributions (lnQ) to Monte-Carlo System Hamiltonian
      */
      double idealHamiltonian_;

      /**
      * Field contribution (HW) to Monte-Carlo System Hamiltonian
      */
      double fieldHamiltonian_;

      /**
      * Has the MC Hamiltonian been computed for the current w and c fields?
      */
      bool hasHamiltonian_;

      /**
      * Has the eigenvector components of the current w fields been computed
      * for the current field?
      */
      bool hasWC_;

      /**
      * Count Monte Carlo step
      */
      long iStep_;

   private:

      /**
      * Projected chi matrix
      *
      * Projected matrix chiP_ = P*chi*P, where P is projection matrix that
      * projects onto the subspace orthogonal to the vector e = [1, ... , 1].
      */
      DMatrix<double> chiP_;

      /**
      * Eigenvectors of the projected chi matrix.
      */
      DMatrix<double> chiEvecs_;

      /**
      * Eigenvalues of the projected chi matrix.
      */
      DArray<double>  chiEvals_;

      /**
      * Pointer to the parent system.
      */
      System<D>* systemPtr_;

      /**
      * State saved during MC simulation.
      */
      mutable McState<D> mcState_;

      /**
      * Manger for Monte Carlo Move.
      */
      McMoveManager<D> mcMoveManager_;

      /**
      * Manger for Monte Carlo Analyzer.
      */
      AnalyzerManager<D> analyzerManager_;

      /**
      * Pointer to a trajectory reader/writer factory.
      */
      Factory<TrajectoryReader<D>>* trajectoryReaderFactoryPtr_;

      // Private member functions

      /**
      * Called at the beginning of the simulation member function.
      */
      void setup();

   };

   // Inline functions

   // Get the Monte-Carlo move manager.
   template <int D>
   inline McMoveManager<D>& McSimulator<D>::mcMoveManager()
   {  return mcMoveManager_; }

   // Get the Monte-Carlo analyzer manager.
   template <int D>
   inline AnalyzerManager<D>& McSimulator<D>::analyzerManager()
   {  return analyzerManager_; }

   // Get the random number generator.
   template <int D>
   inline Random& McSimulator<D>::random()
   {  return random_; }

   // Get the parent System.
   template <int D>
   inline System<D>& McSimulator<D>::system()
   {  return *systemPtr_; }

   // Get the precomputed MC Hamiltonian
   template <int D>
   inline double McSimulator<D>::hamiltonian() const
   {
      UTIL_CHECK(hasHamiltonian_);
      return hamiltonian_;
   }

   // Get the ideal gas component of the precomputed MC Hamiltonian
   template <int D>
   inline double McSimulator<D>::idealHamiltonian() const
   {
      UTIL_CHECK(hasHamiltonian_);
      return idealHamiltonian_;
   }

   // Get the field component of the precomputed MC Hamiltonian
   template <int D>
   inline double McSimulator<D>::fieldHamiltonian() const
   {
      UTIL_CHECK(hasHamiltonian_);
      return fieldHamiltonian_;
   }

   // Has the MC Hamiltonian been computed for the current w fields ?
   template <int D>
   inline bool McSimulator<D>::hasHamiltonian() const
   {  return hasHamiltonian_; }

   // Have eigen-components of current w fields been computed?
   template <int D>
   inline bool McSimulator<D>::hasWC() const
   {  return hasWC_; }

   // Clear all data (eigen-components of w field and Hamiltonian)
   template <int D>
   inline void McSimulator<D>::clearData()
   {
      hasHamiltonian_ = false;
      hasWC_ = false;
   }

   // Get the TrajectoryReaderfactory
   template <int D>
   inline 
   Factory<TrajectoryReader<D>> & McSimulator<D>::trajectoryReaderFactory()
   {
      UTIL_ASSERT(trajectoryReaderFactoryPtr_);
      return *trajectoryReaderFactoryPtr_;
   }

   template <int D>
   inline long McSimulator<D>::iStep()
   {
      return iStep_;
   }

   #ifndef PSPG_MC_SIMULATOR_TPP
   // Suppress implicit instantiation
   extern template class McSimulator<1>;
   extern template class McSimulator<2>;
   extern template class McSimulator<3>;
   #endif

}
}
#endif

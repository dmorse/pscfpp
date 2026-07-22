#ifndef CPC_SIMULATOR_H
#define CPC_SIMULATOR_H

/*
* PSCF - Polymer Self-Consistent Field
*
* Copyright 2015 - 2025, The Regents of the University of Minnesota
* Distributed under the terms of the GNU General Public License.
*/

#include <util/param/ParamComposite.h>     // base class
#include <prdc/field/cpu/CField.h>               // member (template arg)
#include <util/random/Random.h>            // member
#include <util/containers/DArray.h>        // member (template)
#include <util/containers/DMatrix.h>       // member (template)


namespace Pscf {
namespace Cpc {

   // Forward declarations 
   template <int D> class System;
   template <int D> class Step;
   // template <int D> class TrajectoryReader;

   using namespace Util;
   using namespace Prdc;

   /**
   * Simulator for complex Langevin field theoretic simulation.
   *
   * The Simulator class provides functions to compute and diagonalze a 
   * pairwise interction matrix, functions to access components of several 
   * types of fields in a basis of eigenvectors of the interaction matrix,
   * and functions to compute and return contributions to the field theoretic 
   * Hamiltonian.
   *
   * The analyzeU function constructs and diagonalizes the interaction 
   * matrix U.  The computeWc, computeCc and computeDc functions compute 
   * components of various types of multi-component fields (i.e., fields 
   * that are associated with a monomer type index) in a basis of eigenvectors 
   * of the M matrix. Names such as wc, cc and dc that end with a suffix "c" 
   * refer to components of multi-component fields that are defined using 
   * this eigenvector basis. 
   *
   * \ingroup Cpc_Fts_Module
   */
   template <int D>
   class Simulator : public ParamComposite
   {

   public:

      // Suppress automatically generated functions
      Simulator() = delete;
      Simulator(Simulator<D>& ) = delete;

      /**
      * Constructor.
      *
      * \param system parent System
      */
      Simulator(System<D>& system);

      /**
      * Destructor.
      */
      ~Simulator();

      /**
      * Allocate required memory.
      *
      * Values of nMonomer and the mesh dimensions must be defined in
      * Mixture and Domain members of the parent System on entry. This
      * function must be called by the readParameters method of any
      * subclass.
      */
      void allocate();

      /**
      * Read parameters for a simulation.
      *
      * \param in input parameter stream
      */
      virtual void readParameters(std::istream &in);

      /// \name Primary Actions: Simulation and Analysis
      ///@{

      /**
      * Perform a complex Langevin field theoretic simulation.
      *
      * The default implemention is a do-nothing placeholder that throws
      * an error if called, and must be re-implemented by subclasses. A
      * do-nothing placeholder is provided to make this a concrete class
      * simplify unit testing.
      *
      * \param nStep  number of simulation steps
      */
      virtual void simulate(int nStep);

      #if 0
      /**
      * Read and analyze a trajectory file.
      *
      * This function uses an instance of the TrajectoryReader class
      * specified by the "classname" argument to read a trajectory file
      * with the specified filename. The function opens the file,
      * performs the analysis, and closes the file before returning.
      *
      * The default implemention is a do-nothing placeholder that throws
      * an error if called, and must be re-implemented by subclasses.
      *
      * \param min  first frame number
      * \param max  last frame number
      * \param classname  name of TrajectoryReader class
      * \param filename  name of trajectory file
      */
      virtual void analyze(int min, int max,
                           std::string classname,
                           std::string filename);
      #endif

      /**
      * Clear field eigen-components and hamiltonian components.
      *
      * On return from this function, hasHamiltonian(), hasWc(), hasCc(),
      * and hasDc() all return false.
      */
      void clearData();

      ///@}
      /// \name Timers and Counters
      ///@{

      /**
      * Output timing results
      *
      * Empty default implementation.
      *
      * \param out  output stream
      */
      virtual void outputTimers(std::ostream& out) const;

      /**
      * Clear timers.
      *
      * Empty default implementation.
      */
      virtual void clearTimers();

      /**
      * Return the current converged simulation step index.
      */
      long iStep();
      
      /**
      * Return the current simulation step index.
      */
      long iTotalStep();

      ///@}
      /// \name Interaction Matrix
      ///@{

      /**
      * Perform eigenvalue analysis of the interaction matrix U.
      *
      * Uses data obtained from the Interaction member of the parent 
      * System.
      */
      void analyzeInteraction();

      /**
      * Get the interaction matrix U by const reference.
      */
      DMatrix<double> const & u() const;

      /**
      * Get an array of the eigenvalues of the interaction matrix U.
      */
      DArray<double> const & evals() const;

      /**
      * Get a single eigenvalue of the interaction matrix.
      *
      * \param i  index of eigenvalue (0, ... , nMonomer - 1)
      */
      double eval(int i) const;

      /**
      * Is the associated eigenvalue of U positive?
      *
      * \param i  index of eigenvalue (0, ... , nMonomer - 1)
      */
      bool isPositiveEval(int i) const;

      /**
      * Get the matrix of all eigenvectors of the interaction matrix.
      *
      * This function returns the entire nMonomer x nMonomer matrix of the
      * eigenvectors of the interaction matrix, in which each row is an
      * eigenvector. The first (row) index of this matrix thus identifies
      * an eigenvector, while the second (column) index identifies the
      * monomer type associated with one component of an eigen-vector. 
      *
      * Each eigenvector is normalized such that the sum of the squares 
      * of its elements is equal to nMonomer, the number of monomer types. 
      * The sign of each vector is chosen so as to make the first (0) 
      * component non-negative.  The last eigenvector is always the null 
      * vector e = (1,1,...,1).
      * 
      * For the case nMonomer = 2 of an AB system, the resulting two 
      * eigenvectors are (1,-1) and (1,1).
      */
      DMatrix<double> const & evecs() const;

      /**
      * Get one element of an eigenvector of the interaction matrix.
      *
      * See documentation of evecs(), which returns the entire matrix.
      *
      * \param a  eigenvector index (0, ..., nMonomer - 1)
      * \param i  monomoner type index (0, ..., nMonomer - 1)
      */
      double evecs(int a, int i) const;

      #if 0
      /**
      * Get all components of the vector S.
      *
      * The value of component \f$ S_{a} \f$ may be expressed using 
      * Einstein summation convention as
      * \f[
      *     S_{a} \equiv \frac{1}{M^2} v_{ai}\chi_{ij}e_{j}
      * \f]
      * for any \f$ a = 0, \ldots, M - 1 \f$, where M = nMonomer (the
      * number of monomer types), \f$ e_{j} =1 \f$ for any j, and 
      * \f$ v_{ai} \f$ is component associated with monomer type i of
      * eigenvector a of the interaction matrix, with the convention
      * \f$ v_{ia} = e_{i} = 1 \f$ for a = nMonomer - 1.
      */
      DArray<double> const & sc() const;

      /**
      * Get a single component of the S vector.
      *
      * This function retrieves on component of the vector defined in
      * the documentation for function sc().
      * 
      * \param a  eigenvector index (0, ..., nMonomer - 1)
      */
      double sc(int a) const;
      #endif

      ///@}
      /// \name Field Theoretic Hamiltonian 
      ///@{

      /**
      * Compute the Hamiltonian used in PS-FTS.
      */
      void computeHamiltonian();

      /**
      * Get the Hamiltonian used in PS-FTS.
      *
      * This function returns the real, thermodynamically extensive
      * Hamiltonian used in simulations based on partial saddle-point
      * approximation (PS-FTS).
      */
      std::complex<double> hamiltonian() const;

      /**
      * Get ideal gas contribution to the Hamiltonian.
      */
      std::complex<double> idealHamiltonian() const;

      /**
      * Get the quadratic field contribution to the Hamiltonian.
      */
      std::complex<double> fieldHamiltonian() const;

      /**
      * Has the Hamiltonian been computed for current w and c fields?
      */
      bool hasHamiltonian() const;

      ///@}
      /// \name Chemical Potential Field (W Field) Components
      ///@{

      /**
      * Compute eigenvector components of the current w fields.
      *
      * Compute and store the components of the values of the w fields
      * on nodes of a real-space grid (r-grid) in a basis of the
      * eigenvectors of the interaction matrix. The component field
      * \f$ W_{a}({\bf r}) \f$ at grid point \f$ {\bf r} \f$ is given
      * using Einstein summation by
      * \f[
      *    W_{a}({\bf r}) = 
      *    v_{ai} w_{i}({\bf r}) / M
      * \f]
      * where \f$ w_{i}({\bf r}) \f$ is the w-field associated with
      * monomer type \f$ i \f$, \f$ v_{ai} \f$ is eigenvector a of
      * the interaction matrix, and M = nMonomer.
      */
      void computeWc();

      /**
      * Get all eigenvector components of the current w fields.
      *
      * This function returns a DArray of fields in which each field is
      * a chemical field component \f$ W_{a}({\bf r}) \f$ as defined in
      * the documentation of computeWc(), for a = 0, ..., nMonomer - 1.
      */
      DArray< CField<D, CppTp<D> > > const & wc() const;

      /**
      * Get one eigenvector component of the current w fields.
      *
      * See documentation of functions computeWc() and wc() for details.
      *
      * \param a eigenvector index in range 0 , ..., nMonomer -1
      */
      CField<D, CppTp<D> > const & wc(int a) const;

      /**
      * Are eigen-components of current w fields valid ?
      */
      bool hasWc() const;

      ///@}
      /// \name Monomer Concentration Field (C-Field) Components
      ///@{

      /**
      * Compute eigenvector components of the current c fields.
      *
      * Compute and store the components of the values of the c fields
      * on nodes of a real-space grid (r-grid) in a basis of the
      * eigenvectors of the interaction matrix. 
      */
      void computeCc();

      /**
      * Get all eigenvector components of the current c fields.
      *
      * Each component \f$C_{a}({\bf r}) \f$ is a point-wise projection 
      * of the monomer c fields onto a corresponding eigenvector of the 
      * interaction matrix. The resulting value \f$ C_{a}({\bf r}) \f$ 
      * for eigen-component a at grid point \f$ {\bf r} \f$ is given 
      * using Einstein notation as
      * \f[
      *    C_{a}({\bf r}) = v_{ai} c_{i}({\bf r}) 
      * \f]
      * where \f$ c_{i}({\bf r}) \f$ is the concentration / volume 
      * fraction field associated with monomer type i.  
      *
      * Note: The above definition \f$ C_{a} \f$ uses a different 
      * prefactor than that used to define the corresponding w-field 
      * component \f$ W_{a} \f$ given in the documentation of the
      * function wc(), without the prefactor of 1/nMonomer. This is 
      * intentional, and is convenient for other aspects of the 
      * underlying theory.
      */
      DArray< CField<D, CppTp<D> > > const & cc() const;

      /**
      * Get one eigenvector component of the current c fields.
      *
      * This returns a reference to a field \f$ C_{a}({\bf r}) \f$
      * as defined in the documentation of function cc().
      *
      * \param a eigenvector / eigenvalue index
      */
      CField<D, CppTp<D> > const & cc(int a) const;

      /**
      * Are eigen-components of current c fields valid ?
      */
      bool hasCc() const;

      ///@}
      /// \name Functional Derivatives of H[W]
      ///@{

      /**
      * Compute functional derivatives of the Hamiltonian.
      *
      * Compute and store the functional derivatives of the field
      * theoretic Hamiltonian with respect to eigenvector components of 
      * the w fields (i.e., with respect to components of wc).
      */
      void computeDc();

      /**
      * Get all of the current d fields.
      *
      * This function returns an array of fields in which element a
      * is the functional derivative of the Hamiltonian H[W] with
      * respect to the field component \f$ W_{a} \f$ that is returned
      * by the function wc(a).
      */
      DArray< CField<D, CppTp<D> > > const & dc() const;

      /**
      * Get one eigenvector component of the current d fields.
      *
      * \param i eigenvector / eigenvalue index
      */
      CField<D, CppTp<D> > const & dc(int i) const;

      /**
      * Are the current d fields valid ?
      */
      bool hasDc() const;
      
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
      * Does this Simulator have a Step object?
      */
      bool hasStep() const;

      /**
      * Get the Step by reference.
      */
      Step<D>& step();

      #if 0
      /**
      * Get the AnalyzerManager by reference.
      */
      AnalyzerManager<D>& analyzerManager();

      /**
      * Get the trajectory reader factory by reference.
      */
      Factory<TrajectoryReader<D>>& trajectoryReaderFactory();
      #endif

      ///@}

   private:

      #if 0
      /**
      * Manager for Analyzer.
      */
      AnalyzerManager<D> analyzerManager_;
      #endif

      /**
      * Random number generator
      */
      Random random_;

      // Field components in eigenvector basis

      /**
      * Eigenvector components of w fields on a real space grid.
      *
      * Each field component corresponds to a point-wise projection of 
      * the monomer w fields onto an eigenvector of the interaction 
      * matrix. The number of components is equal to the number of
      * monomer types, nMonomer. The last component is a pressure-like
      * field.
      */
      DArray< CField<D, CppTp<D> > > wc_;

      /**
      * Eigenvector components of c fields on a real space grid.
      *
      * Each field component corresponds to a point-wise projection of 
      * the monomer c fields onto an eigenvector of the interaction 
      * matrix. The number of components is equal to the number of 
      * monomer types, nMonomer. The last component must satisfy an
      * incompressibility constraint.
      */
      DArray< CField<D, CppTp<D> > > cc_;

      /**
      * Components of d fields on a real space grid.
      *
      * Each field component is the functional derivative of H[W]
      * with respect to one eigenvector w-field component.
      */
      DArray< CField<D, CppTp<D> > > dc_;
     
      // Interaction matrix and eigen-properties

      /**
      * Projected chi matrix
      */
      DMatrix<double> u_;

      /**
      * Eigenvectors of the interaction matrix, U.
      *
      * Each row (identified by first index) is an eigenvector. 
      * The last eigenvector, with index nMonomer - 1, is always the
      * vector e = [1, 1, ...., 1]. Distinct eigenvectors are orthogonal.
      * Eigenvectors normalized such that the sum of the square of the 
      * elements is equal to nMonomer.
      */
      DMatrix<double> evecs_;

      /**
      * Eigenvalues of the interaction matrix.
      */
      DArray<double>  evals_;

      /**
      * Array of bools - true for positive eigenvalues.
      */
      DArray<bool> isPositiveEval_;

      // Hamiltonian and components
 
      /**
      * Total field theoretic Hamiltonian H[W] (extensive value).
      */
      std::complex<double> hamiltonian_;

      /**
      * Ideal gas contribution (-lnQ) to Hamiltonian H[W]
      */
      std::complex<double> idealHamiltonian_;

      /**
      * Field contribution (H_W) to Hamiltonian
      */
      std::complex<double> fieldHamiltonian_;

      // Pointers to associated and child objects

      /**
      * Pointer to the parent system.
      */
      System<D>* systemPtr_;

      /**
      * Pointer to Brownian dynamics step algorithm.
      */
      Step<D>* stepPtr_;

      /**
      * Pointer to a Step factory.
      */
      Factory< Step<D> >* stepFactoryPtr_;

      #if 0
      /**
      * Pointer to a trajectory reader/writer factory.
      */
      Factory< TrajectoryReader<D> >* trajectoryReaderFactoryPtr_;
      #endif

      // Counters and random seed

      /**
      * Step counter - attempted steps for which compressor converges.
      *
      * Steps for which the compressor fails to converge are returned to
      * the previous state so that another random displacement can be
      * chosen. Attempted MC moves for which the compressor converged 
      * but which are then rejected based on a Metropolis criterion are
      * included in iStep_. The difference iTotalStep_ - iStep_ is the
      * number of moves that failed because the compressor failed to 
      * converge.
      */
      long iStep_;
      
      /**
      * Step counter - total number of attempted BD or MC steps.
      */
      long iTotalStep_;

      /**
      * Random number generator seed.
      */
      long seed_;

      // Boolean status indicators

      /**
      * Has the Hamiltonian been computed for the current w and c fields?
      */
      bool hasHamiltonian_;

      /**
      * Have eigen-components of the current w fields been computed ?
      */
      bool hasWc_;

      /**
      * Have eigen-components of the current c fields been computed ?
      */
      bool hasCc_;

      /**
      * Have functional derivatives of H[W] been computed ?
      */
      bool hasDc_;

      /**
      * Has required memory been allocated?
      */
      bool isAllocated_;

      // Private member functions

      /**
      * Optionally read a random number generator seed.
      *
      * \param in input parameter stream
      */
      void readRandomSeed(std::istream& in);

      /**
      * Called at the beginning of the simulation.
      *
      * \param nStep  number of time steps planned for simulation
      */
      void setup(int nStep);

   };

   // Inline functions

   // Management of owned and associated objects

   // Get the parent System by reference.
   template <int D>
   inline System<D>& Simulator<D>::system()
   {
      UTIL_ASSERT(systemPtr_);  
      return *systemPtr_; 
   }

   // Get the random number generator by reference.
   template <int D>
   inline Random& Simulator<D>::random()
   {  return random_; }

   // Get the Brownian dynamics stepper.
   template <int D>
   inline bool Simulator<D>::hasStep() const
   {  return (bool)stepPtr_; }

   // Get the Brownian dynamics stepper.
   template <int D>
   inline Step<D>& Simulator<D>::step()
   {
      UTIL_CHECK(hasStep());  
      return *stepPtr_; 
   }

   #if 0
   // Get the analyzer manager.
   template <int D>
   inline AnalyzerManager<D>& Simulator<D>::analyzerManager()
   {  return analyzerManager_; }

   // Get the TrajectoryReaderfactory
   template <int D>
   inline
   Factory<TrajectoryReader<D> >& Simulator<D>::trajectoryReaderFactory()
   {
      UTIL_ASSERT(trajectoryReaderFactoryPtr_);
      return *trajectoryReaderFactoryPtr_;
   }
   #endif

   // Interaction Matrix

   // Return the interaction matrix U.
   template <int D>
   inline DMatrix<double> const & Simulator<D>::u() const
   {  return u_; }

   // Return an array of eigenvalues of the interaction matrix U.
   template <int D>
   inline DArray<double> const & Simulator<D>::evals() const
   {  return evals_; }

   // Return a single eigenvalue of the interaction matrix U.
   template <int D>
   inline double Simulator<D>::eval(int a) const
   {  return evals_[a]; }

   // Return a matrix of eigenvectors of the interaction matrix U.
   template <int D>
   inline DMatrix<double> const & Simulator<D>::evecs() const
   {  return evecs_; }

   // Return an element of an eigenvector of the interaction matrix U.
   template <int D>
   inline double Simulator<D>::evecs(int a, int i) const
   {  return evecs_(a, i); }

   // Hamiltonian and its derivatives

   // Get the precomputed Hamiltonian.
   template <int D>
   inline std::complex<double> Simulator<D>::hamiltonian() const
   {
      UTIL_CHECK(hasHamiltonian_);
      return hamiltonian_;
   }

   // Get the ideal gas component of the precomputed Hamiltonian.
   template <int D>
   inline std::complex<double> Simulator<D>::idealHamiltonian() const
   {
      UTIL_CHECK(hasHamiltonian_);
      return idealHamiltonian_;
   }

   // Get the W field component of the precomputed Hamiltonian.
   template <int D>
   inline std::complex<double> Simulator<D>::fieldHamiltonian() const
   {
      UTIL_CHECK(hasHamiltonian_);
      return fieldHamiltonian_;
   }

   // Has the Hamiltonian been computed for the current w fields ?
   template <int D>
   inline bool Simulator<D>::hasHamiltonian() const
   {  return hasHamiltonian_; }

   // Fields

   // Return all eigencomponents of the w fields.
   template <int D>
   inline DArray< CField<D, CppTp<D> > > const & Simulator<D>::wc() const
   {  return wc_; }

   // Return a single eigenvector component of the w fields.
   template <int D>
   inline CField<D, CppTp<D> > const & Simulator<D>::wc(int a) const
   {  return wc_[a]; }

   // Have eigenvector components of current w fields been computed?
   template <int D>
   inline bool Simulator<D>::hasWc() const
   {  return hasWc_; }

   // Return all eigenvector components of the current c fields.
   template <int D>
   inline DArray< CField<D, CppTp<D> > > const & Simulator<D>::cc() const
   {  return cc_; }

   // Return a single eigenvector component of the current c fields.
   template <int D>
   inline CField<D, CppTp<D> > const & Simulator<D>::cc(int a) const
   {  return cc_[a]; }

   // Have eigenvector components of current c fields been computed?
   template <int D>
   inline bool Simulator<D>::hasCc() const
   {  return hasCc_; }

   // Return all eigenvector components of the current d fields.
   template <int D>
   inline DArray< CField<D, CppTp<D> > > const & Simulator<D>::dc() const
   {  return dc_; }

   // Return a single eigenvector component of the current d fields.
   template <int D>
   inline CField<D, CppTp<D> > const & Simulator<D>::dc(int a) const
   {  return dc_[a]; }

   // Have eigenvector components of current d fields been computed?
   template <int D>
   inline bool Simulator<D>::hasDc() const
   {  return hasDc_; }

   // Return the current converged simulation step index.
   template <int D>
   inline long Simulator<D>::iStep()
   {  return iStep_; }
   
   // Return the current simulation step index.
   template <int D>
   inline long Simulator<D>::iTotalStep()
   {  return iTotalStep_; }

   // Explicit instantiation declarations
   extern template class Simulator<1>;
   extern template class Simulator<2>;
   extern template class Simulator<3>;

} // namespace Cpc
} // namespace Pscf
#endif

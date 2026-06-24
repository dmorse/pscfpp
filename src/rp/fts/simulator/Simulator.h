#ifndef RP_SIMULATOR_H
#define RP_SIMULATOR_H

/*
* PSCF - Polymer Self-Consistent Field
*
* Copyright 2015 - 2025, The Regents of the University of Minnesota
* Distributed under the terms of the GNU General Public License.
*/

#include <util/param/ParamComposite.h>     // base class
#include <util/containers/DArray.h>        // member (template)
#include <util/containers/DMatrix.h>       // member (template)
#include <iostream>

// Forward declaration
namespace Util { 
   class Random; 
}
namespace Pscf {
   namespace Rp {
      template <int D, class T> class System;
      template <int D, class T> class SimState;
   }
}

namespace Pscf {
namespace Rp {

   using namespace Util;

   /**
   * Base class for field theoretic PS-FTS simulator.
   *
   * <b> Purpose </b>: A Simulator base class provides tools needed in 
   * field-theoretic simulations that are based on a partial saddle-point 
   * approximation. In the Rpc and Rpg program-level namespaces, subclasses
   * designed for field theoretic Monte Carlo (MC) and Brownian dynamics
   * (BD) simulations, named McSimulator and BdSimulator, provide more
   * specialized algorithms and data structures needed by these two
   * sampling methods.
   *
   * The Simulator class provides functions to compute and diagonalize
   * a projected chi matrix, functions to access components of several
   * types of fields in a basis of eigenvectors of the projected chi
   * matrix, and functions to compute and return contributions to the
   * field theoretic Hamiltonian and its components.
   *
   * The analyzeChi function constructs and diagonalizes the projected
   * chi matrix. This is a singular nMonomer x nMonomer matrix defined
   * by evaluating the orthogonal projection of the chi matrix into the
   * subspace of fluctuations that preserves total monomer concentration.
   * The eigenvalues and eigenvectors of this matrix are accessed via
   * the chiEvals and chiEvecs functions, respectively.
   *
   * The functions computeWc, computeCc and computeDc compute components
   * of various types of multi-component fields (i.e., fields that are
   * associated with a monomer type index) in a basis of eigenvectors of
   * the projected chi matrix. Names such as wc, cc and dc that end with
   * a suffix "c" refer to components of multi-component fields that are
   * defined using this eigenvector basis.
   *
   * <b> Usage </b>: A specialization of Rp::Simulator\<D, T\> serves 
   * as a base class for for each Simulation\<D\> class defined in 
   * namespaces Rpc and Rpg, for D=1, 2 and 3. In this usage, template 
   * parameter T is an instance of a template \<int D\> class Types that
   * is defined in each of these two program-level namespaces, and that 
   * defines a set of typename aliases for classes used in each such
   * namespace, for the relevant value of D.
   *
   * <b> Template parameters and typename aliases </b>:
   *
   *    D - integer dimensionality of space (D=1, 2, or 3)
   *    T - Types class, Rpc::Types<D> or Rpg::Types<D>)
   *
   * \see \ref rp_BdSimulator_page  (manual page)
   * \see \ref rp_McSimulator_page  (manual page)
   * \ingroup Rp_Fts_Simulator_Module
   */
   template <int D, class T>
   class Simulator : public ParamComposite
   {

   public:

      /// Container for a real-valued periodic field
      using RFieldT = typename T::RField;

      /// \name Construction, destruction and initialization
      ///@{

      /**
      * Constructor.
      *
      * \param system  parent System
      * \param simulator  enclosing instance of a subclass
      */
      Simulator(System<D,T>& system);

      /**
      * Destructor.
      */
      ~Simulator();

      // Prohibit copying and assignment.
      Simulator(Simulator<D,T> const &) = delete;
      Simulator<D,T>& operator = (Simulator<D,T> const &) = delete;

      /**
      * Read parameters for a simulation.
      *
      * The default implementation reads an optional random seed, an
      * optional Compressor block, an optional Perturbation, and an
      * optional Ramp, in that order. This default is intended to be
      * used only for unit testing.
      *
      * \param in input parameter stream
      */
      virtual void readParameters(std::istream &in);

      /**
      * Allocate required memory during initialization.
      *
      * This function must be called by the readParameters method of any
      * subclass. Declared as public to also allow use in unit tests.
      *
      * Preconditions: Values of nMonomer and the mesh dimensions must
      * be defined in Mixture and Domain members of the parent System on 
      * entry.
      */
      void allocate();

      ///@}
      /// \name Primary Actions: Simulation and Analysis
      ///@{

      /**
      * Perform a field theoretic Monte-Carlo simulation.
      *
      * Perform a field theoretic simulation of nSteps using the partial
      * saddle-point approximation.
      *
      * The default implemention is a do-nothing placeholder that throws
      * an error if called, and must be overridden by subclasses.
      *
      * \param nStep  number of simulation steps
      */
      virtual void simulate(int nStep);

      /**
      * Read and analyze a trajectory file.
      *
      * This function uses an instance of the TrajectoryReader class
      * specified by the "classname" argument to read a trajectory file
      * with the specified filename. The function opens the file,
      * performs the analysis, and closes the file before returning.
      *
      * The default implemention is a do-nothing placeholder that throws
      * an error if called, and must be overridden by subclasses.
      *
      * \param min  first frame number
      * \param max  last frame number
      * \param classname  name of TrajectoryReader class
      * \param filename  name of trajectory file
      */
      virtual void analyze(int min, int max,
                           std::string classname,
                           std::string filename);

      /**
      * Clear field eigen-components and Hamiltonian components.
      *
      * On return from this function, hasHamiltonian(), hasWc(), hasCc(),
      * and hasDc() will all return false.
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
      * Output MDE counter.
      *
      * Output the number of times the modified diffusion equation has
      * been solved.
      *
      * \param out  output stream
      */
      virtual void outputMdeCounter(std::ostream& out) const;

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
      /// \name Projected Chi Matrix
      ///@{

      /**
      * Perform eigenvalue analysis of projected chi matrix.
      *
      * Uses a chi matrix obtained from the Interaction member of
      * the parent System.
      */
      void analyzeChi();

      /**
      * Get an array of the eigenvalues of the projected chi matrix.
      *
      * The projected chi matrix is given by the matrix product P*chi*P,
      * where P is the symmetric projection matrix that projects onto
      * the subspace orthogonal to the vector e = (1,1,...,1). The
      * projected chi matrix is singular, and has a zero eigenvalue with
      * associated eigenvector e. By convention, this zero eigenvalue
      * and its eigenvector e are listed last, with index nMonomer - 1.
      */
      DArray<double> const & chiEvals() const;

      /**
      * Get a single eigenvalue of the projected chi matrix.
      *
      * \param a  index of eigenvalue (0, ... , nMonomer - 1)
      */
      double chiEval(int a) const;

      /**
      * Get the matrix of all eigenvectors of the projected chi matrix.
      *
      * This function returns the entire nMonomer x nMonomer matrix of the
      * eigenvectors of the projected chi matrix, in which each row is an
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
      DMatrix<double> const & chiEvecs() const;

      /**
      * Get one element of an eigenvector of the projected chi matrix.
      *
      * See documentation of chiEvecs(), which returns the entire matrix.
      *
      * \param a  eigenvector index (0, ..., nMonomer - 1)
      * \param i  monomoner type index (0, ..., nMonomer - 1)
      */
      double chiEvecs(int a, int i) const;

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
      * eigenvector a of the projected chi matrix, with the convention
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
      * approximation (PS-FTS). The value returned by this function 
      * is equal to the sum of values returned by idealHamiltonian(),
      * fieldHamiltonian() and perturbationHamiltonian() functions.
      */
      double hamiltonian() const;

      /**
      * Get an ideal contribution to the Hamiltonian.
      * 
      * The ideal Hamiltonian contribution returned by this function
      * is given by the quantity denoted by \f$ \tilde{H}_{\rm id} \f$, 
      * as defined \ref psfts_psa_pressure_sec "here", given by
      * \f[
      *    H_{\rm id}  = 
      *    \frac{V}{v} \sum_{a} \frac{\overline{\phi}_{a}}{N_{a}}
      *    \left [ 
      *     \ln \left ( \frac{\overline{\phi}_{a}}{Q_{a}} \right ) - 1
      *    \right ]
      *    - \frac{1}{v}\int W_{+}^{*}({\bf r})  \\
      * \f]
      * Here, \f$ a \f$ is an index for molecular species, while 
      * \f$ \overline{\phi}_{a} \f$, \f$ N_{a} \f$ and \f$ Q_{a} \f$ 
      * are the volume fraction, number of monomers per molecule (i.e.,
      * ratio of molecular to monomer volume), and the normalized 
      * single-molecule partition function for polymer or solvent 
      * species \f$ a \f$, respectively. 
      */
      double idealHamiltonian() const;

      /**
      * Get the quadratic field contribution to the Hamiltonian.
      *
      * The field Hamiltonian contribution returned by this function
      * is given by the quantity denoted by \f$ \tilde{H}_{\rm id} \f$, 
      * as defined \ref psfts_psa_pressure_sec "here", given by
      * \f[
      *    H_{\rm f}  =
      *    \frac{1}{v}
      *    \int \! d{\bf r} \; \left \{
      *    \sum_{\alpha=0}^{M-2}
      *    \frac{M (W_{\alpha} - S_{\alpha})^2 }{ 2 |\lambda_{\alpha}|}
      *    + \frac{S_{+}}{2} \right \} \quad,
      * \f]
      * where \f$ \alpha \f$ is an index for eigenvectors of 
      * the projected \f$ \chi \f$ matrix, 
      * \f$ S_{\alpha} = v_{\alpha i} \chi_{ij} e_{j}/M^{2} \f$,
      * and
      * \f$ S_{+} = S_{M-1} = e_{i} \chi_{ij} e_{j}/M^{2} \f$.
      */
      double fieldHamiltonian() const;

      /**
      * Get a perturbation to the standard Hamiltonian.
      *
      * A perturbation to the Hamiltonian, if any, is computed by an
      * associated Perturbation object. 
      *
      * When no perturbation exists, this function returns zero.
      */
      double perturbationHamiltonian() const;

      /**
      * Has the Hamiltonian been computed for the current w and c fields?
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
      * eigenvectors of the projected chi matrix. The component field
      * \f$ W_{a}({\bf r}) \f$ at grid point \f$ {\bf r} \f$ is given
      * using Einstein summation by
      * \f[
      *    W_{a}({\bf r}) =
      *    v_{ai} w_{i}({\bf r}) / M
      * \f]
      * where \f$ w_{i}({\bf r}) \f$ is the w-field associated with
      * monomer type \f$ i \f$, \f$ v_{ai} \f$ is eigenvector a of
      * the projected chi matrix, and M = nMonomer.
      */
      void computeWc();

      /**
      * Get all eigenvector components of the current w fields.
      *
      * This function returns a DArray of fields in which each field is
      * a chemical field component \f$ W_{a}({\bf r}) \f$ as defined in
      * the documentation of computeWc(), for a = 0, ..., nMonomer - 1.
      */
      DArray<RFieldT> const & wc() const;

      /**
      * Get one eigenvector component of the current w fields.
      *
      * See documentation of functions computeWc() and wc() for details.
      *
      * \param a eigenvector index in range 0 , ..., nMonomer -1
      */
      RFieldT const & wc(int a) const;

      /**
      * Are eigen-components of the current w fields valid ?
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
      * eigenvectors of the projected chi matrix.
      */
      void computeCc();

      /**
      * Get all eigenvector components of the current c fields.
      *
      * Each component \f$C_{a}({\bf r}) \f$ is a point-wise projection
      * of the monomer c fields onto a corresponding eigenvector of the
      * projected chi matrix. The resulting value \f$ C_{a}({\bf r}) \f$
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
      DArray<RFieldT> const & cc() const;

      /**
      * Get one eigenvector component of the current c fields.
      *
      * This returns a reference to a field \f$ C_{a}({\bf r}) \f$
      * as defined in the documentation of function cc().
      *
      * \param a eigenvector / eigenvalue index
      */
      RFieldT const & cc(int a) const;

      /**
      * Are eigen-components of the current c fields valid ?
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
      * by the member function wc(a).
      */
      DArray<RFieldT> const & dc() const;

      /**
      * Get one eigenvector component of the current d fields.
      *
      * \param i  eigenvector / eigenvalue index
      */
      RFieldT const & dc(int i) const;

      /**
      * Are the current d fields valid ?
      */
      bool hasDc() const;

      ///@}
      /// \name Save and Restore State
      ///@{

      /**
      * Save a copy of the current system state.
      *
      * This function and restoreState() are intended for use in the
      * implementation of field theoretic moves.  This function stores
      * the current w fields and the corresponding Hamiltonian value.
      * Current cc fields and dc fields are saved based on save policy.
      * This is normally the first step of a Monte-Carlo (MC) move, prior
      * to an attempted modification of the fields stored in the system
      * w field container.
      */
      void saveState();

      /**
      * Restore the system to the saved state.
      *
      * This function and saveState() are intended to be used together
      * in the implementation of FTS moves. If an attempted Monte-Carlo
      * move is rejected, or if the compressor fails to converge after
      * any attempted FTS move, restoreState() is called to restore the
      * fields and Hamiltonian value that were saved by a previous call
      * to the function saveState().
      */
      void restoreState();

      /**
      * Clear the saved copy of the system state.
      *
      * This function, restoreState(), and saveState() are intended to
      * be used together in the implementation of reversible FTS moves.
      * If an attempted move is accepted, clearState() is called to
      * clear the stored state and indicate acceptance.
      */
      void clearState();

      ///@}
      /// \name Miscellaneous (Accessors and Boolean Flags)
      ///@{

      /**
      * Get the parent system by reference.
      */
      System<D,T>& system();

      /**
      * Get the scalar random number generator by reference.
      */
      Random& random();

      /**
      * Get the vector random number generator by reference.
      */
      typename T::VecRandom& vecRandom();

      /**
      * Does this Simulator have a Compressor?
      */
      bool hasCompressor() const;

      /**
      * Get the Compressor by const reference.
      */
      typename T::Compressor const & compressor() const;

      /**
      * Get the Compressor by non-const reference.
      */
      typename T::Compressor& compressor();

      /**
      * Does this Simulator have a Perturbation?
      */
      bool hasPerturbation() const;

      /**
      * Get a Perturbation by const reference.
      */
      typename T::Perturbation const & perturbation() const;

      /**
      * Get a Perturbation by non-const reference.
      */
      typename T::Perturbation& perturbation();

      /**
      * Does this Simulator have a Ramp?
      */
      bool hasRamp() const;

      /**
      * Get a Ramp by const reference.
      */
      typename T::Ramp const & ramp() const;

      /**
      * Get a Ramp by non-const reference.
      */
      typename T::Ramp& ramp();

      ///@}

   protected:

      // Protected member functions

      /**
      * Optionally read a random seed and initialize RNGs.
      *
      * \param in  input parameter stream
      */
      void readRandomSeed(std::istream& in);

      /**
      * Get the Compressor factory by reference.
      */
      typename T::CompressorFactory& compressorFactory();

      /**
      * Optionally read a Compressor parameter file block.
      *
      * If isEnd is true on entry, this function returns without
      * attempting to read the Compressor block.
      *
      * \param in  input parameter stream
      * \param isEnd  Has the end bracket of the Simulator block been read?
      */
      void readCompressor(std::istream& in, bool& isEnd);

      /**
      * Get the Perturbation factory by reference.
      */
      typename T::PerturbationFactory& perturbationFactory();

      /**
      * Optionally read a Perturbation parameter file block.
      *
      * If isEnd is true on entry, this function returns without
      * attempting to read the Perturbation block.
      *
      * \param in  input parameter stream
      * \param isEnd  Has the end bracket of the Simulator block been read?
      */
      void readPerturbation(std::istream& in, bool& isEnd);

      /**
      * Set the associated Perturbation.
      *
      * \param ptr pointer to a new Perturbation
      */
      void setPerturbation(typename T::Perturbation* ptr);

      /**
      * Get the Ramp factory by reference.
      */
      typename T::RampFactory& rampFactory();

      /**
      * Optionally read a Ramp parameter file block.
      *
      * If isEnd is true on entry, this function returns without
      * attempting to read the Ramp block.
      *
      * \param in  input parameter stream
      * \param isEnd  Has the end bracket of the Simulator block been read?
      */
      void readRamp(std::istream& in, bool& isEnd);

      /**
      * Set the associated Ramp.
      *
      * \param ptr pointer to a new Ramp
      */
      void setRamp(typename T::Ramp* ptr);

      /**
      * Get the SimState stored internal state by reference.
      *
      * The T::SimState object is used to store the previous internal
      * state of the system. This allows restoration after a rejected 
      * MC move or failure of the Compressor to converge during either
      * a BD or MC move.
      */
      SimState<D,T>& state();

      // Protected data members

      /**
      * Eigenvector components of w fields on a real space grid.
      *
      * Each field component corresponds to a point-wise projection of
      * the monomer w fields onto an eigenvector of the projected chi
      * matrix. The number of components is equal to the number of
      * monomer types, nMonomer. The last component is a pressure-like
      * field.
      */
      DArray<RFieldT> wc_;

      /**
      * Eigenvector components of c fields on a real space grid.
      *
      * Each field component corresponds to a point-wise projection of
      * the monomer c fields onto an eigenvector of the projected chi
      * matrix. The number of components is equal to the number of
      * monomer types, nMonomer. The last component must satisfy an
      * incompressibility constraint.
      */
      DArray<RFieldT> cc_;

      /**
      * Components of d fields on a real space grid.
      *
      * Each field component is the functional derivative of H[W]
      * with respect to one eigenvector w-field component.
      */
      DArray<RFieldT> dc_;

      /**
      * Previous state saved at the beginning of a step.
      *
      * This data structure is used to restore a previous state if the
      * compressor fails to converge or if a MC move is rejected.
      */
      mutable SimState<D,T> state_;

      /**
      * Total field theoretic Hamiltonian H[W] (extensive value).
      */
      double hamiltonian_;

      /**
      * Ideal gas contribution (-lnQ) to Hamiltonian.
      */
      double idealHamiltonian_;

      /**
      * Quadratic field contribution to Hamiltonian.
      */
      double fieldHamiltonian_;

      /**
      * Perturbation contribution to the Hamiltonian.
      *
      * A perturbation Hamiltonian component, if any, is computed by an
      * associated Perturbation object and added to the ideal and field
      * components to obtain the total hamiltonian_ value.
      */
      double perturbationHamiltonian_;

      /**
      * Step counter - number of steps for which the compressor converged.
      *
      * Steps for which the compressor fails to converge are returned to
      * the previous state so that another random displacement can be
      * chosen. Attempted MC moves for which the compressor converges
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
      * Random number generator seed (input value).
      */
      long seed_;

      // Boolean status flags

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

   private:

      /**
      * Projected chi matrix
      *
      * Projected matrix chiP_ = P*chi*P, where P = I - e e^{T} / M is
      * a projection matrix that projects onto the subspace orthogonal
      * to the vector e = [1, ... , 1]^{T}, where M = nMonomer.
      */
      DMatrix<double> chiP_;

      /**
      * Eigenvectors of the projected chi matrix.
      *
      * Each row (identified by first index) is an eigenvector.
      * The last eigenvector, with index nMonomer - 1, is always the
      * vector e = [1, 1, ...., 1]. Distinct eigenvectors are orthogonal.
      * Eigenvectors are normalized such that the sum of the square of the
      * elements is equal to nMonomer.
      */
      DMatrix<double> chiEvecs_;

      /**
      * Eigenvalues of the projected chi matrix.
      *
      * The last eigenvalue, with index nMonomer - 1, is always zero.
      */
      DArray<double>  chiEvals_;

      /**
      * Components of vector s = chi*e in a basis of eigenvectors.
      *
      * Component sc_[a] is equal to v_{a}^{T} chi e / M^2, where
      * e = [1 1 ... 1]^{T}, v_{a}^{T} is a row vector representation
      * of eigenvector a of the projected chi matrix, given by element
      * a of chiEvecs_, and M = nMonomer.
      */
      DArray<double>  sc_;

      /**
      * Field used as temporary work space.
      */
      mutable RFieldT tmpField_;

      // Pointers to associated objects

      /**
      * Pointer to the parent system.
      */
      System<D,T>* systemPtr_;

      /**
      * Pointer to a scalar random number generator.
      */
      Random* randomPtr_;

      /**
      * Pointer to a vector random number generator.
      */
      typename T::VecRandom* vecRandomPtr_;

      /**
      * Pointer to a Compressor factory.
      */
      typename T::CompressorFactory* compressorFactoryPtr_;

      /**
      * Pointer to a compressor.
      */
      typename T::Compressor* compressorPtr_;

      /**
      * Pointer to a Perturbation factory.
      */
      typename T::PerturbationFactory* perturbationFactoryPtr_;

      /**
      * Pointer to a Perturbation.
      */
      typename T::Perturbation* perturbationPtr_;

      /**
      * Pointer to a Ramp factory.
      */
      typename T::RampFactory* rampFactoryPtr_;

      /**
      * Pointer to a Ramp.
      */
      typename T::Ramp* rampPtr_;

      /**
      * Has required memory been allocated?
      */
      bool isAllocated_;

   };

   // Inline functions

   // Access to associated objects via pointers

   // Get the parent System by reference.
   template <int D, class T> inline 
   System<D,T>& Simulator<D,T>::system()
   {
      UTIL_ASSERT(systemPtr_);
      return *systemPtr_;
   }

   // Get the scalar random number generator by reference.
   template <int D, class T> inline 
   Random& Simulator<D,T>::random()
   {
      UTIL_ASSERT(randomPtr_);
      return *randomPtr_;
   }

   // Get the vector random number generator by reference.
   template <int D, class T> inline 
   typename T::VecRandom& Simulator<D,T>::vecRandom()
   {
      UTIL_ASSERT(vecRandomPtr_);
      return *vecRandomPtr_;
   }

   // Does this Simulator have a Compressor?
   template <int D, class T> inline 
   bool Simulator<D,T>::hasCompressor() const
   {  return (bool)compressorPtr_; }

   // Get the Compressor by non-const reference.
   template <int D, class T> inline 
   typename T::Compressor& Simulator<D,T>::compressor()
   {
      UTIL_CHECK(compressorPtr_);
      return *compressorPtr_;
   }

   // Get the Compressor by const reference.
   template <int D, class T> inline 
   typename T::Compressor const & Simulator<D,T>::compressor() const
   {
      UTIL_CHECK(compressorPtr_);
      return *compressorPtr_;
   }

   // Does this Simulator have an associated Perturbation?
   template <int D, class T> inline 
   bool Simulator<D,T>::hasPerturbation() const
   {  return (bool)perturbationPtr_; }

   // Get a Perturbation by const reference.
   template <int D, class T> inline 
   typename T::Perturbation const & Simulator<D,T>::perturbation() const
   {
      UTIL_CHECK(perturbationPtr_);
      return *perturbationPtr_;
   }

   // Get a Perturbation by non-const reference.
   template <int D, class T> inline 
   typename T::Perturbation& Simulator<D,T>::perturbation()
   {
      UTIL_CHECK(perturbationPtr_);
      return *perturbationPtr_;
   }

   // Does this Simulator have an associated Ramp?
   template <int D, class T> inline
   bool Simulator<D,T>::hasRamp() const
   {  return (bool)rampPtr_; }

   // Get a Ramp by const reference.
   template <int D, class T> inline 
   typename T::Ramp const & Simulator<D,T>::ramp() const
   {
      UTIL_CHECK(rampPtr_);
      return *rampPtr_;
   }

   // Get a Ramp by non-const reference.
   template <int D, class T> inline 
   typename T::Ramp& Simulator<D,T>::ramp()
   {
      UTIL_CHECK(rampPtr_);
      return *rampPtr_;
   }

   // Get the stored internal state by reference.
   template <int D, class T> inline 
   SimState<D,T>& Simulator<D,T>::state()
   {  return state_; }

   // Projected Chi Matrix

   // Return an array of eigenvalues of the projected chi matrix.
   template <int D, class T> inline 
   DArray<double> const & Simulator<D,T>::chiEvals() const
   {  return chiEvals_; }

   // Return a single eigenvalue of the projected chi matrix.
   template <int D, class T> inline 
   double Simulator<D,T>::chiEval(int a) const
   {  return chiEvals_[a]; }

   // Return a matrix of eigenvectors of the projected chi matrix.
   template <int D, class T> inline 
   DMatrix<double> const & Simulator<D,T>::chiEvecs() const
   {  return chiEvecs_; }

   // Return an element of an eigenvector of the projected chi matrix.
   template <int D, class T> inline 
   double Simulator<D,T>::chiEvecs(int a, int i) const
   {  return chiEvecs_(a, i); }

   // Return an array of values of vector S.
   template <int D, class T> inline 
   DArray<double> const & Simulator<D,T>::sc() const
   {  return sc_; }

   // Return one component of vector S.
   template <int D, class T> inline 
   double Simulator<D,T>::sc(int a) const
   {  return sc_[a]; }

   // Hamiltonian and its components

   // Has the Hamiltonian been computed for the current w fields ?
   template <int D, class T> inline 
   bool Simulator<D,T>::hasHamiltonian() const
   {  return hasHamiltonian_; }

   // Get the precomputed total Hamiltonian.
   template <int D, class T> inline 
   double Simulator<D,T>::hamiltonian() const
   {
      UTIL_CHECK(hasHamiltonian_);
      return hamiltonian_;
   }

   // Get the ideal gas component of the precomputed Hamiltonian.
   template <int D, class T> inline 
   double Simulator<D,T>::idealHamiltonian() const
   {
      UTIL_CHECK(hasHamiltonian_);
      return idealHamiltonian_;
   }

   // Get the harmonic field component of the precomputed Hamiltonian.
   template <int D, class T> inline 
   double Simulator<D,T>::fieldHamiltonian() const
   {
      UTIL_CHECK(hasHamiltonian_);
      return fieldHamiltonian_;
   }

   // Get the perturbation component of the precomputed Hamiltonian.
   template <int D, class T> inline 
   double Simulator<D,T>::perturbationHamiltonian() const
   {
      UTIL_CHECK(hasHamiltonian_);
      return perturbationHamiltonian_;
   }

   // Fields

   // Have eigenvector components of the current w fields been computed?
   template <int D, class T> inline 
   bool Simulator<D,T>::hasWc() const
   {  return hasWc_; }

   // Return all eigencomponents of the w fields.
   template <int D, class T> inline 
   DArray<typename T::RField> const & Simulator<D,T>::wc() const
   {  return wc_; }

   // Return a single eigenvector component of the w fields.
   template <int D, class T> inline 
   typename T::RField const & Simulator<D,T>::wc(int a) const
   {  return wc_[a]; }

   // Have eigenvector components of the current c fields been computed?
   template <int D, class T> inline 
   bool Simulator<D,T>::hasCc() const
   {  return hasCc_; }

   // Return all eigenvector components of the current c fields.
   template <int D, class T> inline 
   DArray<typename T::RField> const & Simulator<D,T>::cc() const
   {  return cc_; }

   // Return a single eigenvector component of the current c fields.
   template <int D, class T> inline 
   typename T::RField const & Simulator<D,T>::cc(int a) const
   {  return cc_[a]; }

   // Have eigenvector components of the current d fields been computed?
   template <int D, class T>
   inline bool Simulator<D,T>::hasDc() const
   {  return hasDc_; }

   // Return all eigenvector components of the current d fields.
   template <int D, class T> inline 
   DArray<typename T::RField> const & Simulator<D,T>::dc() const
   {  return dc_; }

   // Return a single eigenvector component of the current d fields.
   template <int D, class T> inline 
   typename T::RField const & Simulator<D,T>::dc(int a) const
   {  return dc_[a]; }

   // Return the current converged simulation step index.
   template <int D, class T>
   inline long Simulator<D,T>::iStep()
   {  return iStep_; }

   // Return the current total simulation step index.
   template <int D, class T>
   inline long Simulator<D,T>::iTotalStep()
   {  return iTotalStep_; }

}
}
#endif

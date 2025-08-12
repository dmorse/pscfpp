#ifndef RPC_SIMULATOR_H
#define RPC_SIMULATOR_H

/*
* PSCF - Polymer Self-Consistent Field
*
* Copyright 2015 - 2025, The Regents of the University of Minnesota
* Distributed under the terms of the GNU General Public License.
*/

#include <util/param/ParamComposite.h>     // base class

#include <rpc/fts/simulator/SimState.h>    // member
#include <prdc/cpu/RField.h>               // member (template arg)
#include <util/random/Random.h>            // member
#include <util/containers/DArray.h>        // member (template)
#include <util/containers/DMatrix.h>       // member (template)

namespace Pscf {
namespace Rpc {

   // Forward declarations 
   template <int D> class System;
   template <int D> class Compressor;
   template <int D> class CompressorFactory;
   template <int D> class Perturbation;
   template <int D> class PerturbationFactory;
   template <int D> class Ramp;
   template <int D> class RampFactory;

   using namespace Util;
   using namespace Prdc;
   using namespace Prdc::Cpu;

   /**
   * Field theoretic simulator (base class).
   *
   * The Simulator base class provides tools needed in field-theoretic
   * simulations that are based on a partial saddle-point approximation.
   * Subclasses designed for field theoretic Monte Carlo (MC) and Brownian
   * dynamics (BD) simulations, named McSimulator and BdSimulator, provide
   * more specialized algorithms and data structures needed by these two
   * sampling methods. 
   *
   * The Simulator class provides functions to compute and diagonalze a 
   * projected chi matrix, functions to access components of several 
   * types of fields in a basis of eigenvectors of the projected chi 
   * matrix, and functions to compute and return contributions to the 
   * field theoretic Hamiltonian.
   *
   * The analyzeChi function constructs and diagonalizes the projected
   * chi matrix. This is a singular nMonomer x nMonomer matrix defined 
   * by evaluating the orthogonal projection of the chi matrix into the 
   * subspace of fluctuations that preserves total monomer concentration. 
   * The eigenvalues and eigenvectors of this matrix are accessed via 
   * the chiEvals and chiEvecs functions, respectively.
   *
   * The computeWc, computeCc and computeDc functions compute components
   * of various types of multi-component fields (i.e., fields that are
   * associated with a monomer type index) in a basis of eigenvectors of 
   * the projected chi matrix. Names such as wc, cc and dc that end with 
   * a suffix "c" refer to components of multi-component fields that are 
   * defined using this eigenvector basis. 
   *
   * \ingroup Rpc_Fts_Module
   */
   template <int D>
   class Simulator : public ParamComposite
   {

   public:

      Simulator() = delete;

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
      * The default implementation reads a Compressor block, an optional 
      * random seed, and optional Perturbation, and an optional Ramp, in
      * that order. This is intended to be used by subclasses designed for 
      * and BD simulations to read the shared initial parts of the 
      * parameter block format, to which blocks for MC moves or a BD step 
      * and analyzers can be added by subclasses.
      *
      * \param in input parameter stream
      */
      virtual void readParameters(std::istream &in);

      /// \name Primary Actions: Simulation and Analysis
      ///@{

      /**
      * Perform a field theoretic Monte-Carlo simulation.
      *
      * Perform a field theoretic simulation of nSteps using the
      * partial saddle-point approximation.
      *
      * The default implemention is a do-nothing placeholder that throws
      * an error if called, and must be re-implemented by subclasses.
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
      * Output MDE counter.
      *
      * Output the number of times the modified diffusion equation has been 
      * solved.
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
      * \param a index of eigenvalue (0, ... , nMonomer - 1)
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
      * approximation (PS-FTS).
      */
      double hamiltonian() const;

      /**
      * Get ideal gas contribution to the Hamiltonian.
      */
      double idealHamiltonian() const;

      /**
      * Get the quadratic field contribution to the Hamiltonian.
      */
      double fieldHamiltonian() const;

      /**
      * Get the perturbation to the standard Hamiltonian (if any).
      *
      * A perturbation to the Hamiltonian, if any, is computed by an 
      * associated Perturbation object. When a perturbation exists, 
      * as indicated by hasPerturbation(), the perturbationHamiltonian
      * component is added to the idealHamiltonian and fieldHamiltonian
      * components to obtain the total value that is returned by the
      * hamiltonian() member function.
      */
      double perturbationHamiltonian() const;

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
      DArray< RField<D> > const & wc() const;

      /**
      * Get one eigenvector component of the current w fields.
      *
      * See documentation of functions computeWc() and wc() for details.
      *
      * \param a eigenvector index in range 0 , ..., nMonomer -1
      */
      RField<D> const & wc(int a) const;

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
      DArray< RField<D> > const & cc() const;

      /**
      * Get one eigenvector component of the current c fields.
      *
      * This returns a reference to a field \f$ C_{a}({\bf r}) \f$
      * as defined in the documentation of function cc().
      *
      * \param a eigenvector / eigenvalue index
      */
      RField<D> const & cc(int a) const;

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
      DArray< RField<D> > const & dc() const;

      /**
      * Get one eigenvector component of the current d fields.
      *
      * \param i eigenvector / eigenvalue index
      */
      RField<D> const & dc(int i) const;

      /**
      * Are the current d fields valid ?
      */
      bool hasDc() const;
      
      ///@}
      /// \name Utilities for moves
      ///@{
      
      /**
      * Save a copy of the fts move state.
      *
      * This function and restoreState() are intended for use 
      * in the implementation of field theoretic moves. 
      * This function stores the current w fields and the corresponding
      * Hamiltonian value. Current cc fields and dc fields are saved 
      * based on save policy. This is normally the first step of a fts
      * move, prior to an attempted modification of the fields stored
      * in the system w field container.
      */
      void saveState();
      
      /**
      * Restore the saved copy of the fts move state.
      *
      * This function and saveState() are intended to be used
      * together in the implementation of fts moves. If an
      * attempted Monte-Carle move is rejected or an fts move 
      * fails to converge restoreState() is called to restore 
      * the fields and Hamiltonian value that were saved
      * by a previous call to the function saveState().
      */
      void restoreState();
      
      /**
      * Clear the saved copy of the fts state.
      *
      * This function, restoreState(), and saveState() are intended
      * to be used together in the implementation of reversible fts moves. 
      * If an attempted move is accepted, clearState() is called to 
      * indicate the need to recompute some quantities.
      */
      void clearState();
      
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
      * Does this Simulator have a Compressor?
      */
      bool hasCompressor() const;

      /**
      * Get the compressor by non-const reference.
      */
      Compressor<D>& compressor();

      /**
      * Get the compressor by const reference.
      */
      Compressor<D> const & compressor() const;

      /**
      * Does this Simulator have a Perturbation?
      */
      bool hasPerturbation() const;

      /**
      * Get the perturbation factory by non-const reference.
      */
      Perturbation<D>& perturbation();

      /**
      * Get the associated Perturbation by const reference.
      */
      Perturbation<D> const & perturbation() const;

      /**
      * Does this Simulator have a Ramp?
      */
      bool hasRamp() const;

      /**
      * Get the associated Ramp by const reference.
      */
      Ramp<D> const & ramp() const;
      
      /**
      * Get the ramp by non-const reference.
      */
      Ramp<D>& ramp();

      ///@}

   protected:

      // Protected member functions

      using Util::ParamComposite::setClassName;

      /**
      * Optionally read a random number generator seed.
      *
      * \param in input parameter stream
      */
      void readRandomSeed(std::istream& in);

      /**
      * Get the compressor factory by reference.
      */
      CompressorFactory<D>& compressorFactory();

      /**
      * Read the compressor block of the parameter file. 
      *
      * If isEnd it true on entry, the function returns immediately. 
      *
      * \param in  input parameter stream
      * \param isEnd  Was the end bracket of the Simulator block read?
      */
      void readCompressor(std::istream& in, bool& isEnd);

      /**
      * Get the perturbation factory by reference.
      */
      PerturbationFactory<D>& perturbationFactory();

      /**
      * Optionally read an associated perturbation.
      *
      * If isEnd it true on entry, this function returns immediately. 
      *
      * \param in  input parameter stream
      * \param isEnd  Was the end bracket of the Simulator block read?
      */
      void readPerturbation(std::istream& in, bool& isEnd);

      /**
      * Set the associated perturbation.
      *
      * \param ptr pointer to a new Perturbation<D> object.
      */
      void setPerturbation(Perturbation<D>* ptr);

      /**
      * Get the ramp factory by reference.
      */
      RampFactory<D>& rampFactory();

      /**
      * Optionally read an associated ramp.
      *
      * If isEnd it true on entry, this function returns immediately. 
      *
      * \param in  input parameter stream
      * \param isEnd  Was the end bracket of the Simulator block read?
      */
      void readRamp(std::istream& in, bool& isEnd);

      /**
      * Set the associated ramp.
      *
      * \param ptr pointer to a new Ramp<D> object.
      */
      void setRamp(Ramp<D>* ptr);

      // Protected data members

      /**
      * Random number generator
      */
      Random random_;

      /**
      * Eigenvector components of w fields on a real space grid.
      *
      * Each field component corresponds to a point-wise projection of 
      * the monomer w fields onto an eigenvector of the projected chi 
      * matrix. The number of components is equal to the number of
      * monomer types, nMonomer. The last component is a pressure-like
      * field.
      */
      DArray< RField<D> > wc_;

      /**
      * Eigenvector components of c fields on a real space grid.
      *
      * Each field component corresponds to a point-wise projection of 
      * the monomer c fields onto an eigenvector of the projected chi 
      * matrix. The number of components is equal to the number of 
      * monomer types, nMonomer. The last component must satisfy an
      * incompressibility constraint.
      */
      DArray< RField<D> > cc_;

      /**
      * Components of d fields on a real space grid.
      *
      * Each field component is the functional derivative of H[W]
      * with respect to one eigenvector w-field component.
      */
      DArray< RField<D> > dc_;
      
      /**
      * Previous state saved during at the beginning of a step.
      * 
      * This data structure is used to restore a previous state if
      * the compressor fails to converge or if a MC move is rejected.
      */
      mutable SimState<D> state_;

      /**
      * Total field theoretic Hamiltonian H[W] (extensive value).
      */
      double hamiltonian_;

      /**
      * Ideal gas contribution (-lnQ) to Hamiltonian H[W]
      */
      double idealHamiltonian_;

      /**
      * Field contribution (H_W) to Hamiltonian
      */
      double fieldHamiltonian_;

      /**
      * Perturbation to the standard Hamiltonian (if any).
      *
      * A perturbation to the Hamiltonian, if any, is computed by an 
      * associated Perturbation object and added to the ideal and field
      * components to obtain the total hamiltonian_ value. 
      */
      double perturbationHamiltonian_;

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
      * Eigenvectors of the projected chi matrix chiP_.
      *
      * Each row (identified by first index) is an eigenvector. 
      * The last eigenvector, with index nMonomer - 1, is always the
      * vector e = [1, 1, ...., 1]. Distinct eigenvectors are orthogonal.
      * Eigenvectors normalized such that the sum of the square of the 
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
      * Component sc_[a] is equal to v_{a} chi e / M^2, where
      * e = [1 1 ... 1]^{T}, v_{a} is a row vector representation of 
      * eigenvector a of the projected chi matrix, given by row a of 
      * chiEvecs_, and M = nMonomer.
      */
      DArray<double>  sc_;
      
      /**
      * Pointer to the parent system.
      */
      System<D>* systemPtr_;

      /**
      * Pointer to the compressor factory object.
      */
      CompressorFactory<D>* compressorFactoryPtr_;

      /**
      * Pointer to a compressor.
      */
      Compressor<D>* compressorPtr_;

      /**
      * Pointer to the perturbation Factory.
      */
      PerturbationFactory<D>* perturbationFactoryPtr_;

      /**
      * Pointer to a perturbation (if any)
      */
      Perturbation<D>* perturbationPtr_;

      /**
      * Pointer to the Ramp Factory.
      */
      RampFactory<D>* rampFactoryPtr_;

      /**
      * Pointer to a Ramp (if any)
      */
      Ramp<D>* rampPtr_;

      /**
      * Has required memory been allocated?
      */
      bool isAllocated_;

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

   // Get the compressor factory by reference.
   template <int D>
   inline CompressorFactory<D>& Simulator<D>::compressorFactory()
   {
      UTIL_CHECK(compressorFactoryPtr_);  
      return *compressorFactoryPtr_; 
   }

   // Does this Simulator have an associated Compressor?
   template <int D>
   inline bool Simulator<D>::hasCompressor() const
   {  return (bool)compressorPtr_; }
   
   // Get the Compressor by non-const reference.
   template <int D>
   inline Compressor<D>& Simulator<D>::compressor()
   {
      UTIL_CHECK(hasCompressor());
      return *compressorPtr_;
   }

   // Get the Compressor by const reference.
   template <int D>
   inline Compressor<D> const& Simulator<D>::compressor() const
   {
      UTIL_CHECK(hasCompressor());
      return *compressorPtr_;
   }

   // Does this Simulator have an associated Perturbation?
   template <int D>
   inline bool Simulator<D>::hasPerturbation() const
   {  return (bool)perturbationPtr_; }
   
   // Get the perturbation (if any) by const reference.
   template <int D>
   inline Perturbation<D> const & Simulator<D>::perturbation() const
   {
      UTIL_CHECK(perturbationPtr_);  
      return *perturbationPtr_; 
   }

   // Get the perturbation (if any) by non-const reference.
   template <int D>
   inline Perturbation<D>& Simulator<D>::perturbation()
   {
      UTIL_CHECK(perturbationPtr_);  
      return *perturbationPtr_; 
   }

   // Get the perturbation factory by reference.
   template <int D>
   inline PerturbationFactory<D>& Simulator<D>::perturbationFactory()
   {
      UTIL_CHECK(perturbationFactoryPtr_);  
      return *perturbationFactoryPtr_; 
   }

   // Does this Simulator have an associated Ramp?
   template <int D>
   inline bool Simulator<D>::hasRamp() const
   {  return (bool)rampPtr_; }
   
   // Get the ramp (if any) by const reference.
   template <int D>
   inline Ramp<D> const & Simulator<D>::ramp() const
   {
      UTIL_CHECK(rampPtr_);  
      return *rampPtr_; 
   }

   // Get the ramp (if any) by non-const reference.
   template <int D>
   inline Ramp<D>& Simulator<D>::ramp()
   {
      UTIL_CHECK(rampPtr_);  
      return *rampPtr_; 
   }

   // Get the ramp factory.
   template <int D>
   inline RampFactory<D>& Simulator<D>::rampFactory()
   {
      UTIL_CHECK(rampFactoryPtr_);  
      return *rampFactoryPtr_; 
   }

   // Projected Chi Matrix

   // Return an array of eigenvalues of the projected chi matrix.
   template <int D>
   inline DArray<double> const & Simulator<D>::chiEvals() const
   {  return chiEvals_; }

   // Return a single eigenvalue of the projected chi matrix.
   template <int D>
   inline double Simulator<D>::chiEval(int a) const
   {  return chiEvals_[a]; }

   // Return a matrix of eigenvectors of the projected chi matrix.
   template <int D>
   inline DMatrix<double> const & Simulator<D>::chiEvecs() const
   {  return chiEvecs_; }

   // Return an element of an eigenvector of the projected chi matrix.
   template <int D>
   inline double Simulator<D>::chiEvecs(int a, int i) const
   {  return chiEvecs_(a, i); }

   // Return the vector S in eigenvector basis.
   template <int D>
   inline DArray<double> const & Simulator<D>::sc() const
   {  return sc_; }

   // Return one component of vector S in an eigenvector basis.
   template <int D>
   inline double Simulator<D>::sc(int a) const
   {  return sc_[a]; }

   // Hamiltonian and its derivatives

   // Get the precomputed Hamiltonian.
   template <int D>
   inline double Simulator<D>::hamiltonian() const
   {
      UTIL_CHECK(hasHamiltonian_);
      return hamiltonian_;
   }

   // Get the ideal gas component of the precomputed Hamiltonian.
   template <int D>
   inline double Simulator<D>::idealHamiltonian() const
   {
      UTIL_CHECK(hasHamiltonian_);
      return idealHamiltonian_;
   }

   // Get the W field component of the precomputed Hamiltonian.
   template <int D>
   inline double Simulator<D>::fieldHamiltonian() const
   {
      UTIL_CHECK(hasHamiltonian_);
      return fieldHamiltonian_;
   }

   // Get the perturbation component of the precomputed Hamiltonian.
   template <int D>
   inline double Simulator<D>::perturbationHamiltonian() const
   {
      UTIL_CHECK(hasHamiltonian_);
      return perturbationHamiltonian_;
   }

   // Has the Hamiltonian been computed for the current w fields ?
   template <int D>
   inline bool Simulator<D>::hasHamiltonian() const
   {  return hasHamiltonian_; }

   // Fields

   // Return all eigencomponents of the w fields.
   template <int D>
   inline DArray< RField<D> > const & Simulator<D>::wc() const
   {  return wc_; }

   // Return a single eigenvector component of the w fields.
   template <int D>
   inline RField<D> const & Simulator<D>::wc(int a) const
   {  return wc_[a]; }

   // Have eigenvector components of current w fields been computed?
   template <int D>
   inline bool Simulator<D>::hasWc() const
   {  return hasWc_; }

   // Return all eigenvector components of the current c fields.
   template <int D>
   inline DArray< RField<D> > const & Simulator<D>::cc() const
   {  return cc_; }

   // Return a single eigenvector component of the current c fields.
   template <int D>
   inline RField<D> const & Simulator<D>::cc(int a) const
   {  return cc_[a]; }

   // Have eigenvector components of current c fields been computed?
   template <int D>
   inline bool Simulator<D>::hasCc() const
   {  return hasCc_; }

   // Return all eigenvector components of the current d fields.
   template <int D>
   inline DArray< RField<D> > const & Simulator<D>::dc() const
   {  return dc_; }

   // Return a single eigenvector component of the current d fields.
   template <int D>
   inline RField<D> const & Simulator<D>::dc(int a) const
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

   #ifndef RPC_SIMULATOR_TPP
   // Suppress implicit instantiation
   extern template class Simulator<1>;
   extern template class Simulator<2>;
   extern template class Simulator<3>;
   #endif

}
}
#endif

#ifndef RP_SYSTEM_H
#define RP_SYSTEM_H

/*
* PSCF - Polymer Self-Consistent Field
*
* Copyright 2015 - 2025, The Regents of the University of Minnesota
* Distributed under the terms of the GNU General Public License.
*/

#include <util/param/ParamComposite.h>   // base class
#include <pscf/chem/PolymerModel.h>      // member

// Forward declarations
namespace Util {
   class FileMaster;
   template <typename E, int N> class FSArray;
}
namespace Pscf {
   class Interaction;
   namespace Prdc {
      template <int D> class UnitCell;
      class Environment;
   }
   namespace Rp {
      template <int D, class T> class WFields;
      template <int D, class T> class CFields;
      template <int D, class T> class Mask;
      template <int D, class T> class Domain;
   }
}

namespace Pscf {
namespace Rp {

   // Namespaces that may be used implicitly
   using namespace Util;
   using namespace Pscf::Prdc;

   /**
   * Base class template for classes that represent a complete system.
   *
   * <b> Template parameters and typename aliases</b>:
   *
   *    D - integer dimensionality of space (D=1, 2, or 3)
   *    T - "Types" class collection of aliases for other classes
   *
   * <b> Usage </b>: A specialization of Rp::System\<D, T\> is used as
   * a base class for each System\<D\> class defined in the Rpc and Rpg
   * program-level namespaces, for D=1, 2, or 3.  In this usage, template
   * parameter T is an instance of a class template named Types that
   * is defined in each of these two namespaces.  For example, in the
   * Pscf::Rpc namespace, for each value of D, class Rpc::System\<D\>
   * is derived from the class Prdc::System\< D, Rpc::Types\<D\> >.
   *
   * For each such specialization, class T = Types\<D\> defines a set
   * of typename aliases for classes used in the relevant program-level
   * namespace, for the specified value of D.  For example, for each value
   * of D, the typename Rpc::Types\<D\>::Mixture is an alias for the type
   * Rpc::Mixture<D> that used to represent a mixture in the Rpc namespace
   * for systems of spatial dimension D. See the definitions of Rpc::Types
   * and Rpg::Types (src/rpc/system/Types.h and src/rpg/system/Types.h)
   * for lists of all of the typenames defined in these class templates.
   *
   * In the remainder of this documentation for the Rp::System template,
   * unqualified names such as "Mixture", "Iterator", etc. are used as
   * shorthand for typename aliases such as T::Mixture or T::Iterator
   * that are defined in the types class T (i.e., in Rpc::Types\<D\> or
   * Rpg::Types\<D\>), which are aliases for class names such as
   * Rpc::Mixture<D> or Rpg::Iterator<D>.
   *
   * <b> Class Components </b>:
   * A System object has (among other components):
   *
   *    - a Mixture (container for polymer and solvent solvers)
   *    - an %Interaction (list of binary interaction parameters)
   *    - a Domain (description of unit cell and discretization)
   *    - a CFields container of monomer concentration (c) fields
   *    - a WFields container of monomer chemical potential (w) fields
   *    - a WFields container of external (h) fields
   *    - a Mask field that defines an inhomogeneous density constraint
   *
   * The container of external fields and the Mask are only needed
   * to describe systems that are subjected to inhomogeneous imposed
   * environments (such as in thin films) and are otherwise left empty
   * and unused. The h().hasData() and mask().hasData() bool functions
   * may be queried to determine if these components are actually in
   * use.
   *
   * A System may also optionally have:
   *
   *    - an %Environment
   *    - an Iterator
   *    - a Sweep,
   *    - a Simulator
   *
   * Each optional component is constructed if and only if the parameter
   * file that is used to initialize the System contains a corresponding
   * optional parameter file block. An %Environment is used to generate
   * external and mask fields to describe inhomogeneous environments, and
   * is omitted in standard calculations for structures formed in a
   * homogeneous environment.  An Iterator is used for SCFT calculations.
   * A Sweep is used for "sweep" calculations that solve a sequence of
   * SCFT problems along a path through parameter space.  A Simulator is
   * only used for PS-FTS calculations,  i.e., field theoretic simulations
   * based on a partial saddle-point approximation. The Iterator and Sweep
   * objects is thus normally omitted for PS-FTS calculations, while the
   * Simulator object is normally omitted for SCFT calculations. The
   * hasEnvironment(), hasIterator(), hasSweep(), and hasSimulator() bool
   * member functions may be queried after processing of the parameter
   * file to determine which of these optional components exist.
   *
   * See also:
   * <ul>
   *  <li> \ref scft_param_pc_page   "Parameter File: SCFT" </li>
   *  <li> \ref psfts_param_page     "Parameter File: PS-FTS" </li>
   *  <li> \ref rp_System_page       "Parameter File: Full Format" </li>
   *  <li> \ref scft_command_pc_page "Command File Format" </li>
   * </ul>
   *
   * \ingroup Rp_System_Module
   */
   template <int D, class T>
   class System : public ParamComposite
   {

   public:

      // Suppress several compiler-generated member functions
      System() = delete;
      System(System<D,T> const &) = delete;
      System<D,T>& operator = (System<D,T> const & ) = delete;

      /// \name Lifetime Actions
      ///@{

      /**
      * Process command line options.
      *
      * This function takes the same arguments as any C/C++ main program
      * function. The arguments of the "main" function should be passed
      * to this function unaltered, to allow this function to process the
      * command line options passed to the program executable.
      *
      * \param argc number of command line arguments
      * \param argv array of pointers to command line arguments
      */
      void setOptions(int argc, char **argv);

      /**
      * Read input parameters (with opening and closing lines).
      *
      * \param in input parameter stream
      */
      virtual void readParam(std::istream& in);

      /**
      * Read input parameters from default param file.
      *
      * This function reads the parameter file set by the -p command
      * line option.
      */
      void readParam();

      /**
      * Read body of parameter block (without opening and closing lines).
      *
      * \param in input parameter stream
      */
      virtual void readParameters(std::istream& in);

      /**
      * Read and process commands from an input stream.
      *
      * \param in command script input stream
      */
      void readCommands(std::istream& in);

      /**
      * Read and process commands from the default command file.
      *
      * This function reads the parameter file set by the -c command
      * line option.
      */
      void readCommands();

      ///@}
      /// \name Field Theory Computations
      ///@{

      /**
      * Solve the modified diffusion equation once, without iteration.
      *
      * This function calls the Mixture::compute() function to solve
      * the statistical mechanics problem for a non-interacting system
      * subjected to the currrent system chemical potential fields. This
      * requires solution of the modified diffusion equation for all
      * polymers, computation of Boltzmann weights for all solvents,
      * computation of molecular partition functions for all species,
      * computation of concentration fields for all blocks and solvents,
      * and computation of overall concentrations for all monomer types.
      * This function does not compute the canonical (Helmholtz) free
      * energy or grand-canonical free energy (i.e., pressure).
      *
      * This function also computes the stress, by calling computeStress(),
      * if and only if the argument needStress is true.
      *
      * \pre  w().hasData() == true
      * \post c().hasData() == true
      * \post hasStress() == true iff needStress == true
      *
      * \param needStress  true if stress is needed, false otherwise
      */
      void compute(bool needStress = false);

      /**
      * Compute SCFT stress.
      *
      * This function computes the standard definition of stress maintained
      * by the Mixture class. If an %Environment exists, it also allows the
      * %Environment to compute a modified definition of the stress.
      *
      * \pre w().hasData() == true
      * \pre c().hasData() == true
      * \post hasStress() == true
      */
      void computeStress();

      /**
      * Iteratively solve a SCFT problem.
      *
      * This function calls the iterator to solve the SCFT problem for
      * the current system parameters, using the current chemical
      * potential fields and unit cell parameters as initial guesses.
      * Upon exit, c().hasData() == true whether or not convergence is
      * obtained to within the desired tolerance, but the SCFT Helmholtz
      * free energy and pressure are computed only if convergence is
      * successful.
      *
      * \pre hasIterator() == true
      * \pre w().hasData() == true
      * \pre w().isSymmetric() == true if the iterator is symmetric
      * \post c().hasData() == true
      * \post scft().hasData() == true upon successful convergence
      *
      * \return returns 0 for successful convergence, 1 for failure
      *
      * \param isContinuation  true if a continuation within a sweep
      */
      int iterate(bool isContinuation = false);

      /**
      * Sweep in parameter space, solving an SCF problem at each point.
      *
      * This function uses a Sweep object that was initialized in the
      * parameter file to solve the SCFT problem at a sequence of points
      * along a path in parameter space. The nature of this sequence
      * is determined by implementation of a subclass of Sweep and the
      * parameters passed to the sweep object in the parameter file.
      *
      * \pre hasSweep() == true
      * \pre All preconditions of the iterate() function must be satisfied
      */
      void sweep();

      /**
      * Perform a field theoretic simulation (PS-FTS).
      *
      * Perform a field theoretic simulation using the partial saddle-point
      * approximation (PS-FTS). The type of simulation (BD or MC) is
      * determined by the type of Simulator (BdSimulator or McSimulator)
      * that is created in the parameter file. The number of BD steps or
      * attempted MC moves to be performed is given by the parameter
      * "nStep".
      *
      * \pre Function hasSimulator() == true
      * \pre Function w().hasData() == true
      *
      * \param nStep  number of simulation (BD or MC) steps
      */
      void simulate(int nStep);

      /**
      * Mark c-fields and free energy as outdated or invalid.
      *
      * This function should be called whenever any of the inputs to the
      * solution of the modified diffusion equation are modified, including
      * the w fields, unit cell parameters, external fields, or mask. Upon
      * return, c().hasData(), scft().hasData(), and mixture().hasStress()
      * will all return false. If the system has an %Environment, then
      * environment().needsUpdate() will also return true.
      */
      void clearCFields();

      ///@}
      /// \name Unit Cell Modifiers
      ///@{

      /**
      * Set parameters of the associated unit cell.
      *
      * The lattice (i.e., lattice system type) set in the UnitCell<D>
      * input parameter must agree with a lattice enum value that was
      * set previously in the parameter file.
      *
      * If a space group has been declared but a basis has not yet been
      * initialized, then a symmetry-adapted basis will be constructed.
      *
      * \param unitCell  new UnitCell<D> (with new parameters)
      */
      void setUnitCell(UnitCell<D> const & unitCell);

      /**
      * Set parameters of the associated unit cell.
      *
      * The lattice type must have been set before this function is
      * called. The logical size of the FSArray<double, 6> "parameters"
      * array must match the expected number of parameters for the
      * current lattice type.
      *
      * If a space group has been declared but a basis has not yet been
      * initialized, then a symmetry-adapted basis will be constructed.
      *
      * \param parameters  array of new unit cell parameters
      */
      void setUnitCell(FSArray<double, 6> const & parameters);

      /**
      * Notify System members that unit cell parameters have been modified.
      *
      * This function should be called whenever the unit cell parameters
      * are modified. It calls functions mixture().clearUnitCellData(),
      * wavelist().clearUnitCellData(), clearCFields(), and, if an
      * %Environment exists, environment().reset().
      */
      void clearUnitCellData();

      ///@}
      /// \name Field Containers
      ///@{

      /**
      * Get the monomer concentration (c) fields (const).
      */
      CFields<D,T> const & c() const;

      /**
      * Get the chemical potential (w) fields (const).
      */
      WFields<D,T> const & w() const;

      /**
      * Get the chemical potential (w) fields (non-const).
      */
      WFields<D,T>& w();

      /**
      * Get the external potential (h) fields (const).
      */
      WFields<D,T> const & h() const;

      /**
      * Get the external potential (h) fields (non-const).
      */
      WFields<D,T>& h();

      /**
      * Get the mask (const).
      */
      Mask<D,T> const & mask() const;

      /**
      * Get the mask (non-const).
      */
      Mask<D,T>& mask();

      ///@}
      /// \name Component Object Accessors
      ///@{

      /**
      * Get the Mixture (const).
      */
      typename T::Mixture const & mixture() const;

      /**
      * Get the MixtureModifier (non-const).
      */
      typename T::MixtureModifier& mixtureModifier();

      /**
      * Get the %Interaction (const).
      */
      Interaction const & interaction() const;

      /**
      * Get the %Interaction (non-const).
      */
      Interaction& interaction();

      /**
      * Get the Domain (const).
      */
      Domain<D,T> const & domain() const;

      /**
      * Get the WaveList (non-const).
      *
      * This function provides direct non-const access to the WaveList.
      * The WaveList is owned by the Domain, and const (i.e., read-only)
      * access is also provided via the domain() member function. This
      * function exists to allow access to non-const member functions
      * that update data structures that are maintained by the WaveList
      * via a reference to the parent system.
      */
      typename T::WaveList & waveList();

      /**
      * Does this system have an %Environment?
      */
      bool hasEnvironment() const;

      /**
      * Get the %Environment (const).
      */
      Environment const & environment() const;

      /**
      * Get the %Environment (non-const).
      */
      Environment& environment();

      /**
      * Get the ScftThermo object (const).
      */
      typename T::ScftThermo const & scft() const;

      /**
      * Get the ScftThermo object (non-const).
      */
      typename T::ScftThermo& scft();

      /**
      * Does this system have an Iterator?
      */
      bool hasIterator() const;

      /**
      * Get the Iterator (const).
      */
      typename T::Iterator const & iterator() const;

      /**
      * Get the Iterator (non-const).
      */
      typename T::Iterator& iterator();

      /**
      * Does this system have a Sweep?
      */
      bool hasSweep() const;

      /**
      * Does this system have a Simulator?
      */
      bool hasSimulator() const;

      /**
      * Get the Simulator (const).
      */
      typename T::Simulator const & simulator() const;

      /**
      * Get the Simulator (non-const).
      */
      typename T::Simulator& simulator();

      /**
      * Get the FileMaster (const).
      */
      FileMaster const & fileMaster() const;

      /**
      * Get the FileMaster (non-const).
      *
      * Access (non-const) is used in some unit tests.
      */
      FileMaster& fileMaster();

      ///@}
      /// \name Property Output
      ///@{

      /**
      * Write partial parameter file to an ostream.
      *
      * This function writes the Mixture, %Interaction, and Domain blocks
      * of a parameter file, as well as any %Environment and Iterator
      * blocks, but omits any Sweep or Simulator blocks. The intent is
      * to produce an output during an SCFT sweep that only refers to
      * parameters relevant to a single state point, in a form that could
      * be used as a parameter file for a single SCFT calculation.
      *
      * \param out  output stream
      */
      void writeParamNoSweep(std::ostream& out) const;

      ///@}
      /// \name Timers
      ///@{

      /**
      * Write timer information to an output stream.
      *
      * \param out  output stream
      */
      void writeTimers(std::ostream& out) const;

      /**
      * Clear all timers.
      */
      void clearTimers();

      ///@}

   protected:

      /// \name Construction and Destruction
      ///@{

      /**
      * Constructor.
      *
      * When a specialization of System\<D,T\> is used as a base class
      * for a subclass defined in the Rpc or Rpg program-level namespace,
      * such as Rpc::System\<D\>, typename T::System is an alias for the
      * System subclass defined in Rpc or Rpg. In the constructor such a
      * derived class, the relevant instance of the derived class must
      * be passed to the Rp::System<D,T> base class constructor via the
      * standard "this" pointer. The address of this T::System subclass
      * instance is retained in the Rp::System base class instance by a
      * private member variable named systemPtr_ of type T::System*.
      * See definitions of the constructors for the Rpc::System and
      * Rpc::System class templates for examples of this usage.
      *
      * \param system  instance of System subclass
      */
      System(typename T::System& system);

      /**
      * Destructor.
      */
      ~System();

      ///@}

   private:

      /**
      * Pointer to enclosing instance of System subclass.
      */
      typename T::System* systemPtr_;

      // Pointers to associated sub-objects (owned by System)

      /**
      * Monomer concentration / volume fraction fields.
      */
      CFields<D,T>* cPtr_;

      /**
      * Chemical potential fields.
      */
      WFields<D,T>* wPtr_;

      /**
      * External potential fields.
      */
      WFields<D,T>* hPtr_;

      /**
      * Field to which the total density is constrained.
      */
      Mask<D,T>* maskPtr_;

      /**
      * Pointer to Mixture object (solves MDE for all species).
      */
      typename T::Mixture* mixturePtr_;

      /**
      * Pointer to MixtureModifier (public non-const interface for Mixture).
      */
      typename T::MixtureModifier* mixtureModifierPtr_;

      /**
      * Pointer to Domain object (unit cell, mesh, fft, group, basis).
      */
      Domain<D,T>* domainPtr_;

      /**
      * Pointer to %Interaction (excess free energy model).
      */
      Interaction* interactionPtr_;

      /**
      * Pointer to SCFT property calculator.
      */
      typename T::ScftThermo* scftPtr_;

      /**
      * Pointer to an %Environment.
      */
      Environment* environmentPtr_;

      /**
      * Pointer to an %Environment factory object.
      */
      typename T::EnvironmentFactory* environmentFactoryPtr_;

      /**
      * Pointer to an SCFT Iterator.
      */
      typename T::Iterator* iteratorPtr_;

      /**
      * Pointer to an Iterator factory object.
      */
      typename T::IteratorFactory* iteratorFactoryPtr_;

      /**
      * Pointer to an SCFT Sweep object.
      */
      typename T::Sweep* sweepPtr_;

      /**
      * Pointer to a sweep factory object.
      */
      typename T::SweepFactory* sweepFactoryPtr_;

      /**
      * Pointer to a Simulator.
      */
      typename T::Simulator* simulatorPtr_;

      /**
      * Pointer to a simulator factory object.
      */
      typename T::SimulatorFactory* simulatorFactoryPtr_;

      /**
      * Filemaster (holds path prefixes for input and output files).
      */
      FileMaster* fileMasterPtr_;

      /**
      * Pointer to mutable unit cell (work space).
      */
      UnitCell<D>* tmpUnitCellPtr_;

      // Boolean and enum state variables

      /**
      * Polymer model enumeration (thread or bead), read from file.
      */
      PolymerModel::Type polymerModel_;

      /**
      * Has memory been allocated for fields in FFT grid formats?
      */
      bool isAllocatedGrid_;

      /**
      * Has memory been allocated for fields in basis format?
      */
      bool isAllocatedBasis_;

      /**
      * Has the mixture been initialized?
      */
      bool hasMixture_;

      // Private member functions

      /**
      * Get the Mixture by non-const reference (private).
      */
      typename T::Mixture & mixture_();

      /**
      * Get the Domain by non-const reference (private).
      */
      Domain<D,T>& domain_();

      /**
      * Get the concentration (c) fields by non-const reference (private).
      */
      CFields<D,T>& c_();

      /**
      * Allocate memory for fields in grid formats (private).
      */
      void allocateFieldsGrid();

      /**
      * Allocate memory for fields in basis format (private).
      */
      void allocateFieldsBasis();

      /**
      * Read a string and echo to log file.
      *
      * Used to read filenames in readCommands.
      *
      * \param in  input stream (i.e., input file)
      * \param string  string to read and echo
      */
      void readEcho(std::istream& in, std::string& string) const;

      /**
      * Read a floating point number and echo to log file.
      *
      * Used to read numerical values in readCommands.
      *
      * \param in  input stream (i.e., input file)
      * \param value  number to read and echo
      */
      void readEcho(std::istream& in, double& value) const;

   };

   // Inline member functions

   // Get the Mixture (const).
   template <int D, class T> inline
   typename T::Mixture const & System<D,T>::mixture() const
   {
      UTIL_ASSERT(mixturePtr_);
      return *mixturePtr_;
   }

   // Get the MixtureModifier (non-const).
   template <int D, class T> inline
   typename T::MixtureModifier& System<D,T>::mixtureModifier()
   {
      UTIL_ASSERT(mixtureModifierPtr_);
      return *mixtureModifierPtr_;
   }

   // Get the %Interaction (const).
   template <int D, class T> inline
   Interaction const & System<D,T>::interaction() const
   {
      UTIL_ASSERT(interactionPtr_);
      return *interactionPtr_;
   }

   // Get the %Interaction (non-const).
   template <int D, class T> inline
   Interaction& System<D,T>::interaction()
   {
      UTIL_ASSERT(interactionPtr_);
      return *interactionPtr_;
   }

   // Get the Domain (const).
   template <int D, class T> inline
   Domain<D,T> const & System<D,T>::domain() const
   {
      UTIL_ASSERT(domainPtr_);
      return *domainPtr_;
   }

   // Get the WaveList (non-const).
   template <int D, class T> inline
   typename T::WaveList& System<D,T>::waveList()
   {
      UTIL_ASSERT(domainPtr_);
      return domainPtr_->waveList();
   }

   // Get the FileMaster (const).
   template <int D, class T> inline
   FileMaster const & System<D,T>::fileMaster() const
   {
      UTIL_ASSERT(fileMasterPtr_);
      return *fileMasterPtr_;
   }

   // Get the FileMaster (non-const).
   template <int D, class T> inline
   FileMaster& System<D,T>::fileMaster()
   {
      UTIL_ASSERT(fileMasterPtr_);
      return *fileMasterPtr_;
   }

   // Accessors for field containers

   // Get the container of c fields (const).
   template <int D, class T> inline
   CFields<D,T> const & System<D,T>::c() const
   {  return *cPtr_; }

   // Get the container of w fields (const).
   template <int D, class T> inline
   WFields<D,T> const & System<D,T>::w() const
   {  return *wPtr_; }

   // Get the container of w fields (non-const).
   template <int D, class T> inline
   WFields<D,T>& System<D,T>::w()
   {  return *wPtr_; }

   // Get the container of external fields (const).
   template <int D, class T> inline
   WFields<D,T> const & System<D,T>::h() const
   {  return *hPtr_; }

   // Get the container of external fields (non-const).
   template <int D, class T> inline
   WFields<D,T>& System<D,T>::h()
   {  return *hPtr_; }

   // Get the mask field (const).
   template <int D, class T> inline
   Mask<D,T> const & System<D,T>::mask() const
   {  return *maskPtr_; }

   // Get the mask field (non-const).
   template <int D, class T> inline
   Mask<D,T>& System<D,T>::mask()
   {  return *maskPtr_; }

   // Accessors for optional elements

   // Does this system have an %Environment?
   template <int D, class T> inline
   bool System<D,T>::hasEnvironment() const
   {  return (environmentPtr_); }

   // Get the %Environment (const).
   template <int D, class T> inline
   Environment const & System<D,T>::environment() const
   {
      UTIL_ASSERT(environmentPtr_);
      return *environmentPtr_;
   }

   // Get the %Environment (non-const).
   template <int D, class T> inline
   Environment & System<D,T>::environment()
   {
      UTIL_ASSERT(environmentPtr_);
      return *environmentPtr_;
   }

   // Get the Scft thermodynamics calculator (const).
   template <int D, class T> inline
   typename T::ScftThermo const & System<D,T>::scft() const
   {
      UTIL_ASSERT(scftPtr_);
      return *scftPtr_;
   }

   // Get the Scft thermodynamics calculator (non-const).
   template <int D, class T> inline
   typename T::ScftThermo & System<D,T>::scft()
   {
      UTIL_ASSERT(scftPtr_);
      return *scftPtr_;
   }

   // Does this system have an SCFT Iterator?
   template <int D, class T> inline
   bool System<D,T>::hasIterator() const
   {  return (iteratorPtr_); }

   // Get the SCFT Iterator (const).
   template <int D, class T> inline
   typename T::Iterator const & System<D,T>::iterator() const
   {
      UTIL_ASSERT(iteratorPtr_);
      return *iteratorPtr_;
   }

   // Get the SCFT Iterator (non-const).
   template <int D, class T> inline
   typename T::Iterator& System<D,T>::iterator()
   {
      UTIL_ASSERT(iteratorPtr_);
      return *iteratorPtr_;
   }

   // Does this system have an SCFT Sweep?
   template <int D, class T> inline
   bool System<D,T>::hasSweep() const
   {  return (sweepPtr_); }

   // Does this system have a Simulator?
   template <int D, class T> inline
   bool System<D,T>::hasSimulator() const
   {  return (simulatorPtr_); }

   // Get the Simulator (const).
   template <int D, class T> inline
   typename T::Simulator const & System<D,T>::simulator() const
   {
      UTIL_ASSERT(simulatorPtr_);
      return *simulatorPtr_;
   }

   // Get the Simulator (non-const).
   template <int D, class T> inline
   typename T::Simulator& System<D,T>::simulator()
   {
      UTIL_ASSERT(simulatorPtr_);
      return *simulatorPtr_;
   }

   // Private inline functions:

   // Get the Mixture (non-const).
   template <int D, class T> inline
   typename T::Mixture & System<D,T>::mixture_()
   {  return *mixturePtr_; }

   // Get the Domain (non-const).
   template <int D, class T> inline
   Domain<D,T> & System<D,T>::domain_()
   {  return *domainPtr_; }

   // Get the CFields container (non-const).
   template <int D, class T> inline
   CFields<D,T> & System<D,T>::c_()
   {  return *cPtr_; }

} // namespace Rp
} // namespace Pscf
#endif

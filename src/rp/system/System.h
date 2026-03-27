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
   namespace Prdc {
      template <int D> class UnitCell;
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
   * <b> Usage </b>: A specialization of Rp::System\<D, T\> is used as a
   * base class for each System\<D\> class defined in namespaces Rpc and 
   * Rpg, for D=1, 2, or 3.  In this use, template parameter T is taken 
   * to be an instance of a class template named Types that is defined in 
   * each of these two program-level namespaces. For example, in the Rpc 
   * namespace, for each value of D, class Rpc::System\<D\> is derived 
   * from the class Prdc::System\< D, Rpc::Types\<D\> >. For each such 
   * specialization, class Types\<D\> defines a set of typename aliases 
   * for classes used in the relevant program-level namespace, for the 
   * specified value of D.  For example, for each value of D, the typename 
   * Rpc::Types\<D\>::Mixture is an alias for the type Rpc::Mixture<D>
   * used to represent a mixture in the Rpc namespace for systems of 
   * dimension D. See the definitions of Rpc::Types and Rpg::Types (files
   * src/rpc/system/Types.h and src/rpg/system/Types.h) for lists of all 
   * of the typenames defined in each these class templates.
   *
   * In the remainder of this documentation for the Rp::System template, 
   * unqualified names such as "Mixture", "Iterator", etc. are used 
   * as shorthand for typename aliases such as T::Mixture or T::Iterator 
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
   * to describe systems that are subjected to inhomgeneous imposed 
   * environments (such as in thin films) and are otherwise left empty 
   * and unused. The bool h().hasData() and mask().hasData() functions
   * may be queried to determine if these components are in use.
   *
   * A System may also optionally have:
   *
   *    - an %Environment
   *    - an Iterator
   *    - a Sweep,
   *    - a Simulator
   *
   * Each optional component is constructed if and only if the parameter
   * file contains a corresponding optional parameter file block. The
   * %Environment is only used to generate external and mask fields to
   * describe inhomogeneous environments, and is omitted in standard 
   * calculations for structures formed in a homogeneous environment. 
   * An Iterator is only used for SCFT calculations. A Sweep is only
   * used for "sweep" calculations that solve a sequence of SCFT problems
   * along a path through parameter space.  A Simulator is only used for 
   * PS-FTS calculations, for, i.e., field theoretic simulations based on 
   * a partial saddle-point approximation. The Iterator and Sweep objects
   * may thus be omitted for PS-FTS calculations, while the Simulator
   * object may be omitted for SCFT calculations. The hasEnvironment(), 
   * hasIterator(), hasSweep(), and hasSimulator() member functions may 
   * be queried after processing of the parameter file to determine 
   * which of these optional components have been constructed.
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

      // Public type name aliases
      using MixtureT = typename T::Mixture;
      using InteractionT = typename T::Interaction;
      using DomainT = typename T::Domain;
      using WFieldsT = typename T::WFields;
      using CFieldsT = typename T::CFields;
      using MaskT = typename T::Mask;
      using RFieldT = typename T::RField;

      /// \name Construction and Destruction
      ///@{

      /**
      * Constructor.
      *
      * When a specialization of System<D,T> is used as a base class for 
      * a subclass defined in the Rpc or Rpg program-level namespace,
      * such as Rpc::System\<D\>, the typename T::System is an alias for 
      * the name of the subclass defined in Rpc or Rpg. In the constructor
      * for such a subclass, the associated instance of the subclass must
      * be passed to the Rp::System<D,T> base class constructor as *this,
      * and the address of this T::System subclass instance is retained in 
      * the Rp::System base class instance by a private member variable 
      * named systemPtr_ that is of type T::System*.  See definitions of 
      * the constructors for the Rpc::System and Rpc::System class 
      * templates for examples of this usage. 
      *
      * \param system  instance of System subclass
      */
      System(typename T::System& system);

      /**
      * Destructor.
      */
      ~System();

      // Suppress several compiler-generated member functions
      System() = delete;
      System(System<D,T> const &) = delete;
      System<D,T>& operator = (System<D,T> const & ) = delete;

      ///@}
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
      * will all return false. If the system has an %Environment,
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
      * domain().wavelist().clearUnitCellData(), clearCFields(), and, if
      * an %Environment exists, environment().reset().
      */
      void clearUnitCellData();

      ///@}
      /// \name Field Containers
      ///@{

      /**
      * Get the monomer concentration (c) fields (const).
      */
      typename T::CFields const & c() const;

      /**
      * Get the chemical potential (w) fields (non-const).
      */
      typename T::WFields& w();

      /**
      * Get the chemical potential (w) fields (const).
      */
      typename T::WFields const & w() const;

      /**
      * Get the external potential (h) fields (non-const).
      */
      typename T::WFields& h();

      /**
      * Get the external potential (h) fields (const).
      */
      typename T::WFields const & h() const;

      /**
      * Get the mask (non-const).
      */
      typename T::Mask& mask();

      /**
      * Get the mask (const).
      */
      typename T::Mask const & mask() const;

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
      * Get the %Interaction (non-const).
      */
      typename T::Interaction& interaction();

      /**
      * Get the %Interaction (const).
      */
      typename T::Interaction const & interaction() const;

      /**
      * Get the Domain (const).
      */
      typename T::Domain const & domain() const;

      /**
      * Does this system have an %Environment?
      */
      bool hasEnvironment() const;

      /**
      * Get the %Environment (non-const).
      */
      typename T::Environment& environment();

      /**
      * Get the %Environment (const).
      */
      typename T::Environment const & environment() const;

      /**
      * Get the ScftThermo object (non-const).
      */
      typename T::ScftThermo& scft();

      /**
      * Get the ScftThermo object (const).
      */
      typename T::ScftThermo const & scft() const;

      /**
      * Does this system have an Iterator?
      */
      bool hasIterator() const;

      /**
      * Get the Iterator (non-const).
      */
      typename T::Iterator& iterator();

      /**
      * Get the Iterator (const).
      */
      typename T::Iterator const & iterator() const;

      /**
      * Does this system have a Sweep?
      */
      bool hasSweep() const;

      /**
      * Does this system have a Simulator?
      */
      bool hasSimulator() const;

      /**
      * Get the Simulator (non-const).
      */
      typename T::Simulator& simulator();

      /**
      * Get the Simulator (const).
      */
      typename T::Simulator const & simulator() const;

      /**
      * Get the FileMaster (non-const).
      *
      * Access (non-const) is used in some unit tests.
      */
      FileMaster& fileMaster();

      /**
      * Get the FileMaster (const).
      */
      FileMaster const & fileMaster() const;

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
      * Clear timers
      */
      void clearTimers();

      ///@}

   protected:

      /**
      * Set the number of threads given as a command line argument.
      *
      * This function is called in the setOpts function that processes
      * command line arguments. The argument nThread may be passed to the 
      * main program as the argument of the -t option. This value gives
      * the number of threads in a threaded CPU implementation or an
      * explicit choice for the maximum number of threads per block in 
      * GPU code. 
      *
      * The do-nothing default implementation is used by CPU code that
      * has not implemented threading (the current status). 
      *
      * \param nThread  thread count
      */
      virtual void setThreadCount(int nThread)
      {};

   private:

      // Component objects

      /**
      * Chemical potential fields.
      */
      typename T::WFields w_;

      /**
      * Monomer concentration / volume fraction fields.
      */
      typename T::CFields c_;

      /**
      * External potential fields.
      */
      typename T::WFields h_;

      /**
      * Field to which the total density is constrained.
      */
      typename T::Mask mask_;

      // Pointers to associated objects

      /**
      * Pointer to enclosing instance of System subclass.
      */
      typename T::System* systemPtr_;

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
      typename T::Domain* domainPtr_;

      /**
      * Pointer to %Interaction (excess free energy model).
      */
      typename T::Interaction* interactionPtr_;

      /**
      * Pointer to an %Environment.
      */
      typename T::Environment* environmentPtr_;

      /**
      * Pointer to an %Environment factory object.
      */
      typename T::EnvironmentFactory* environmentFactoryPtr_;

      /**
      * Pointer to SCFT property calculator.
      */
      typename T::ScftThermo* scftPtr_;

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

      /**
      * Polymer model enumeration (thread or bead), read from file.
      */
      PolymerModel::Type polymerModel_;

      // Boolean state variables

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

      #if 0
      // Mutable work space

      mutable UnitCell<D> tmpUnitCell_;
      #endif

      // Private member functions

      /**
      * Get the Mixture by non-const reference (private interface).
      */
      typename T::Mixture & mixture_();

      /**
      * Get the Domain by non-const reference (private interface).
      */
      typename T::Domain& domain_();

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
   {  return *mixturePtr_; }

   // Get the MixtureModifier (non-const).
   template <int D, class T> inline 
   typename T::MixtureModifier& System<D,T>::mixtureModifier()
   {
      UTIL_ASSERT(mixtureModifierPtr_);
      return *mixtureModifierPtr_;
   }

   // Get the %Interaction (non-const).
   template <int D, class T> inline 
   typename T::Interaction& System<D,T>::interaction()
   {
      UTIL_ASSERT(interactionPtr_);
      return *interactionPtr_;
   }

   // Get the %Interaction (const).
   template <int D, class T> inline 
   typename T::Interaction const & System<D,T>::interaction() const
   {
      UTIL_ASSERT(interactionPtr_);
      return *interactionPtr_;
   }

   // Get the Domain (const).
   template <int D, class T> inline 
   typename T::Domain const & System<D,T>::domain() const
   {  return *domainPtr_; }

   // Does this system have an %Environment?
   template <int D, class T> inline 
   bool System<D,T>::hasEnvironment() const
   {  return (environmentPtr_); }

   // Get the %Environment (non-const).
   template <int D, class T> inline 
   typename T::Environment & System<D,T>::environment()
   {
      UTIL_ASSERT(environmentPtr_);
      return *environmentPtr_;
   }

   // Get the %Environment (const).
   template <int D, class T> inline 
   typename T::Environment const & System<D,T>::environment() const
   {
      UTIL_ASSERT(environmentPtr_);
      return *environmentPtr_;
   }

   // Get the Scft calculator (non-const).
   template <int D, class T> inline 
   typename T::ScftThermo & System<D,T>::scft()
   {
      UTIL_ASSERT(scftPtr_);
      return *scftPtr_;
   }

   // Get the Scft calculator (const).
   template <int D, class T> inline 
   typename T::ScftThermo const & System<D,T>::scft() const
   {
      UTIL_ASSERT(scftPtr_);
      return *scftPtr_;
   }

   // Does this system have an Iterator?
   template <int D, class T> inline 
   bool System<D,T>::hasIterator() const
   {  return (iteratorPtr_); }

   // Get the Iterator (non-const).
   template <int D, class T> inline 
   typename T::Iterator& System<D,T>::iterator()
   {
      UTIL_ASSERT(iteratorPtr_);
      return *iteratorPtr_;
   }

   // Get the Iterator (const).
   template <int D, class T> inline 
   typename T::Iterator const & System<D,T>::iterator() const
   {
      UTIL_ASSERT(iteratorPtr_);
      return *iteratorPtr_;
   }

   // Does this system have a Sweep?
   template <int D, class T> inline 
   bool System<D,T>::hasSweep() const
   {  return (sweepPtr_); }

   // Does this system have a Simulator?
   template <int D, class T> inline 
   bool System<D,T>::hasSimulator() const
   {  return (simulatorPtr_); }

   // Get the Simulator (non-const).
   template <int D, class T> inline 
   typename T::Simulator& System<D,T>::simulator()
   {
      UTIL_ASSERT(simulatorPtr_);
      return *simulatorPtr_;
   }

   // Get the Simulator (const).
   template <int D, class T> inline 
   typename T::Simulator const & System<D,T>::simulator() const
   {
      UTIL_ASSERT(simulatorPtr_);
      return *simulatorPtr_;
   }

   // Get the FileMaster (non-const).
   template <int D, class T> inline 
   FileMaster& System<D,T>::fileMaster()
   {  return *fileMasterPtr_; }

   // Get the FileMaster (const).
   template <int D, class T> inline 
   FileMaster const & System<D,T>::fileMaster() const
   {  return *fileMasterPtr_; }

   // Get the container of c fields (const).
   template <int D, class T> inline
   typename T::CFields const & System<D,T>::c() const
   {  return c_; }

   // Get the container of w fields (non-const).
   template <int D, class T> inline
   typename T::WFields& System<D,T>::w()
   {  return w_; }

   // Get the container of w fields (const).
   template <int D, class T> inline
   typename T::WFields const & System<D,T>::w() const
   {  return w_; }

   // Get the container of external fields (non-const).
   template <int D, class T> inline 
   typename T::WFields& System<D,T>::h()
   {  return h_; }

   // Get the container of external fields (const).
   template <int D, class T> inline 
   typename T::WFields const & System<D,T>::h() const
   {  return h_; }

   // Get the mask field (non-const).
   template <int D, class T> inline 
   typename T::Mask& System<D,T>::mask()
   {  return mask_; }

   // Get the mask field (const).
   template <int D, class T> inline 
   typename T::Mask const & System<D,T>::mask() const
   {  return mask_; }

   // Private inline functions:

   // Get the Mixture (non-const).
   template <int D, class T> inline 
   typename T::Mixture & System<D,T>::mixture_() 
   {  return *mixturePtr_; }

   // Get the Domain (non-const).
   template <int D, class T> inline 
   typename T::Domain & System<D,T>::domain_() 
   {  return *domainPtr_; }

} // namespace Rp
} // namespace Pscf
#endif

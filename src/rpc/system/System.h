#ifndef RPC_SYSTEM_H
#define RPC_SYSTEM_H

/*
* PSCF - Polymer Self-Consistent Field Theory
*
* Copyright 2015 - 2025, The Regents of the University of Minnesota
* Distributed under the terms of the GNU General Public License.
*/

// Header file includes
#include <util/param/ParamComposite.h>   // base class
#include <rpc/solvers/Mixture.h>         // member
#include <rpc/field/Domain.h>            // member
#include <rpc/field/WFieldContainer.h>   // member
#include <rpc/field/CFieldContainer.h>   // member
#include <rpc/field/Mask.h>              // member
#include <prdc/crystal/UnitCell.h>       // member
#include <pscf/chem/PolymerModel.h>      // member
#include <util/misc/FileMaster.h>        // member

// Forward declarations
namespace Util {
   template <typename T, int N> class FSArray;
}
namespace Pscf {
   class Interaction;
   namespace Prdc {
      class Environment;
      namespace Cpu {
         template <int D> class RField;
      }
   }
   namespace Rpc {
      template <int D> class EnvironmentFactory;
      template <int D> class ScftThermo;
      template <int D> class Iterator;
      template <int D> class IteratorFactory;
      template <int D> class Sweep;
      template <int D> class SweepFactory;
      template <int D> class Simulator;
      template <int D> class SimulatorFactory;
   }
}

namespace Pscf {
namespace Rpc {

   // Namespaces that may be used implicitly 
   using namespace Util;
   using namespace Prdc;
   using namespace Prdc::Cpu;

   /**
   * Main class, representing one complete system.
   *
   * A System has (among other components):
   *
   *    - a Mixture (container for polymer and solvent solvers)
   *    - an Interaction (list of binary interaction parameters)
   *    - a Domain (description of unit cell and discretization)
   *    - a container of monomer chemical potential (w) fields 
   *    - a container of monomer concentration fields (c) fields
   *
   * A System may also optionally have Environment Iterator, Sweep, and 
   * Simulator (BdSimulator or McSimulator) components. Iterator and Sweep
   * objects are only used for SCFT calculations. A Simulator object is 
   * only used for PS-FTS calculations (i.e., field theoretic simulations
   * that use a partial saddle-point approximation).
   *
   * System is a class template with an integer template parameter
   * D = 1, 2, or 3 that represents the dimension of space. Many related
   * components are also class templates of the same type. Names such as
   * System, Mixture, Domain, etc. mentioned above are thus names of
   * class templates, while actual class names are of the form
   * Mixture\<D\>, Domain\<D\>, etc. with D=1, 2, or 3.
   *
   * Usage of a System\<D\> object in the pscf_pc main program looks
   * something like this:
   * \code
   *    System<D> system;
   *    system.setOptions(argc, argv);
   *    system.readParam();
   *    system.readCommands();
   * \endcode
   * where argc, and argv are parameters containing information about
   * command line arguments that must be passed from the main program.
   * This is implemented as function template Pscf::Rpc::run in the
   * file src/rpc/pscf_pc.cpp.
   *
   * See also:
   * <ul>
   *  <li> \ref scft_param_pc_page   "Parameter File: SCFT" </li>
   *  <li> \ref psfts_param_page     "Parameter File: PS-FTS" </li>
   *  <li> \ref rpc_System_page      "Parameter File: Full Format" </li>
   *  <li> \ref scft_command_pc_page "Command File Format" </li>
   * </ul>
   *
   * \ingroup Pscf_Rpc_Module
   */
   template <int D>
   class System : public ParamComposite
   {

   public:

      // Public type name aliases
      using MixtureT = Mixture<D>;
      using InteractionT = Interaction;
      using DomainT = Domain<D>;
      using WFieldContainerT = WFieldContainer<D>;
      using CFieldContainerT = CFieldContainer<D>;
      using MaskT = Mask<D>;
      using FieldT = RField<D>;

      /// \name Construction and Destruction
      ///@{

      /**
      * Constructor.
      */
      System();

      /**
      * Destructor.
      */
      ~System();

      ///@}
      /// \name Lifetime (Actions in Main Program)
      ///@{

      /**
      * Process command line options.
      *
      * This function takes the same arguments as any C/C++ main program
      * function. The arguments of the main function should d be passed
      * to this function unaltered, to allow this function to process the
      * command line options.
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
      /// \name Unit Cell Modifiers
      ///@{

      /**
      * Set parameters of the associated unit cell.
      *
      * The lattice (i.e., lattice system type) set in the UnitCell<D>
      * unitCell input parameter must agree with any lattice enum value
      * that was set previously in the parameter file.
      *
      * If a space group has been set but a basis has not yet been
      * constructed, then this and the other setUnitCell member function
      * will initialize a symmetry-adpated basis as a side effect.
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
      * See documentation of setUnitCell(UnitCell<D> const &) regarding
      * possible construction of a basis as a side effect.
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
      * an Environment exists, environment().reset(). 
      */
      void clearUnitCellData();

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
      * If argument needStress is true, then this function also calls
      * computeStress() to compute the stress.
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
      * by the Mixture class. If an Environment exists, it also allows the
      * Environment to compute a modified definition of the stress.
      *
      * \pre w().hasData() == true
      * \pre c().hasData() == true
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
      * \pre w().isSymmetric() == true if iterator is symmetric
      * \post c().hasData() == true
      * \post scft().hasData() == true on successful convergence
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
      * along a contour in parameter space. The nature of this sequence
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
      * \pre Function hasSimulator() must return true
      * \pre Function w().hasData() must return true
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
      * all return false; if the system has an Environment, 
      * environment().needsUpdate() will return true.
      */
      void clearCFields();

      ///@}
      /// \name Property Output
      ///@{

      /**
      * Write partial parameter file to an ostream.
      *
      * This function writes the Mixture, Interaction, and Domain blocks
      * of a parameter file, as well as any Environment and Iterator 
      * blocks, but omits any Sweep or Simulator blocks. The intent is 
      * to produce an output during an SCFT sweep that only refers to 
      * parameters relevant to a single state point, in a form that could 
      * be used as a parameter file for a single SCFT calculation.
      *
      * \param out  output stream
      */
      void writeParamNoSweep(std::ostream& out) const;

      ///@}
      /// \name Field Containers
      ///@{

      /**
      * Get monomer concentration (c) field container (const reference).
      */
      CFieldContainer<D> const & c() const;

      /**
      * Get chemical potential (w) field container (non-const reference).
      */
      WFieldContainer<D>& w();

      /**
      * Get chemical potential (w) field container (const reference).
      */
      WFieldContainer<D> const & w() const;

      /**
      * Get external potential field container (non-const reference).
      */
      WFieldContainer<D>& h();

      /**
      * Get external potential (h) field container (const reference).
      */
      WFieldContainer<D> const & h() const;

      /**
      * Get the mask field container (non-const reference).
      */
      Mask<D>& mask();

      /**
      * Get the mask field container (const reference).
      */
      Mask<D> const & mask() const;

      ///@}
      /// \name Member Object Accessors
      ///@{

      /**
      * Get the Mixture (non-const reference).
      */
      Mixture<D>& mixture();

      /**
      * Get the Mixture (const reference).
      */
      Mixture<D> const & mixture() const;

      /**
      * Get the Interaction (non-const reference).
      */
      Interaction& interaction();

      /**
      * Get the Interaction (const reference).
      */
      Interaction const & interaction() const;

      /**
      * Get the Domain (const reference).
      */
      Domain<D> const & domain() const;

      /**
      * Does this system have an Environment?
      */
      bool hasEnvironment() const;

      /**
      * Get the Environment (non-const reference).
      */
      Environment& environment();

      /**
      * Get the Environment (const reference).
      */
      Environment const & environment() const;

      /**
      * Get the ScftThermo<D> object (non-const reference).
      */
      ScftThermo<D>& scft();

      /**
      * Get the ScftThermo<D> object (const reference).
      */
      ScftThermo<D> const & scft() const;

      /**
      * Does this system have an Iterator?
      */
      bool hasIterator() const;

      /**
      * Get the Iterator (non-const reference).
      */
      Iterator<D>& iterator();

      /**
      * Get the Iterator (const reference).
      */
      Iterator<D> const & iterator() const;

      /**
      * Does this system have a Sweep?
      */
      bool hasSweep() const;

      /**
      * Does this system have a Simulator?
      */
      bool hasSimulator() const;

      /**
      * Get the Simulator (non-const reference).
      */
      Simulator<D>& simulator();

      /**
      * Get the Simulator (const reference).
      */
      Simulator<D> const & simulator() const;

      /**
      * Get the FileMaster (non-const reference).
      *
      * Access (non-const reference) is used in some unit tests.
      */
      FileMaster& fileMaster();

      /**
      * Get the FileMaster (const reference).
      */
      FileMaster const & fileMaster() const;

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

   private:

      // Component objects

      /**
      * Mixture object (solves MDE for all species).
      */
      Mixture<D> mixture_;

      /**
      * Domain object (unit cell, mesh, fft, space group, and basis).
      */
      Domain<D> domain_;

      /**
      * Filemaster (holds path prefixes for input and output files).
      */
      FileMaster fileMaster_;

      /**
      * Chemical potential fields.
      */
      WFieldContainer<D> w_;

      /**
      * Monomer concentration / volume fraction fields.
      */
      CFieldContainer<D> c_;

      /**
      * External potential fields.
      */
      WFieldContainer<D> h_;

      /**
      * Field to which the total density is constrained.
      */
      Mask<D> mask_;

      // Pointers to dynamic objects owned by this System

      /**
      * Pointer to Interaction (excess free energy model).
      */
      Interaction* interactionPtr_;

      /**
      * Pointer to an Environment.
      */
      Environment* environmentPtr_;

      /**
      * Pointer to an Environment factory object.
      */
      EnvironmentFactory<D>* environmentFactoryPtr_;

      /**
      * Pointer to SCFT property calculator.
      */
      ScftThermo<D>* scftPtr_;

      /**
      * Pointer to an SCFT iterator.
      */
      Iterator<D>* iteratorPtr_;

      /**
      * Pointer to an iterator factory object.
      */
      IteratorFactory<D>* iteratorFactoryPtr_;

      /**
      * Pointer to an SCFT Sweep object.
      */
      Sweep<D>* sweepPtr_;

      /**
      * Pointer to a sweep factory object.
      */
      SweepFactory<D>* sweepFactoryPtr_;

      /**
      * Pointer to a Simulator.
      */
      Simulator<D>* simulatorPtr_;

      /**
      * Pointer to a simulator factory object.
      */
      SimulatorFactory<D>* simulatorFactoryPtr_;

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

      // Mutable work space

      mutable UnitCell<D> tmpUnitCell_;

      // Private member functions

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

   // Get the Mixture by non-const reference.
   template <int D>
   inline Mixture<D>& System<D>::mixture()
   {  return mixture_; }

   // Get the Mixture by const reference.
   template <int D>
   inline Mixture<D> const & System<D>::mixture() const
   {  return mixture_; }

   // Get the Interaction (excess free energy) by non-const reference.
   template <int D>
   inline Interaction& System<D>::interaction()
   {
      UTIL_ASSERT(interactionPtr_);
      return *interactionPtr_;
   }

   // Get the Interaction (excess free energy) by const reference.
   template <int D>
   inline Interaction const & System<D>::interaction() const
   {
      UTIL_ASSERT(interactionPtr_);
      return *interactionPtr_;
   }

   // Get the Domain by const reference.
   template <int D>
   inline Domain<D> const & System<D>::domain() const
   {  return domain_; }

   // Does this system have an Environment?
   template <int D>
   inline bool System<D>::hasEnvironment() const
   {  return (environmentPtr_); }

   // Get the Environment by non-const reference.
   template <int D>
   inline Environment & System<D>::environment()
   {
      UTIL_ASSERT(environmentPtr_);
      return *environmentPtr_;
   }

   // Get the Environment by const reference.
   template <int D>
   inline Environment const & System<D>::environment() const
   {
      UTIL_ASSERT(environmentPtr_);
      return *environmentPtr_;
   }

   // Get the Scft calculator by non-const reference.
   template <int D>
   inline ScftThermo<D> & System<D>::scft()
   {
      UTIL_ASSERT(scftPtr_);
      return *scftPtr_;
   }

   // Get the Scft calculator by const reference.
   template <int D>
   inline ScftThermo<D> const & System<D>::scft() const
   {
      UTIL_ASSERT(scftPtr_);
      return *scftPtr_;
   }

   // Does this system have an Iterator?
   template <int D>
   inline bool System<D>::hasIterator() const
   {  return (iteratorPtr_); }

   // Get the Iterator by non-const reference.
   template <int D>
   inline Iterator<D>& System<D>::iterator()
   {
      UTIL_ASSERT(iteratorPtr_);
      return *iteratorPtr_;
   }

   // Get the Iterator by const reference.
   template <int D>
   inline Iterator<D> const & System<D>::iterator() const
   {
      UTIL_ASSERT(iteratorPtr_);
      return *iteratorPtr_;
   }

   // Does this system have a Sweep?
   template <int D>
   inline bool System<D>::hasSweep() const
   {  return (sweepPtr_); }

   // Does this system have a Simulator?
   template <int D>
   inline bool System<D>::hasSimulator() const
   {  return (simulatorPtr_); }

   // Get the Simulator by non-const reference.
   template <int D>
   inline Simulator<D>& System<D>::simulator()
   {
      UTIL_ASSERT(simulatorPtr_);
      return *simulatorPtr_;
   }

   // Get the Simulator by const reference.
   template <int D>
   inline Simulator<D> const & System<D>::simulator() const
   {
      UTIL_ASSERT(simulatorPtr_);
      return *simulatorPtr_;
   }

   // Get the FileMaster by non-const reference.
   template <int D>
   inline FileMaster& System<D>::fileMaster()
   {  return fileMaster_; }

   // Get the FileMaster by const reference.
   template <int D>
   inline FileMaster const & System<D>::fileMaster() const
   {  return fileMaster_; }

   // Get the container of w fields (non-const reference).
   template <int D>
   inline
   WFieldContainer<D>& System<D>::w()
   {  return w_; }

   // Get the container of w fields (const reference).
   template <int D>
   inline
   WFieldContainer<D> const & System<D>::w() const
   {  return w_; }

   // Get the container of c fields (const reference).
   template <int D>
   inline
   CFieldContainer<D> const & System<D>::c() const
   {  return c_; }

   // Get the container of external fields (non-const reference).
   template <int D>
   inline WFieldContainer<D>& System<D>::h()
   {  return h_; }

   // Get the container of external fields (const reference).
   template <int D>
   inline WFieldContainer<D> const & System<D>::h() const
   {  return h_; }

   // Get the mask field (non-const reference).
   template <int D>
   inline Mask<D>& System<D>::mask()
   {  return mask_; }

   // Get the mask field (const reference).
   template <int D>
   inline Mask<D> const & System<D>::mask() const
   {  return mask_; }

   #ifndef RPC_SYSTEM_TPP
   // Suppress implicit instantiation
   extern template class System<1>;
   extern template class System<2>;
   extern template class System<3>;
   #endif

} // namespace Rpc
} // namespace Pscf
#endif

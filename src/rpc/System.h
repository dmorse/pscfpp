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
#include <pscf/chem/PolymerModel.h>      // member
#include <util/misc/FileMaster.h>        // member


// Forward references
namespace Util {
   template <typename T> class DArray;
   template <typename T, int N> class FSArray;
}
namespace Pscf {
   class Interaction;
   namespace Prdc {
      template <int D> class UnitCell;
      namespace Cpu {
         template <int D> class RField;
         template <int D> class RFieldDft;
      }
   }
   namespace Rpc {
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

   // Namespaces that are used implicitly, without name qualification
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
   *    - a container of monomer chemical potential fields (w fields)
   *    - a container of monomer concentration fields (c fields)
   *
   * A System may also optionally own Iterator, Sweep, and Simulator
   * (BdSimulator or McSimulator) components. Iterator and Sweep objects
   * are only used for SCFT calculations. A Simulator object is only used
   * for PS-FTS calculations (i.e., field theoretic simulations that use
   * a partial saddle-point approximation).
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
      /// \name W Field (Chemical Potential Field) Modifiers
      ///@{

      /**
      * Read chemical potential fields in symmetry-adapted basis format.
      *
      * This function opens and reads the file named by the "filename"
      * string parameter, which must contain chemical potential fields
      * in symmetry-adapted basis format. This function sets the system
      * w fields equal to those given in this file, by copying the
      * coefficients of the basis expansion and computing values on a
      * real-space grid (r-grid format). System unit cell parameters
      * are also set to values read from the field file header. Upon
      * exit, both w().basis() and w().rgrid() are set, w().hasData()
      * and w().isSymmetric() are true, while hasCFields() and
      * hasFreeEnergy() are false.
      *
      * If a space group has been set but a basis has not yet been
      * constructed, then this and every other member function that reads
      * unit cell parameters from file and/or sets values for unit cell
      * parameters will construct a symmetry-adapted basis and allocate
      * memory for fields stored in basis form. Member functions that
      * may construct a basis as a side effect include this function,
      * readWRGrid, estimateWfromC, and all of the overloaded setUnitCell
      * functions.
      *
      * SCFT calculations that use an iterator that preserves space group
      * symmetry must set an initial field using a function that creates
      * fields in symmetry-adapted basis form, such as this function,
      * setWBasis, or estimateWFromC.
      *
      * \param filename  name of input w field file in basis format
      */
      void readWBasis(std::string const & filename);

      /**
      * Read chemical potential fields in real space grid (r-grid) format.
      *
      * This function opens and reads the file with the name given by the
      * "filename" string, which must contain chemical potential fields
      * in real space grid (r-grid) format. The function sets values for
      * the system w fields in r-grid format. It does not attempt to set
      * field values in symmetry-adapted basis format, because it cannot
      * assume that the r-grid field exhibits the declared space group
      * symmetry. Upon exit, w().rgrid() is set and w().hasData() is
      * true, while w().isSymmetric(), hasCFields(), and hasFreeEnergy()
      * are false. System unit cell parameters are set to values read
      * from the field file header.
      *
      * Chemical potential fields for field theoretic simulations are
      * normally initialized using a function that sets the fields in
      * r-grid format, such as this function or setWRGrid.
      *
      * \param filename  name of input w field file in r-grid format
      */
      void readWRGrid(std::string const & filename);

      /**
      * Set chemical potential fields in symmetry-adapted basis format.
      *
      * This function sets values for w fields in both symmetry-adapted
      * and r-grid format by copying coefficient values provided in the
      * "fields" container that is passed as an argument, and computing
      * values on a real-space grid. Upon return, values of both
      * w().basis() and w().rgrid() are set, while w().hasData() and
      * w().isSymmetric() are true, and hasCFields() and hasFreeEnergy()
      * are false. System unit cell parameters are not modified.
      *
      * \param fields  array of new w fields in basis format
      */
      void setWBasis(DArray< DArray<double> > const & fields);

      /**
      * Set chemical potential fields in real-space (r-grid) format.
      *
      * This function set values for w fields in r-grid format, but
      * does not set components for symmetry-adapted basis format. Upon
      * return, w().rgrid() is set and w().hasData() is true, while
      * hasCFields(), hasFreeEnergy(), and w().isSymmetric() are false.
      * System unit cell parameters are not modified.
      *
      * \param fields  array of new w fields in r-grid form
      */
      void setWRGrid(DArray< RField<D> > const & fields);

      /**
      * Construct trial w fields from c fields in basis form.
      *
      * This function reads concentration fields in symmetrized basis
      * form and constructs an initial guess for corresponding chemical
      * potential fields by setting the Lagrange multiplier pressure field
      * to zero. The result is stored in the system w field container.
      *
      * Upon return, w().hasData() and w().isSymmetric() are true, while
      * hasCFields() and hasFreeEnergy() are false. System unit cell
      * parameters are set to those read from the c field file header.
      *
      * \param filename  name of input c field file (basis format)
      */
      void estimateWfromC(std::string const & filename);

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
      * constructed, then this and the other setUnitCell member functions
      * will all construct a symmetry-adapted basis and then allocate
      * memory for fields stored in a symmetry-adapted basis format.
      *
      * \param unitCell  new UnitCell<D> (i.e., new parameters)
      */
      void setUnitCell(UnitCell<D> const & unitCell);

      /**
      * Set parameters of the associated unit cell.
      *
      * The "lattice" (lattice system) enumeration parameter must agree
      * with any lattice value that was set previously in the parameter
      * file. The logical size of the "parameters" array must agree with
      * the expected number of lattice parameters for this lattice type.
      *
      * See documentation of setUnitCell(UnitCell<D> const &) regarding
      * possible construction of a basis as a side effect.
      *
      * \param lattice  lattice system
      * \param parameters  array of new unit cell parameters
      */
      void setUnitCell(typename UnitCell<D>::LatticeSystem lattice,
                       FSArray<double, 6> const & parameters);

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
      * energy or grand-canonical free energy (i.e., pressure). Upon
      * return, hasCFields() is true.
      *
      * If argument needStress is true, then this function also calls
      * Mixture<D>::computeStress() to compute the stress.
      *
      * \pre The w().hasData() flag must be true on entry.
      *
      * \param needStress  true if stress is needed, false otherwise
      */
      void compute(bool needStress = false);

      /**
      * Iteratively solve a SCFT problem.
      *
      * This function calls the iterator to solve the SCFT problem for
      * the current system parameters, using the current chemical
      * potential fields and unit cell parameters as initial guesses.
      * Upon exit, hasCFields is set true whether or not convergence
      * is obtained to within the desired tolerance. The Helmholtz free
      * energy and pressure are computed only if convergence is obtained.
      *
      * \pre Function hasIterator() must return true
      * \pre Function w().hasData() must return true
      * \pre Function w().isSymmetric() must return true if the iterator
      * uses a symmetry-adapted basis.
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
      * \pre Function hasSweep() must return true
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

      ///@}
      /// \name SCFT Thermodynamic Properties
      ///@{

      /**
      * Compute SCFT free energy density and pressure for current fields.
      *
      * This function should be called after a successful call of
      * iterate(). Resulting values are retrieved by the fHelmholtz(),
      * fIdeal(), fInter(), fExt(), and pressure() accessor functions.
      *
      * \pre w().hasData() must return true
      * \pre hasCFields() must return true
      */
      void computeFreeEnergy();

      /**
      * Get total Helmholtz free energy per monomer / kT.
      *
      * This function retrieves a value computed by computeFreeEnergy().
      */
      double fHelmholtz() const;

      /**
      * Get the ideal gas contribution to fHelmholtz.
      *
      * This function retrieves a value computed by computeFreeEnergy().
      */
      double fIdeal() const;

      /**
      * Get the interaction contribution to fHelmholtz.
      *
      * This function retrieves a value computed by computeFreeEnergy().
      */
      double fInter() const;

      /**
      * Get the external field contribution to fHelmholtz.
      *
      * This function retrieves a value computed by computeFreeEnergy().
      */
      double fExt() const;

      /**
      * Get the precomputed pressure times monomer volume / kT.
      *
      * This function retrieves a value computed by computeFreeEnergy().
      * The value is -1 times the grand-canonical free energy per monomer
      * divided by kT.
      */
      double pressure() const;

      ///@}
      /// \name SCFT Thermodynamic Data Output
      ///@{

      /**
      * Write partial parameter file to an ostream.
      *
      * This function writes the Mixture, Interaction, and Domain blocks
      * of a parameter file, as well as any Iterator block, but omits any
      * Sweep or Simulator blocks. The intent is to produce an output
      * during an SCFT sweep that only refers to parameters relevant to a
      * single state point, in a form that could be used as a parameter
      * file for a single SCFT calculation.
      *
      * \param out  output stream
      */
      void writeParamNoSweep(std::ostream& out) const;

      /**
      * Write SCFT thermodynamic properties to a file.
      *
      * This function outputs Helmholtz free energy per monomer, pressure
      * (in units of kT per monomer volume), the volume fraction and
      * chemical potential of each species, and all unit cell parameters.
      *
      * If parameter "out" is a file that already exists, this function
      * will append to the end of that file, rather than overwriting it.
      * Calling writeParamNoSweep and writeThermo in succession with the
      * same output stream will thus produce a single file containing both
      * input parameters and calculated thermodynamic properties.
      *
      * \param out  output stream
      */
      void writeThermo(std::ostream& out);

      /**
      * Write SCFT stress to a file.
      *
      * This function outputs derivatives of SCFT free energy w/ respect
      * to each unit cell parameter.
      *
      * If parameter "out" is a file that already exists, this function
      * will append to the end of that file, rather than overwriting it.
      * The log output created by an Iterator upon convergence should call
      * writeStress after writeThermo if and only if the iterator is not
      * flexible.
      *
      * \param out  output stream
      */
      void writeStress(std::ostream& out);

      ///@}
      /// \name Field Output
      ///@{

      /**
      * Write chemical potential fields in symmetry-adapted basis format.
      *
      * \param filename  name of output file
      */
      void writeWBasis(std::string const & filename) const;

      /**
      * Write chemical potential fields in real space grid (r-grid) format.
      *
      * \param filename  name of output file
      */
      void writeWRGrid(std::string const & filename) const;

      /**
      * Write concentration fields in symmetry-adapted basis format.
      *
      * \param filename  name of output file
      */
      void writeCBasis(std::string const & filename) const;

      /**
      * Write concentration fields in real space grid (r-grid) format.
      *
      * \param filename  name of output file
      */
      void writeCRGrid(std::string const & filename) const;

      /**
      * Write c fields for all blocks and solvents in r-grid format.
      *
      * Writes concentrations for all blocks of all polymers and all
      * solvent species in r-grid format. Columns associated with blocks
      * appear ordered by polymer id and then by block id, with blocks of
      * the same polymer listed sequentially, followed by columns
      * associated with solvent species ordered by solvent id.
      *
      * \param filename  name of output file
      */
      void writeBlockCRGrid(std::string const & filename) const;

      ///@}
      /// \name Propagator Output
      ///@{

      /**
      * Write one slice of a propagator at fixed s in r-grid format.
      *
      * \param filename  name of output file
      * \param polymerId  integer id of the polymer
      * \param blockId  integer id of the block within the polymer
      * \param directionId  integer id of the direction (0 or 1)
      * \param segmentId  integer integration step index
      */
      void writeQSlice(std::string const & filename,
                       int polymerId, int blockId,
                       int directionId, int segmentId)  const;

      /**
      * Write the final slice of a propagator in r-grid format.
      *
      * \param filename  name of output file
      * \param polymerId  integer id of the polymer
      * \param blockId  integer id of the block within the polymer
      * \param directionId  integer id of the direction (0 or 1)
      */
      void writeQTail(std::string const & filename, int polymerId,
                      int blockId, int directionId)  const;

      /**
      * Write the complete propagator for one block, in r-grid format.
      *
      * \param filename  name of output file
      * \param polymerId  integer id of the polymer
      * \param blockId  integer id of the block within the polymer
      * \param directionId  integer id of the direction (0 or 1)
      */
      void writeQ(std::string const & filename, int polymerId,
                  int blockId, int directionId)  const;

      /**
      * Write all propagators of all blocks, each to a separate file.
      *
      * Write all propagators for both directions for all blocks
      * of all polymers, with each propagator in a separate file.
      * The function writeQ is called internally for each propagator,
      * and is passed an automatically generated file name. The file
      * name for each propagator is given by a string of the form
      * (basename)_(ip)_(ib)_(id), where (basename) denotes the value
      * of the std::string function parameter basename, and where
      * (ip), (ib), and (id) denote the string representations of
      * a polymer indiex ip, a block index ib, and direction index id,
      * with id = 0 or 1. For example, if basename == "out/q", then
      * the file name of the propagator for direction 1 of block 2
      * of polymer 0 would be "out/q_0_2_1".
      *
      * \param basename  common prefix for output file names
      */
      void writeQAll(std::string const & basename);

      ///@}
      /// \name Timers
      ///@{

      /**
      * Write timer information to an output stream.
      *
      * \param out  output stream
      */
      void writeTimers(std::ostream& out);

      /**
      * Clear timers
      */
      void clearTimers();

      ///@}
      /// \name Field Accessors
      ///@{

      /**
      * Get container of chemical potential fields (w fields).
      */
      WFieldContainer<D> const & w() const;

      /**
      * Get container of monomer concentration fields (c fields).
      */
      CFieldContainer<D> const & c() const;

      /**
      * Get container of external potential fields (non-const reference).
      */
      WFieldContainer<D>& h();

      /**
      * Get container of external potential fields (const reference).
      */
      WFieldContainer<D> const & h() const;

      /**
      * Get the mask (field to which total density is constrained).
      */
      Mask<D>& mask();

      /**
      * Get the mask by const reference.
      */
      Mask<D> const & mask() const;

      ///@}
      /// \name Member Object Accessors
      ///@{

      /**
      * Get the Mixture by non-const reference.
      */
      Mixture<D>& mixture();

      /**
      * Get the Mixture by const reference.
      */
      Mixture<D> const & mixture() const;

      /**
      * Get the Interaction (excess free energy) by non-const reference.
      */
      Interaction& interaction();

      /**
      * Get the Interaction (excess free energy) by const reference.
      */
      Interaction const & interaction() const;

      /**
      * Get the Domain by const reference.
      */
      Domain<D> const & domain() const;

      /**
      * Get the Iterator by non-const reference.
      */
      Iterator<D>& iterator();

      /**
      * Get the Iterator by const reference.
      */
      Iterator<D> const & iterator() const;

      /**
      * Get the Simulator by non-const reference.
      */
      Simulator<D>& simulator();

      /**
      * Get the FileMaster by non-const reference.
      *
      * Access by non-const reference is used in some unit tests.
      */
      FileMaster& fileMaster();

      /**
      * Get the FileMaster by const reference.
      */
      FileMaster const & fileMaster() const;

      ///@}
      /// \name Boolean Flags
      ///@{

      /**
      * Does this system have an Iterator?
      */
      bool hasIterator() const;

      /**
      * Does this system have a Sweep?
      */
      bool hasSweep() const;

      /**
      * Does this system have a Simulator?
      */
      bool hasSimulator() const;

      /**
      * Does this system have external potential fields?
      */
      bool hasExternalFields() const;

      /**
      * Does this system have a mask (inhomogeneous density constraint)?
      */
      bool hasMask() const;

      /**
      * Are c fields current, consistent with the current w fields?
      */
      bool hasCFields() const;

      /**
      * Is the SCFT free energy current, consistent with current w fields?
      */
      bool hasFreeEnergy() const;

      /**
      * Mark c-fields and free energy as outdated or invalid.
      *
      * Upon return, hasCfields() and hasFreeEnergy() both return false.
      * This function should be called by functions that modify any of
      * inputs to the solution of the modified diffusion equation and
      * calculation of c-fields and free energy, including the w fields,
      * unit cell parameters, external fields or mask. 
      */
      void clearCFields();

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
      * Pointer to an iterator.
      */
      Iterator<D>* iteratorPtr_;

      /**
      * Pointer to an iterator factory object.
      */
      IteratorFactory<D>* iteratorFactoryPtr_;

      /**
      * Pointer to a Sweep object.
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

      // SCFT Thermodynamic properties

      /**
      * Helmholtz free energy per monomer / kT.
      */
      double fHelmholtz_;

      /**
      * Ideal gas contribution to fHelmholtz_.
      *
      * This includes the internal energy and entropy of
      * non-interacting molecules in the current w fields.
      */
      double fIdeal_;

      /**
      * Interaction contribution to fHelmholtz_.
      */
      double fInter_;

      /**
      * External field contribution to fHelmholtz_ (if any).
      */
      double fExt_;

      /**
      * Pressure times monomer volume / kT.
      *
      * This is -1 times the grand-canonical free energy per monomer,
      * divided by kT.
      */
      double pressure_;

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

      /**
      * Have c fields been computed for the current w fields?
      *
      * When hasCFields_ is true, the c fields for all individual blocks
      * and solvent species in the Mixture and the total concentration
      * fields for all monomer types in the System::c_ container are  all
      * values computed from the current w fields in System::w_ container,
      * using the current unit cell parameters. This should be set true
      * when all c fields are computed, and set false whenever the system
      * w fields, the mask or external fields, or the system unit cell
      * parameters are modified.
      */
      bool hasCFields_;

      /**
      * Has SCFT free energy been computed for the current w and c fields?
      *
      * This is set true in the computeFreeEnergy function, and is set
      * false whenever the system w fields, mask or external fields, or
      * system unit cell parameters are modified.
      */
      bool hasFreeEnergy_;

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
      * Read a field file header, make the basis if not done previously.
      *
      * Used to peek at a file header to get initial unit cell parameters
      * that can be used to initialize a basis.
      *
      * \param filename  name of field file
      */
      void readFieldHeader(std::string const & filename);

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

   // Get the Simulator by non-const reference.
   template <int D>
   inline Simulator<D>& System<D>::simulator()
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

   // Get the Helmholtz free energy per monomer / kT.
   template <int D>
   inline double System<D>::fHelmholtz() const
   {
      UTIL_CHECK(hasFreeEnergy_);
      return fHelmholtz_;
   }

   // Get the ideal gas contribution to fHelmholtz.
   template <int D>
   inline double System<D>::fIdeal() const
   {
      UTIL_CHECK(hasFreeEnergy_);
      return fIdeal_;
   }

   // Get the interaction contribution to fHelmholtz.
   template <int D>
   inline double System<D>::fInter() const
   {
      UTIL_CHECK(hasFreeEnergy_);
      return fInter_;
   }

   // Get the external field contribution to fHelmholtz.
   template <int D>
   inline double System<D>::fExt() const
   {
      UTIL_CHECK(hasFreeEnergy_);
      return fExt_;
   }

   // Get the precomputed pressure (units of kT / monomer volume).
   template <int D>
   inline double System<D>::pressure() const
   {
      UTIL_CHECK(hasFreeEnergy_);
      return pressure_;
   }

   // Does this system have an Iterator?
   template <int D>
   inline bool System<D>::hasIterator() const
   {  return (iteratorPtr_); }

   // Does this system have a Sweep?
   template <int D>
   inline bool System<D>::hasSweep() const
   {  return (sweepPtr_); }

   // Does this system have external potential fields?
   template <int D>
   inline bool System<D>::hasExternalFields() const
   {  return h_.hasData(); }

   // Does this system have a mask?
   template <int D>
   inline bool System<D>::hasMask() const
   {  return mask_.hasData(); }

   // Does this system have a Simulator?
   template <int D>
   inline bool System<D>::hasSimulator() const
   {  return (simulatorPtr_); }

   // Have the c fields been computed for the current w fields?
   template <int D>
   inline bool System<D>::hasCFields() const
   {  return hasCFields_; }

   // Has the free energy been computed for the current w fields?
   template <int D>
   inline bool System<D>::hasFreeEnergy() const
   {  return hasFreeEnergy_; }

   #ifndef RPC_SYSTEM_TPP
   // Suppress implicit instantiation
   extern template class System<1>;
   extern template class System<2>;
   extern template class System<3>;
   #endif

} // namespace Rpc
} // namespace Pscf
#endif

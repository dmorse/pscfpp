#ifndef PSPG_SYSTEM_H
#define PSPG_SYSTEM_H

/*
* PSCF - Polymer Self-Consistent Field Theory
*
* Copyright 2016 - 2022, The Regents of the University of Minnesota
* Distributed under the terms of the GNU General Public License.
*/

#include <util/param/ParamComposite.h>     // base class

#include <rpg/solvers/Mixture.h>          // member
#include <rpg/field/Domain.h>             // member
#include <rpg/field/FieldIo.h>            // member
#include <rpg/field/WFieldContainer.h>    // member
#include <rpg/field/CFieldContainer.h>    // member

#include <prdc/cuda/RField.h>              // member
#include <prdc/cuda/RFieldDft.h>           // member

#include <pscf/homogeneous/Mixture.h>      // member

#include <util/misc/FileMaster.h>          // member
#include <util/containers/DArray.h>        // member template
#include <util/containers/FSArray.h>       // member template

namespace Pscf {
   class Interaction;
   namespace Prdc {
      template <int D> class UnitCell;
   }
}

namespace Pscf {
namespace Rpg {

   template <int D> class Iterator;
   template <int D> class IteratorFactory;
   template <int D> class Sweep;
   template <int D> class SweepFactory;
   template <int D> class Compressor;
   template <int D> class CompressorFactory;
   template <int D> class Simulator;
   template <int D> class SimulatorFactory;

   using namespace Util;
   using namespace Pscf::Prdc;
   using namespace Pscf::Prdc::Cuda;

   /**
   * Main class in SCFT simulation of one system.
   *
   * A System has (among other components):
   *
   *    - a Mixture (a container for polymer and solvent solvers)
   *    - an Interaction (list of chi parameters)
   *    - a Domain (description of unit cell and discretization)
   *    - a container of monomer chemical potential fields (w fields)
   *    - a container of monomer concentration fields (c fields)
   *
   * A System may also optionally contain Iterator and Sweep objects.
   *
   * Typical usage of a System<D> object looks something like this:
   * \code
   *    System<D> system;
   *    system.setOptions(argc, argv);
   *    system.readParam();
   *    system.readCommands();
   * \endcode
   * where argc, and argv are parameters containing information about
   * command line arguments that must be passed from the main program.
   * This is implemented as function template Pscf::Rpg::run in the
   * file src/rpg/pscf_pg.cpp.
   *
   * \ref user_param_pg_page "Parameter File Format"
   * \ref user_command_pc_page "Command File Format"
   * \ingroup Pscf_Rpg_Module
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
      /// \name Lifetime (Actions)
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
      * Read body of parameter file (without opening, closing lines).
      *
      * \param in input parameter stream
      */
      virtual void readParameters(std::istream& in);

      /**
      * Read command script.
      *
      * \param in command script file.
      */
      void readCommands(std::istream& in);

      /**
      * Read commands from default command file.
      *
      * This function reads the parameter file set by the -c command
      * line option.
      */
      void readCommands();

      ///@}
      /// \name W Field Modifiers 
      ///@{

      /**
      * Read chemical potential fields in symmetry adapted basis format.
      *
      * This function opens and reads the file with the name given by the
      * "filename" string, which must contain chemical potential fields
      * in symmetry-adapted basis format. The function copies these
      * fields to set new values for the system w fields in basis format, 
      * and also computes and resets the system w fields in r-space 
      * format. On exit, both w().basis() and w().rgrid() have been 
      * reset, w().hasData and w().isSymmetric() are true, and 
      * hasCFields() is false.
      *
      * \param filename name of input w-field basis file
      */
      void readWBasis(const std::string & filename);

      /**
      * Read chemical potential fields in real space grid (r-grid) format.
      *
      * This function opens and reads the file with the name given by the
      * "filename" string, which must contains chemical potential fields
      * in real space grid (r-grid) format. The function sets values for 
      * system w fields in r-grid format. It does not set attempt to set
      * field values in symmetry-adapted basis format, because it cannot
      * be known whether the r-grid field exhibits the declared space
      * group symmetry.  On exit, w().rgrid() is reset and w().hasData() 
      * is true, while w().isSymmetric() and hasCFields() are false.
      *
      * \param filename  name of input w-field basis file
      */
      void readWRGrid(const std::string & filename);

      /**
      * Set chemical potential fields, in symmetry-adapted basis format.
      *
      * This function sets values for w fields in both symmetry adapted 
      * and r-grid format.  On exit, values of both w().basis() and 
      * w().rgrid() are reset, w().hasData() and w().isSymmetric() are 
      * true, and hasCFields() is false. 
      *
      * \param fields  array of new w (chemical potential) fields
      */
      void setWBasis(DArray< DArray<double> > const & fields);

      /**
      * Set new w fields, in real-space (r-grid) format.
      *
      * This function set values for w fields in r-grid format, but does
      * not set components the symmetry-adapted basis format.  On exit, 
      * w.rgrid() is reset, w().hasData() is true, hasCFields() is false 
      * and w().isSymmetric() is false.
      *
      * \param fields  array of new w (chemical potential) fields
      */
      void setWRGrid(DArray< RField<D> > const & fields);

      /**
      * Set new w fields, in unfolded real-space (r-grid) format.
      *
      * The function parameter "fields" is an unfolded array containing
      * r-grid fields for all monomer types in a single array, with the
      * field for monomer 0 first, followed by the field for monomer 1, 
      * etc. This function sets w().rgrid() but does not set w().basis().
      * On exit, w().hasData is true, while w().isSymmetric() and
      * hasCFields() are both false.
      *
      * \param fields  unfolded array of new chemical potential fields
      */
      void setWRGrid(Field<cudaReal> & fields);

      /**
      * Symmetrize r-grid w-fields, compute basis components.
      *
      * On exit, w().hasData() and w().isSymmetric() are true, while
      * hasCFields() is false.
      */
      void symmetrizeWFields();

      /**
      * Construct trial w-fields from c-fields.
      *
      * This function reads concentration fields in symmetrized basis 
      * format and constructs an initial guess for corresponding chemical 
      * potential fields by setting the Lagrange multiplier field xi to 
      * zero. The result is stored in the system w field container.
      *
      * Upon return, w().hasData() and w().isSymmetric() are true, while 
      * hasCFields() is false. 
      *
      * \param filename  name of input c-field file (basis format)
      */
      void estimateWfromC(const std::string& filename);
   
      ///@}
      /// \name Unit Cell Modifiers 
      ///@{

      /**
      * Set parameters of the associated unit cell.
      *
      * The lattice system declared within the input unitCell must agree 
      * with Domain::lattice() on input if Domain::lattice() is not null. 
      *
      * \param unitCell  new UnitCell<D> (i.e., new parameters)
      */
      void setUnitCell(UnitCell<D> const & unitCell);

      /**
      * Set state of the associated unit cell.
      *
      * The lattice argument must agree with Domain::lattice() on input 
      * if Domain::lattice() is not null, and the size of the parameters
      * array must agree with the expected number of lattice parameters.
      *  
      * \param lattice  lattice system
      * \param parameters  array of new unit cell parameters.
      */
      void setUnitCell(typename UnitCell<D>::LatticeSystem lattice,
                       FSArray<double, 6> const & parameters);

      /**
      * Set parameters of the associated unit cell.
      *
      * The size of the FSArray<double, 6> parameters must agree with the 
      * expected number of parameters for the current lattice type.
      *
      * \param parameters  array of new unit cell parameters.
      */
      void setUnitCell(FSArray<double, 6> const & parameters);

      ///@}
      /// \name Primary SCFT Computations
      ///@{

      /**
      * Solve the modified diffusion equation once, without iteration.
      *
      * This function calls the Mixture::compute() function to solve
      * the statistical mechanics problem for a non-interacting system
      * subjected to the currrent chemical potential fields (wFields
      * and wFieldRGrid). This requires solution of the modified
      * diffusion equation for all polymers, computation of Boltzmann
      * weights for all solvents, computation of molecular partition
      * functions for all species, and computation of concentration
      * fields for blocks and solvents, and computation of overall
      * concentrations for all monomer types. This function does not
      * compute the canonical (Helmholtz) free energy or grand-canonical
      * free energy (i.e., pressure). Upon return, the flag hasCFields
      * is set true.
      *
      * \pre The boolean w().hasData() must be true on entry
      * \post hasCFields is true on exit
      *
      * If argument needStress == true, then this function also calls
      * Mixture<D>::computeStress() to compute the stress.
      *
      * \param needStress  true if stress is needed, false otherwise
      */
      void compute(bool needStress = false);

      /**
      * Iteratively solve a SCFT problem.
      *
      * This function calls the iterator to attempt to solve the SCFT
      * problem for the current mixture and system parameters, using
      * the current chemical potential fields (wFields and wFieldRGrid)
      * and current unit cell parameter values as initial guesses.
      * Upon exist, hasCFields is set true whether or not convergence
      * is obtained to within the desired tolerance. The Helmholtz free
      * energy and pressure are computed if and only if convergence is
      * obtained.
      *
      * \pre The w().hasData() flag must be true on entry, to confirm
      * that chemical potential fields have been set. 
      *
      * \pre The w().isSymmetric() flag must be set true if the chosen 
      * iterator uses a basis representation, and thus requires this.
      *
      * \param isContinuation  true iff continuation within a sweep
      * \return returns 0 for successful convergence, 1 for failure
      */
      int iterate(bool isContinuation = false);

      /**
      * Sweep in parameter space, solving an SCF problem at each point.
      *
      * This function uses a Sweep object that was initialized in the
      * parameter file to solve the SCF problem at a sequence of points
      * along a line in parameter space. The nature of this sequence
      * is determined by implementation of a subclass of Sweep and the
      * parameters passed to the sweep object in the parameter file.
      *
      * \pre The w().hasData() flag must be true on entry.
      * \pre A sweep object must have been created in the parameter file.
      */
      void sweep();
      
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
      * Write timer file to an ostream
      *
      * \param out output stream
      */
      void writeTimers(std::ostream& out);
      
      /**
      * Clear timers 
      */
      void clearTimers();

      ///@}
      /// \name Thermodynamic Properties
      ///@{

      /**
      * Compute free energy density and pressure for current fields.
      *
      * This function should be called after a successful call of
      * iterator().solve(). Resulting values are returned by the
      * freeEnergy() and pressure() accessor functions.
      *
      * \pre w().hasData must be true.
      * \pre hasCFields must be true.
      */
      void computeFreeEnergy();

      /**
      * Get precomputed Helmoltz free energy per monomer / kT.
      *
      * The value retrieved by this function is pre-computed by the
      * computeFreeEnergy() function.
      */
      double fHelmholtz() const;

      /**
      * Get precomputed pressure times monomer volume / kT.
      *
      * The value retrieved by this function is pre-computed by the
      * computeFreeEnergy() function.
      */
      double pressure() const;

      ///@}
      /// \name Thermodynamic Data Output 
      ///@{

      /**
      * Write parameter file to an ostream, omitting the sweep block.
      *
      * This function omits the Sweep block of the parameter file, if
      * any, in order to allow the output produced during a sweep to refer
      * only to parameters relevant to a single state point, and to be
      * rerunnable as a parameter file for a single SCFT calculation.
      */
      void writeParamNoSweep(std::ostream& out) const;

      /**
      * Write thermodynamic properties to a file.
      *
      * This function outputs Helmholtz free energy per monomer, pressure
      * (in units of kT per monomer volume), the volume fraction and 
      * chemical potential of each species, and unit cell parameters.
      *
      * If parameter "out" is a file that already exists, this function
      * will append this information to the end of the file, rather than 
      * overwriting that file. Calling writeParamNoSweep and writeThermo 
      * in succession with the same file will thus produce a single file
      * containing both input parameters and resulting thermodynanic
      * properties.
      *
      * \param out output stream
      */
      void writeThermo(std::ostream& out);

      ///@}
      /// \name Field Output 
      ///@{

      /**
      * Write chemical potential fields in symmetry adapted basis format.
      *
      * \param filename name of output file
      */
      void writeWBasis(const std::string & filename);

      /**
      * Write chemical potential fields in real space grid (r-grid) format.
      *
      * \param filename name of output file
      */
      void writeWRGrid(const std::string & filename) const;

      /**
      * Write concentrations in symmetry-adapted basis format.
      *
      * \param filename name of output file
      */
      void writeCBasis(const std::string & filename);

      /**
      * Write concentration fields in real space grid (r-grid) format.
      *
      * \param filename name of output file
      */
      void writeCRGrid(const std::string & filename) const;

      /**
      * Write c fields for all blocks and solvents in r-grid format.
      *
      * Writes concentrations for all blocks of all polymers and all
      * solvent species in r-grid format. Columns associated with blocks
      * appear ordered by polymer id and then by block id, followed by
      * solvent species ordered by solvent id. 
      *
      * \param filename name of output file
      */
      void writeBlockCRGrid(const std::string & filename) const;

      ///@}
      /// \name Propagator Output 
      ///@{

      /**
      * Write specified slice of a propagator at fixed s in r-grid format.
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
      * Write one propagator for one block, in r-grid format.
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
      * and is passed an automatically generated file name.  The file 
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
      /// \name Crystallographic Data Output 
      ///@{

      /**
      * Output information about stars and symmetrized basis functions.
      *
      * This function opens a file with the specified filename and then
      * calls Basis::outputStars.
      *
      * \param filename name of output file
      */
      void writeStars(const std::string & filename) const;
   
      /**
      * Output information about waves.
      *
      * This function opens a file with the specified filename and then
      * calls Basis::outputWaves.
      *
      * \param filename name of output file 
      */
      void writeWaves(const std::string & filename) const;
   
      /**
      * Output all elements of the space group.
      *
      * \param filename name of output file 
      */
      void writeGroup(std::string const & filename) const;

      ///@}
      /// \name Field File Operations
      ///@{

      /**
      * Convert a field from symmetry-adapted basis to r-grid format.
      *
      * This and other field conversion functions do not change the w 
      * or c fields stored by this System - all required calculations 
      * are performed using temporary or mutable memory. 
      *
      * \param inFileName name of input file
      * \param outFileName name of output file
      */
      void basisToRGrid(const std::string & inFileName,
                        const std::string & outFileName);

      /**
      * Convert a field from real-space grid to symmetrized basis format.
      *
      * This function checks if the input fields have the declared space
      * group symmetry, and prints a warning if it detects deviations
      * that exceed some small threshhold, but proceeds to attempt the
      * conversion even if such an error is detected. Converting a field 
      * that does not have the declared space group symmetry to basis
      * format is a destructive operation that modifies the fields.
      *
      * \param inFileName name of input file
      * \param outFileName name of output file
      */
      void rGridToBasis(const std::string & inFileName,
                        const std::string & outFileName);

      /**
      * Convert fields from Fourier (k-grid) to real-space (r-grid) format.
      *
      * \param inFileName name of input file
      * \param outFileName name of output file
      */
      void kGridToRGrid(const std::string& inFileName, 
                        const std::string& outFileName);
   
      /**
      * Convert fields from real-space (r-grid) to Fourier (k-grid) format.
      *
      * \param inFileName name of input file
      * \param outFileName name of output file
      */
      void rGridToKGrid(const std::string & inFileName, 
                        const std::string & outFileName);

      /**
      * Convert fields from Fourier (k-grid) to symmetrized basis format.
      *
      * This function checks if the input fields have the declared space
      * group symmetry, and prints a warning if it detects deviations
      * that exceed some small threshhold, but proceeds to attempt the
      * conversion even if such an error is detected. Converting a field 
      * that does not have the declared space group symmetry to basis
      * format is a destructive operation that modifies the fields.
      *
      * \param inFileName name of input file
      * \param outFileName name of output file
      */
      void kGridToBasis(const std::string& inFileName, 
                        const std::string& outFileName);
   
      /**
      * Convert fields from symmetrized basis to Fourier (k-grid) format.
      *
      * \param inFileName name of input file (basis format)
      * \param outFileName name of output file (k-grid format)
      */
      void basisToKGrid(const std::string & inFileName, 
                        const std::string & outFileName);
  
      ///@}
      /// \name Member Object Accessors
      ///@{

      /**
      * Get container of chemical potential fields.
      */
      WFieldContainer<D> const & w() const;

      /**
      * Get container of monomer concentration / volume fration fields.
      */
      CFieldContainer<D> const & c() const;

      /**
      * Get Mixture by reference.
      */
      Mixture<D>& mixture();

      /**
      * Get interaction (i.e., excess free energy model) by reference.
      */
      Interaction & interaction();
	
      /**
      * Get interaction (i.e., excess free energy model) by const ref.
      */
      Interaction const & interaction() const;
	
      /**
      * Get Domain by const reference.
      */
      Domain<D> const & domain() const;
      
      /**
      * Get Domain by nonconst reference.
      */
      Domain<D> & domain();

      /**
      * Get the iterator by reference.
      */
      Iterator<D>& iterator();
      
      /**
      * Get the compressor by reference.
      */
      Compressor<D>& compressor();
      
      /**
      * Get Simulator container for Monte Carlo moves.
      */
      Simulator<D>& simulator();

      /**
      * Get crystal UnitCell by const reference.
      */
      UnitCell<D> const & unitCell() const;

      /**
      * Get spatial discretization Mesh by const reference.
      */
      Mesh<D> const & mesh() const;

      /**
      * Get the Basis by reference.
      */
      Basis<D> const & basis() const;

      /**
      * Get the FieldIo by const reference.
      */
      FieldIo<D> const & fieldIo() const;

      /**
      * Get the FFT object by reference.
      */
      FFT<D> const & fft() const ;

      /**
      * Get homogeneous mixture (for reference calculations).
      */
      Homogeneous::Mixture& homogeneous();

      /**
      * Get FileMaster by reference.
      */
      FileMaster& fileMaster();

      ///@}
      /// \name Queries
      ///@{

      /**
      * Have monomer concentration fields (c fields) been computed?
      *
      * Returns  true if and only if monomer concentration fields have
      * been computed by solving the modified diffusion equation for the
      * current w fields.
      */
      bool hasCFields() const;

      /**
      * Does this system have an associated Sweep object?
      */
      bool hasSweep() const;
      
      /**
      * Does this system have a Compressor object?
      */
      bool hasCompressor() const;
      
      /**
      * Does this system have an initialized Simulator?
      */
      bool hasSimulator() const;

      ///@}

   private:

      /**
      * Mixture object (solves MDE for all species).
      */
      Mixture<D> mixture_;

      /**
      * Domain object (crystallography and mesh).
      */
      Domain<D> domain_;
      
      /**
      * Filemaster (holds paths to associated I/O files).
      */
      FileMaster fileMaster_;

      /**
      * Homogeneous mixture, for reference.
      */
      Homogeneous::Mixture homogeneous_;

      /**
      * Pointer to Interaction (excess free energy model).
      */
      Interaction* interactionPtr_;
	
      /**
      * Pointer to an iterator.
      */
      Iterator<D>* iteratorPtr_;

      /**
      * Pointer to iterator factory object
      */
      IteratorFactory<D>* iteratorFactoryPtr_;

      /**
      * Pointer to an Sweep object
      */
      Sweep<D>* sweepPtr_;

      /**
      * Pointer to SweepFactory object
      */
      SweepFactory<D>* sweepFactoryPtr_;
      
      /**
      * Pointer to a compressor.
      */
      Compressor<D>* compressorPtr_;

      /**
      * Pointer to compressor factory object
      */
      CompressorFactory<D>* compressorFactoryPtr_;

      /**
      * Pointer to a simulator.
      */
      Simulator<D>* simulatorPtr_;
      
      /**
      * Pointer to Simulator factory object
      */
      SimulatorFactory<D>* simulatorFactoryPtr_;

      /**
      * Chemical potential fields.
      */
      WFieldContainer<D> w_;

      /**
      * Monomer concentration / volume fraction fields.
      */
      CFieldContainer<D> c_;

      /**
      * Work array of field coefficients for all monomer types.
      *
      * Indexed by monomer typeId, size = nMonomer.
      */
      mutable DArray< DArray<double> > tmpFieldsBasis_;

      /**
      * Work array of fields on real space grid.
      *
      * Indexed by monomer typeId, size = nMonomer.
      */
      mutable DArray< RField<D> > tmpFieldsRGrid_;

      /**
      * Work array of fields on Fourier grid (k-grid).
      *
      * Indexed by monomer typeId, size = nMonomer.
      */
      mutable DArray<RFieldDft<D> > tmpFieldsKGrid_;

      /**
      * Helmholtz free energy per monomer / kT.
      */
      double fHelmholtz_;

      /**
      * Ideal gas contribution to fHelmholtz_. 
      * 
      * This encompasses the internal energy and entropy of 
      * non-interacting free chains in their corresponding 
      * potential fields defined by w_.
      */
      double fIdeal_;

      /**
      * Multi-chain interaction contribution to fHelmholtz_.
      */
      double fInter_;


      /**
      * Pressure times monomer volume / kT.
      *
      * This is -1 times the grand-canonical free energy per monomer,
      * divided by kT.
      */
      double pressure_;

      /**
      * Has the mixture been initialized?
      */
      bool hasMixture_;
      
      /**
      * Does a Simulator exist?
      */
      bool hasSimulator_;

      /**
      * Has memory been allocated for fields in FFT grid formats?
      */
      bool isAllocatedGrid_;

      /**
      * Has memory been allocated for fields in basis format?
      */
      bool isAllocatedBasis_;

      /**
      * Have C fields been computed for the current w fields?
      */
      bool hasCFields_;
      
      /**
      * Has fHelmholtz been computed for the current w and c fields?
      */ 
      bool hasFreeEnergy_;

      /**
      * Dimemsions of the k-grid (discrete Fourier transform grid).
      */
      IntVec<D> kMeshDimensions_;

      /**
      * Work array for r-grid field.
      */
      RField<D> workArray_;

      // Private member functions
      
      /**
      * Allocate memory for fields in grid formats (private).
      *
      * Can be called when mesh is known.
      */
      void allocateFieldsGrid();

      /**
      * Allocate memory for fields in basis formats (private).
      *
      * Can be called only after the basis is initialized.
      */
      void allocateFieldsBasis();

      /**
      * Read field file header and initialize unit cell and basis.
      *
      * \param filename name of field file
      */
      void readFieldHeader(std::string filename);

      /**
      * Read a filename string and echo to log file.
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
      * Used to read filenames in readCommands.
      *
      * \param in  input stream (i.e., input file)
      * \param value  number to read and echo
      */
      void readEcho(std::istream& in, double& value) const;

      /**
      * Initialize Homogeneous::Mixture object (private).
      */
      void initHomogeneous();

   };

   // Inline member functions

   // Get the associated Mixture<D> object.
   template <int D>
   inline Mixture<D>& System<D>::mixture()
   { return mixture_; }

   // Get the Interaction (excess free energy model).
   template <int D>
   inline Interaction & System<D>::interaction()
   {
      UTIL_ASSERT(interactionPtr_);
      return *interactionPtr_;
   }

   // Get the Interaction (excess free energy model).
   template <int D>
   inline Interaction const & System<D>::interaction() const
   {
      UTIL_ASSERT(interactionPtr_);
      return *interactionPtr_;
   }

   template <int D>
   inline Domain<D> const & System<D>::domain() const
   { return domain_; }
   
   // Get the Domain by nonconst reference.
   template <int D>
   inline Domain<D> & System<D>::domain() 
   { return domain_; }

   template <int D>
   inline UnitCell<D> const & System<D>::unitCell() const
   { return domain_.unitCell(); }

   // get the const Mesh<D> object.
   template <int D>
   inline Mesh<D> const & System<D>::mesh() const
   { return domain_.mesh(); }

   // Get the Basis<D> object.
   template <int D>
   inline Basis<D> const & System<D>::basis() const
   {  return domain_.basis(); }

   // Get the FFT<D> object.
   template <int D>
   inline FFT<D> const & System<D>::fft() const
   {  return domain_.fft(); }

   // Get the const FieldIo<D> object.
   template <int D>
   inline FieldIo<D> const & System<D>::fieldIo() const
   {  return domain_.fieldIo(); }
   
   // Get the Iterator.
   template <int D>
   inline Iterator<D>& System<D>::iterator()
   {
      UTIL_ASSERT(iteratorPtr_);
      return *iteratorPtr_;
   }
   
   // Get the Compressor
   template <int D>
   inline Compressor<D>& System<D>::compressor()
   {
      UTIL_ASSERT(compressorPtr_);
      return *compressorPtr_;
   }

   // Get the Simulator.
   template <int D>
   inline Simulator<D>& System<D>::simulator()
   {  
      UTIL_ASSERT(simulatorPtr_);
      return *simulatorPtr_; 
   }

   // Get the FileMaster.
   template <int D>
   inline FileMaster& System<D>::fileMaster()
   {  return fileMaster_; }

   // Get the Homogeneous::Mixture object.
   template <int D>
   inline
   Homogeneous::Mixture& System<D>::homogeneous()
   {  return homogeneous_; }

   // Get container of w fields (const reference)
   template <int D>
   inline
   WFieldContainer<D> const & System<D>::w() const
   {  return w_; }

   // Get container of c fields (const reference)
   template <int D>
   inline
   CFieldContainer<D> const & System<D>::c() const
   {  return c_; }

   // Get precomputed Helmoltz free energy per monomer / kT.
   template <int D>
   inline double System<D>::fHelmholtz() const
   {
      UTIL_CHECK(hasFreeEnergy_);  
      return fHelmholtz_; 
   }

   // Get precomputed pressure (units of kT / monomer volume).
   template <int D>
   inline double System<D>::pressure() const
   {
      UTIL_CHECK(hasFreeEnergy_);  
      return pressure_; 
   }

   // Have the c fields been computed for the current w fields?
   template <int D>
   inline bool System<D>::hasCFields() const
   {  return hasCFields_; }

   // Does this system have an associated Sweep object?
   template <int D>
   inline bool System<D>::hasSweep() const
   {  return (sweepPtr_ != 0); }
   
   // Does the system have a Compressor object?
   template <int D>
   inline bool System<D>::hasCompressor() const
   {  return (compressorPtr_ != 0); }
   
   // Does the system have an initialized Simulator ?
   template <int D>
   inline bool System<D>::hasSimulator() const
   {  return (simulatorPtr_ != 0); }

   #ifndef PSPG_SYSTEM_TPP
   // Suppress implicit instantiation
   extern template class System<1>;
   extern template class System<2>;
   extern template class System<3>;
   #endif

} // namespace Rpg
} // namespace Pscf
//#include "System.tpp"
#endif

#ifndef PSPG_SYSTEM_H
#define PSPG_SYSTEM_H

/*
* PSCF - Polymer Self-Consistent Field Theory
*
* Copyright 2016 - 2019, The Regents of the University of Minnesota
* Distributed under the terms of the GNU General Public License.
*/

#include <util/param/ParamComposite.h>     // base class

#include <pspg/solvers/Mixture.h>          // member
#include <pspg/field/Domain.h>             // member
#include <pspg/field/FieldIo.h>            // member
#include <pspg/solvers/WaveList.h>         // member
#include <pspg/field/DField.h>             // typedef
#include <pspg/field/RDField.h>            // typedef
#include <pspg/field/RDFieldDft.h>         // typedef

#include <pscf/homogeneous/Mixture.h>      // member

#include <util/misc/FileMaster.h>          // member
#include <util/containers/DArray.h>        // member template
#include <util/containers/FSArray.h>       // member template

namespace Pscf {

   class ChiInteraction;

namespace Pspg
{

   template <int D> class Iterator;
   template <int D> class IteratorFactory;
   template <int D> class Sweep;
   template <int D> class SweepFactory;

   using namespace Util;

   /**
   * Main class in SCFT simulation of one system.
   *
   * \ingroup Pscf_Pspg_Module
   */
   template <int D>
   class System : public ParamComposite
   {

   public:

      // Field on a real-space grid (r-grid)
      typedef RDField<D> Field;

      /**
      * Constructor.
      */
      System();

      /**
      * Destructor.
      */
      ~System();

      /// \name Lifetime (Actions)
      //@{

      /**
      * Process command line options.
      *
      * \param argc number of command line arguments
      * \param argv array of argument strings
      */
      void setOptions(int argc, char **argv);

      /**
      * Set number of blocks and number of threads
      *
      * \param nBlocks number of blocks
      * \param nThreads number of threads per block
      */
      void setGpuResources (int nBlocks, int nThreads);

      /**
      * Read input parameters (with opening and closing lines).
      *
      * \param in input parameter stream
      */
      virtual void readParam(std::istream& in);

      /**
      * Read input parameters from default param file.
      */
      void readParam();

      /**
      * Read body of input parameters block (without opening and closing lines).
      *
      * \param in input parameter stream
      */
      virtual void readParameters(std::istream& in);

      /**
      * Write parameter file to an ostream, omitting the sweep block. 
      */
      void writeParam(std::ostream& out);

      /**
      * Read command script.
      *
      * \param in command script file.
      */
      void readCommands(std::istream& in);

      /**
      * Read commands from default command file.
      */
      void readCommands();

      //@}
      /// \name Thermodynamic Properties
      //@{

      /**
      * Compute free energy density and pressure for current fields.
      *
      * This function should be called after a successful call of
      * iterator().solve(). Resulting values are returned by the
      * freeEnergy() and pressure() accessor functions.
      */
      void computeFreeEnergy();

      /**
      * Output thermodynamic properties to a file.
      *
      * This function outputs Helmholtz free energy per monomer,
      * pressure (in units of kT per monomer volume), and the
      * volume fraction and chemical potential of each species.
      *
      * \param out output stream
      */
      void outputThermo(std::ostream& out);

      /**
      * Get precomputed Helmoltz free energy per monomer / kT.
      *
      * The value retrieved by this function is computed by the
      * computeFreeEnergy() function.
      */
      double fHelmholtz() const;

      /**
      * Get precomputed pressure x monomer volume kT.
      *
      * The value retrieved by this function is computed by the
      * computeFreeEnergy() function.
      */
      double pressure() const;

      //@}
      /// \name State Modification Functions (Modify w Fields and Unit Cell)
      ///@{

      /**
      * Read chemical potential fields in symmetry adapted basis format.
      *
      * This function opens and reads the file with name "filename",
      * which must contain chemical potential fields in symmetry-adapted
      * basis format, stores these fields in the system wFields array,
      * converts these fields to real-space grid format and stores the
      * result in the wFieldsRGrid array. On exit hasWFields is set true
      * and hasCFields is false.
      *
      * \param filename name of input w-field basis file
      */
      void readWBasis(const std::string & filename);

      /**
      * Set new w fields, in symmetrized Fourier format.
      *
      * On exit hasWFields and hasSymmetric fields are set true, while
      * hasCFields is set false.
      *
      * \param fields  array of new w (chemical potential) fields
      */  
      void setWBasis(DArray< DArray<double> > const & fields);
 
      /**
      * Set new w fields, in real-space (r-grid) format.
      *
      * On exit, hasWFields is set true while hasSymmetricFields and 
      * hasCFields are both false.
      *
      * \param fields  array of new w (chemical potential) fields
      */  
      void setWRGrid(DArray<Field> const & fields);
 
      /**
      * Set new w fields, in unfoled real-space (r-grid) format.
      *
      * The array fields is an unfolded array that contains fields for
      * all monomer types, with the field for monomer 0 first, etc. 
      *
      * On exit hasWFields is set true while hasSymmetricFields and 
      * hasCFields are both false.
      *
      * \param fields  unfolded array of new w (chemical potential) fields
      */  
      void setWRGrid(DField<cudaReal> & fields);
 
      /**
      * Symmetrize r-grid w-fields, compute basis components.
      *
      * Use this function after setting or reading symmetric w fields 
      * in r-grid format to set corresponding symmetrized fields. The 
      * function assumes that the user knows that the wFieldsRgrid 
      * fields are symmetric, and does not check this.
      *
      * On entry hasWFields must be true.
      * On exit, hasWFields and hasSymmetricFields are true, and 
      * hasCFields is false.
      */  
      void symmetrizeWFields();
 
      /**
      * Set parameters of the associated unit cell.
      *
      * \param unitCell  new UnitCell<D> (i.e., new parameters)
      */
      void setUnitCell(UnitCell<D> const & unitCell);

      /**
      * Set parameters of the associated unit cell.
      *
      * \param parameters  array of new unit cell parameters.
      */
      void setUnitCell(FSArray<double, 6> const & parameters);

      ///@}
      /// \name Chemical Potential Field (W Field) Accessors
      //@{

      /**
      * Get array of chemical potential fields, in a basis.
      *
      * This function returns an array in which each element is an
      * array containing the coefficients of the chemical potential
      * field (w field) in a symmetry-adapted basis for one monomer
      * type. The array capacity is the number of monomer types.
      */
      DArray< DArray<double> >& wFieldsBasis();

      /**
      * Get chemical potential field for one monomer type, in a basis.
      *
      * This function returns an array containing coefficients of the
      * chemical potential field (w field) in a symmetry-adapted basis
      * for a specified monomer type.
      *
      * \param monomerId integer monomer type index
      */
      DArray<double>& wFieldBasis(int monomerId);

      /**
      * Get array of all w fields, in r-space format.
      *
      * This function returns an array in which each element is a
      * WField object containing values of the chemical potential field
      * (w field) on a regular grid for one monomer type. The array
      * capacity is the number of monomer types.
      */
      DArray<Field>& wFieldsRGrid();

      /**
      * Get array of all w fields, in r-space format, by const ref.
      *
      * This function returns an array in which each element is a
      * WField object containing values of the chemical potential field
      * (w field) on a regular grid for one monomer type. The array
      * capacity is the number of monomer types.
      */
      DArray<Field> const & wFieldsRGrid() const;

      /**
      * Get the w field for one monomer type, in r-space format.
      *
      * \param monomerId integer monomer type index
      */
      Field& wFieldRGrid(int monomerId);

      /**
      * Get w field for one monomer type in r-space format by const ref.
      *
      * \param monomerId integer monomer type index
      */
      Field const & wFieldRGrid(int monomerId) const;

      /**
      * Get array of chemical potential fields, in Fourier space.
      *
      * The array capacity is equal to the number of monomer types.
      */
      DArray<RDFieldDft<D> >& wFieldsKGrid();

      /**
      * Get the chemical potential field for one monomer, in Fourier space.
      *
      * \param monomerId integer monomer type index
      */
      RDFieldDft<D>& wFieldKGrid(int monomerId);

      //@{
      /// \name Monomer Concentration / Volume Fraction Fields (C Fields)
      //@{

      /**
      * Get an array of all monomer concentration fields, in a basis
      *
      * This function returns an array in which each element is an
      * array containing the coefficients of the monomer concentration
      * field (cfield) for one monomer type in a symmetry-adapted basis.
      * The array capacity is equal to the number of monomer types.
      */
      DArray< DArray<double> >& cFieldsBasis();

      /**
      * Get the concentration field for one monomer type, in a basis.
      *
      * This function returns an array containing the coefficients of
      * the monomer concentration / volume fraction field (c field)
      * for a specific monomer type.
      *
      * \param monomerId integer monomer type index
      */
      DArray<double>& cFieldBasis(int monomerId);

      /**
      * Get array of all concentration fields (c fields), on a grid.
      *
      * This function returns an array in which each element is the
      * monomer concentration field for one monomer type on a regular
      * grid (an r-grid).
      */
      DArray<Field>& cFieldsRGrid();

      /**
      * Get the concentration (c field) for one monomer type, on a grid.
      *
      * \param monomerId integer monomer type index
      */
      Field& cFieldRGrid(int monomerId);

      /**
      * Get all monomer concentration fields, in Fourier space (k-grid).
      *
      * This function returns an arrray in which each element is the
      * discrete Fourier transform (DFT) of the concentration field
      * (c field) for on monomer type.
      */
      DArray<RDFieldDft<D> >& cFieldsKGrid();

      /**
      * Get the c field for one monomer type, in Fourier space (k-grid).
      *
      * This function returns the discrete Fourier transform (DFT) of the
      * concentration field (c field) for monomer type index monomerId.
      *
      * \param monomerId integer monomer type index
      */
      RDFieldDft<D>& cFieldKGrid(int monomerId);

      //@}
      /// \name Accessors (access objects by reference)
      //@{

      /**
      * Get Mixture by reference.
      */
      Mixture<D>& mixture();

      /**
      * Get spatial discretization Mesh by reference.
      */
      Mesh<D> & mesh();

      /**
      * Get spatial discretization Mesh by const reference.
      */
      Mesh<D> const & mesh() const;

      /**
      * Get crystal UnitCell by reference.
      */
      UnitCell<D> & unitCell();

      /**
      * Get crystal UnitCell by const reference.
      */
      UnitCell<D> const & unitCell() const;

      /**
      * Get Domain by const reference.
      */
      Domain<D> const & domain() const;

      /**
      * Get interaction (i.e., excess free energy model) by reference.
      */
      ChiInteraction& interaction();
	
      /**
      * Get the iterator by reference.
      */
      Iterator<D>& iterator();

      /**
      * Get the Basis by reference.
      */
      Basis<D> & basis();

      /**
      * Get container for wavevector data.
      */
      WaveList<D>& wavelist();

      /**
      * Get the FieldIo by const reference.
      */
      FieldIo<D> const & fieldIo() const;

      /**
      * Get the FFT object by reference.
      */
      FFT<D> & fft();

      /**
      * Get homogeneous mixture (for reference calculations).
      */
      Homogeneous::Mixture& homogeneous();

      /**
      * Get FileMaster by reference.
      */
      FileMaster& fileMaster();

      //@}
      /// \name Accessors (return by value)
      //@{

      /**
      * Get the group name string.
      */
      std::string groupName();

      /**
      * Have monomer chemical potential fields (w fields) been set?
      *
      * A true value is returned if and only if values have been set on a
      * real space grid. The READ_W_BASIS command must immediately convert
      * from a basis to a grid to satisfy this requirement.
      */
      bool hasWFields() const;

      /**
      * Have monomer concentration fields (c fields) been computed?
      *
      * A true value is returned if and only if monomer concentration fields
      * have been computed by solving the modified diffusion equation for the
      * current w fields, and values are known on a grid (cFieldsRGrid).
      */
      bool hasCFields() const;

      /** 
      * Are w-fields symmetric under all elements of the space group?
      *
      * This is true iff the fields were originally input in basis format.
      */
      bool hasSymmetricFields() const;

      ///@}
      /// \name Commands (each corresponds to a command file command)
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
      * If argument needStress == true, then this function also calls
      * Mixture<D>::computeStress() to compute the stress.
      *
      * \param needStress true if stress is needed, false otherwise
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
      * is obtained to within the desired tolerance.  The Helmholtz free
      * energy and pressure are computed if and only if convergence is
      * obtained.
      *
      * \pre The hasWFields flag must be true on entry.
      * \return returns 0 for successful convergence, 1 for failure.
      */
      int iterate();

      /**
      * Sweep in parameter space, solving an SCF problem at each point.
      *
      * This function uses a Sweep object that was initialized in the 
      * parameter file to solve the SCF problem at a sequence of points
      * along a line in parameter space. The nature of this sequence 
      * is determined by implementation of a subclass of Sweep and the
      * parameters passed to the sweep object in the parameter file.
      *
      * An Exception is thrown if sweep() is called when no Sweep has 
      * been created.
      */
      void sweep();

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
      * Write concentration fields in real space (r-grid) format, for each
      * block (or solvent) individually rather than for each species.
      *
      * \param filename name of output file
      */
      void writeBlockCRGrid(const std::string & filename) const;

      /**
      * Write last contour slice of the propagator in real space grid format.
      *
      * \param filename  name of output file
      * \param polymerID  integer index of polymer species
      * \param blockID  integer index of block within polymer
      */
      void writePropagatorRGrid(const std::string & filename,
                                int polymerID, int blockID);

      /**
      * Write all data associated with the converged solution. This
      * includes the full param file, as well as the thermodynamic
      * data (free energy, pressure, phi and mu for each species).
      * 
      * \param filename name of output file
      */
      void writeData(const std::string & filename);

      /**
      * Convert a field from symmetry-adapted basis to r-grid format.
      *
      * This function uses the arrays that stored monomer concentration
      * fields for temporary storage, and thus corrupts any previously
      * stored values. As a result, flag hasCFields is false on output.
      *
      * \param inFileName name of input file
      * \param outFileName name of output file
      */
      void basisToRGrid(const std::string & inFileName,
                        const std::string & outFileName);

      /**
      * Convert a field from real-space grid to symmetrized basis format.
      *
      * This function uses the arrays that stored monomer concentration
      * fields for temporary storage, and thus corrupts any previously
      * stored values. As a result, flag hasCFields is false on return.
      *
      * \param inFileName name of input file
      * \param outFileName name of output file
      */
      void rGridToBasis(const std::string & inFileName,
                        const std::string & outFileName);

      ///@}

      #if 0
      // Additional functions for field-theoretic Monte-Carlo

      RDField<D>& pressureField();

      RDField<D>& compositionField();
      #endif

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
      ChiInteraction* interactionPtr_;
	
      /**
      * Pointer to an iterator.
      */
      Iterator<D>* iteratorPtr_;

      /**
      * Pointer to iterator factory object
      */
      IteratorFactory<D>* iteratorFactoryPtr_;

      /**
      * Container for wavevector data.
      */
      WaveList<D>* wavelistPtr_;

      /**
      * Pointer to an Sweep object
      */
      Sweep<D>* sweepPtr_;

      /**
      * Pointer to SweepFactory object
      */
      SweepFactory<D>* sweepFactoryPtr_;

      /**
      * Array of chemical potential fields for monomer types.
      *
      * Indexed by monomer typeId, size = nMonomer.
      */
      DArray< DArray<double> > wFieldsBasis_;

      /**
      * Array of chemical potential fields for monomer types.
      *
      * Indexed by monomer typeId, size = nMonomer.
      */
      DArray<Field> wFieldsRGrid_;

      /**
      * Work space for chemical potential fields
      */
      DArray<RDFieldDft<D> > wFieldsKGrid_;

      /**
      * Array of concentration fields for monomer types, basis format.
      *
      * Indexed by monomer typeId, size = nMonomer.
      */
      DArray< DArray<double> > cFieldsBasis_;

      /**
      * Array of concentration fields for monomer types, r-grid format.
      */
      DArray<Field> cFieldsRGrid_;

      /**
      * Array of concentration fields for monomer types, k-grid format.
      */
      DArray<RDFieldDft<D> > cFieldsKGrid_;

      /**
      * Work array (size = # of grid points).
      */
      DArray<double> f_;

      /**
      * Work array (size = # of monomer types).
      */
      DArray<double> c_;

      /**
      * Helmholtz free energy per monomer / kT.
      */
      double fHelmholtz_;

      /**
      * Pressure times monomer volume / kT.
      */
      double pressure_;

      /**
      * Has the mixture been initialized?
      */
      bool hasMixture_;

      /**
      * Has memory been allocated for fields?
      */
      bool isAllocated_;

      /**
      * Have W fields been set?
      *
      * True iff wFieldsRGrid_ has been set.
      */
      bool hasWFields_;

      /**
      * Have C fields been computed by solving the MDE ?
      *
      * True iff cFieldsRGrid_ has been set.
      */
      bool hasCFields_;

      /**
      * Does this system have symmetric fields ?
      * 
      * Set true iff WFields are set and were input using the symmetry
      * adapated basis format, and are thus invariant under all elements 
      * of the specified space group. 
      */
      bool hasSymmetricFields_;

      /**
      * Does this system have a Sweep object?
      */
      bool hasSweep_;

      /**
      * Work array of field coefficients for all monomer types.
      *
      * Indexed by monomer typeId, size = nMonomer.
      */
      DArray< DArray<double> > tmpFields_;

      /**
      * Work array of fields on real space grid.
      *
      * Indexed by monomer typeId, size = nMonomer.
      */
       DArray<Field> tmpFieldsRGrid_;

      /**
      * Work array of fields on Fourier grid (k-grid).
      *
      * Indexed by monomer typeId, size = nMonomer.
      */
      DArray<RDFieldDft<D> > tmpFieldsKGrid_;

      /**
      * Dimemsions of the k-grid (discrete Fourier transform grid).
      */
      IntVec<D> kMeshDimensions_;

      /**
      * Work array for r-grid field.
      */
      RDField<D> workArray;

      cudaReal* d_kernelWorkSpace_;

      cudaReal* kernelWorkSpace_;

      /**
      * Allocate memory for fields (private)
      */
      void allocate();

      /**
      * Initialize Homogeneous::Mixture object (private).
      */
      void initHomogeneous();

      /**
      * Read a filename string and echo to log file.
      *
      * Used to read filenames in readCommands.
      *
      * \param in  input stream (i.e., input file)
      * \param string  string to read and echo
      */
      void readEcho(std::istream& in, std::string& string) const;

      #if 0
      // Additional member variables for field-theoretic Monte Carlo

      RDField<D> compositionField_; //rField version

      RDFieldDft<D> compositionKField_; //kField

      RDField<D> pressureField_;

      // Free energy of the new configuration due to random change
      double fHelmholtzOld_;
      #endif

   };

   // Inline member functions

   template <int D>
   inline Domain<D> const & System<D>::domain() const
   { return domain_; }

   template <int D>
   inline UnitCell<D> & System<D>::unitCell()
   { return domain_.unitCell(); }

   template <int D>
   inline UnitCell<D> const & System<D>::unitCell() const
   { return domain_.unitCell(); }

   // get the Mesh<D> object.
   template <int D>
   inline Mesh<D> & System<D>::mesh()
   { return domain_.mesh(); }

   // get the const Mesh<D> object.
   template <int D>
   inline Mesh<D> const & System<D>::mesh() const
   { return domain_.mesh(); }

   // Get the Basis<D> object.
   template <int D>
   inline Basis<D> & System<D>::basis()
   {  return domain_.basis(); }

   // Get the FFT<D> object.
   template <int D>
   inline FFT<D> & System<D>::fft()
   {  return domain_.fft(); }

   // Get the const FieldIo<D> object.
   template <int D>
   inline FieldIo<D> const & System<D>::fieldIo() const
   {  return domain_.fieldIo(); }

   // Get the groupName string.
   template <int D>
   inline std::string System<D>::groupName()
   { return domain_.groupName(); }

   // Get the associated Mixture<D> object.
   template <int D>
   inline Mixture<D>& System<D>::mixture()
   { return mixture_; }

   // Get the Interaction (excess free energy model).
   template <int D>
   inline ChiInteraction& System<D>::interaction()
   {
      UTIL_ASSERT(interactionPtr_);
      return *interactionPtr_;
   }

   // Get the Iterator.
   template <int D>
   inline Iterator<D>& System<D>::iterator()
   {
      UTIL_ASSERT(iteratorPtr_);
      return *iteratorPtr_;
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

   // Get the WaveList<D> object.
   template <int D>
   inline WaveList<D>& System<D>::wavelist()
   {  return *wavelistPtr_; }

   // Get all w field, in a basis.
   template <int D>
   inline
   DArray< DArray<double> >& System<D>::wFieldsBasis()
   {  return wFieldsBasis_; }

   // Get a single monomer w field in basis format.
   template <int D>
   inline
   DArray<double>& System<D>::wFieldBasis(int id)
   {  return wFieldsBasis_[id]; }

   // Get all w fields in r-grid format
   template <int D>
   inline
   DArray< typename System<D>::Field >& System<D>::wFieldsRGrid()
   {  return wFieldsRGrid_; }

   // Get all w fields in r-grid format, by const reference.
   template <int D>
   inline
   DArray<typename System<D>::Field> const & System<D>::wFieldsRGrid() 
   const
   {  return wFieldsRGrid_; }

   // Get a single w field in r-grid format.
   template <int D>
   inline
   typename System<D>::Field& System<D>::wFieldRGrid(int id)
   {  return wFieldsRGrid_[id]; }

   // Get a single w field in r-grid format, by const reference.
   template <int D>
   inline
   typename System<D>::Field const & System<D>::wFieldRGrid(int id) const
   {  return wFieldsRGrid_[id]; }

   template <int D>
   inline
   DArray<RDFieldDft<D> >& System<D>::wFieldsKGrid()
   {  return wFieldsKGrid_; }

   template <int D>
   inline
   RDFieldDft<D>& System<D>::wFieldKGrid(int id)
   {  return wFieldsKGrid_[id]; }

   // Get all monomer concentration fields, in a basis.
   template <int D>
   inline
   DArray< DArray<double> >& System<D>::cFieldsBasis()
   {  return cFieldsBasis_; }

   // Get a single monomer concentration field, in a basis.
   template <int D>
   inline
   DArray<double>& System<D>::cFieldBasis(int monomerId)
   {  return cFieldsBasis_[monomerId]; }

   // Get all monomer concentration fields, on a grid.
   template <int D>
   inline
   DArray< typename System<D>::Field >& System<D>::cFieldsRGrid()
   {  return cFieldsRGrid_; }

   // Get a single monomer concentration field, on a grid.
   template <int D>
   inline typename System<D>::Field& System<D>::cFieldRGrid(int id)
   {  return cFieldsRGrid_[id]; }

   template <int D>
   inline
   DArray<RDFieldDft<D> >& System<D>::cFieldsKGrid()
   {  return cFieldsKGrid_; }

   template <int D>
   inline
   RDFieldDft<D>& System<D>::cFieldKGrid(int id)
   {  return cFieldsKGrid_[id]; }

   // Get precomputed Helmoltz free energy per monomer / kT.
   template <int D>
   inline double System<D>::fHelmholtz() const
   {  return fHelmholtz_; }

   // Get precomputed pressure (units of kT / monomer volume).
   template <int D>
   inline double System<D>::pressure() const
   {  return pressure_; }

   // Have the w fields been set?
   template <int D>
   inline bool System<D>::hasWFields() const
   {  return hasWFields_; }

   // Have the c fields been computed for the current w fields?
   template <int D>
   inline bool System<D>::hasCFields() const
   {  return hasCFields_; }

   #if 0
   // Additional functions for field-theoretic Monte-Carlo

   template <int D>
   inline RDField<D>& System<D>::pressureField()
   { return pressureField_;}

   template <int D>
   inline RDField<D>& System<D>::compositionField()
   { return compositionField_;}
   #endif

   #ifndef PSPG_SYSTEM_TPP
   // Suppress implicit instantiation
   extern template class System<1>;
   extern template class System<2>;
   extern template class System<3>;
   #endif

} // namespace Pspg
} // namespace Pscf
//#include "System.tpp"
#endif

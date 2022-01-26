#ifndef PSPC_SYSTEM_H
#define PSPC_SYSTEM_H

/*
* PSCF - Polymer Self-Consistent Field Theory
*
* Copyright 2016 - 2019, The Regents of the University of Minnesota
* Distributed under the terms of the GNU General Public License.
*/

#include <util/param/ParamComposite.h>     // base class

#include <pspc/solvers/Mixture.h>          // member
#include <pspc/field/Domain.h>             // member
#include <pspc/field/RField.h>             // typedef

#include <pscf/homogeneous/Mixture.h>      // member

#include <util/misc/FileMaster.h>          // member
#include <util/containers/DArray.h>        // member template
#include <util/containers/Array.h>         // function parameter

namespace Pscf { class ChiInteraction; }

namespace Pscf {
namespace Pspc
{

   template <int D> class Iterator;
   template <int D> class IteratorFactory;
   template <int D> class Sweep;
   template <int D> class SweepFactory;

   using namespace Util;

   /**
   * Main class for SCFT simulation of one system.
   *
   * A System has (among other components):
   *
   *    - a Mixture (a container for polymer and solvent solvers)
   *    - a Domain (a description of the crystal domain and discretization)
   *    - an Iterator
   *    - Monomer chemical fields in both basis and grid formats
   *    - Monomer concentration fields in both basis and grid formats
   *
   * \ingroup Pscf_Pspc_Module
   */
   template <int D>
   class System : public ParamComposite
   {

   public:

      /// Base class for WField and CField
      typedef RField<D> Field;

      /// Monomer chemical potential field type.
      typedef typename Propagator<D>::WField WField;

      /// Monomer concentration / volume fraction field type.
      typedef typename Propagator<D>::CField CField;

      /// \name Construction and Destruction
      //@{

      /**
      * Constructor.
      */
      System();

      /**
      * Destructor.
      */
      ~System();

      //@}
      /// \name Lifetime (Primary Actions)
      //@{

      /**
      * Process command line options.
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
      */
      void readParam();

      /**
      * Read body of parameter block (without opening and closing lines).
      *
      * \param in input parameter stream
      */
      virtual void readParameters(std::istream& in);

      /**
      * Read command script from a file.
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
      * Iterator::solve(). Resulting values are stored and then
      * accessed by the fHelmholtz() and pressure() functions.
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
      void outputThermo(std::ostream& out) const;

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
      /// \name Setter functions
      //@{

      /**
      * Set new w fields, in symmetrized Fourier format.
      *
      * \param fields  array of new w (chemical potential) fields
      */  
      void setWBasis(DArray< DArray<double> > const & fields);
 
      /**
      * Set new w fields, in real-space (r-grid) format.
      *
      * \param fields  array of new w (chemical potential) fields
      */  
      void setWRGrid(DArray<WField> const & fields);
 
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

      //@}
      /// \name Chemical Potential Field (w-Field) Accessor Functions
      //@{

      /**
      * Get array of all chemical potential fields expanded in a basis.
      *
      * The array capacity is equal to the number of monomer types.
      */
      DArray< DArray<double> > const & wFields() const;

      /**
      * Get chemical potential field for a specific monomer type.
      *
      * \param monomerId integer monomer type index
      */
      DArray<double> const & wField(int monomerId) const;
      
      /**
      * Get array of all chemical potential fields on an r-space grid.
      *
      * The array capacity is equal to the number of monomer types.
      */
      DArray<WField> const & wFieldsRGrid() const;

      /**
      * Get chemical potential field for one monomer type on r-space grid.
      *
      * \param monomerId integer monomer type index
      */
      WField const & wFieldRGrid(int monomerId) const;

      //@}
      /// \name Concentration Field (c-Field) Accessor Functions
      //@{

      /**
      * Get array of all concentration fields expanded in a basis.
      *
      * The array capacity is equal to the number of monomer types.
      */
      DArray< DArray<double> > const & cFields() const;

      /**
      * Get concentration field for one monomer type expanded in a basis.
      *
      * \param monomerId integer monomer type index
      */
      DArray<double> const & cField(int monomerId) const;

      /**
      * Get array of all concentration fields on r-space grid.
      *
      * The array capacity is equal to the number of monomer types.
      */
      DArray<CField> const & cFieldsRGrid() const;

      /**
      * Get concentration field for one monomer type on r-space grid.
      *
      * \param monomerId integer monomer type index
      */
      CField const & cFieldRGrid(int monomerId) const;

      //@}
      /// \name Miscellaneous Accessors 
      //@{

      /**
      * Get UnitCell (i.e., type and parameters) by const reference.
      */
      UnitCell<D> const & unitCell() const;


      /**
      * Get the spatial discretization mesh by const reference.
      */
      Mesh<D> const & mesh() const;

      /** 
      * Get group name.
      */  
      std::string groupName() const;

      /**
      * Get associated Basis object by reference.
      */
      Basis<D> const & basis() const;

      /**
      * Get associated FFT object.
      */
      FFT<D> const & fft() const;

      /**
      * Get associated FieldIo object.
      */
      FieldIo<D> const & fieldIo() const;

      /**
      * Get the Mixture by reference.
      */
      Mixture<D>& mixture();

      /**
      * Get Interaction (i.e., excess free energy model) by reference.
      */
      ChiInteraction& interaction();

      /**
      * Get the Iterator by reference.
      */
      Iterator<D>& iterator();

      /**
      * Get homogeneous mixture (for reference calculations).
      */
      Homogeneous::Mixture& homogeneous();

      /**
      * Get FileMaster by reference.
      */
      FileMaster& fileMaster();

      /** 
      * Have monomer chemical potential fields (w fields) been set?
      */
      bool hasWFields() const;

      /** 
      * Have monomer concentration fields (c fields) been computed?
      *
      * A true value is returned iff monomer concentration fields have
      * been computed by solving the modified diffusion equation for the
      * current w fields.
      */  
      bool hasCFields() const;

      //@}
      /// \name Commands (correspond to command file commands)
      //@{

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
      * Read chemical potential fields in real-space grid (r-grid) format.
      *
      * This function opens and reads the file with name filename, 
      * which must contain chemical potential fields in real space grid
      * (r-grid) format, stores these fields in the system wFieldsGrid
      * array, converts these fields to symmetrized basis format and
      * stores the result in the wFields array. On exit hasWFields is
      * true and hasCFields is false. 
      *
      * \param filename name of input w-field file in r-grid format
      */
      void readWRGrid(const std::string & filename);

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
      * On exit, hasCFields is set true whether or not convergence is
      * obtained to within the desired tolerance.  The Helmholtz free 
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
      * along a line in parameter space. The nature of this sequence of
      * points is determined by implementation of a subclass of Sweep
      * and the parameters passed to the sweep object in the parameter 
      * file.  The Iterator that is initialized in the parameter file 
      * is called at each state point.
      *
      * An Exception is thrown if this is called when no Sweep has been 
      * created (i.e., if hasSweep_ == false).
      */
      void sweep();

      /**
      * Compare two basis function format fields and output their maximum 
      * and root-mean-squared difference. Requires a parameter file 
      * to set up the system object. 
      */
      void compare(const DArray< DArray<double> > field1, 
                   const DArray< DArray<double> > field2);

      /**
      * Compare two real-space fields and output their maximum and
      * root-mean-squared difference. Requires a parameter file
      * to set up the system object. 
      */
      void compare(const DArray< RField<D> > field1, 
                   const DArray< RField<D> > field2);

      /**
      * Write chemical potential fields in symmetry adapted basis format.
      *
      * \param filename name of output file
      */
      void writeWBasis(const std::string & filename) const;
   
      /**
      * Write chemical potential fields in real space grid (r-grid) format.
      *
      * \param filename name of output file
      */
      void writeWRGrid(const std::string & filename) const;
   
      /**
      * Write concentration fields in symmetry-adapted basis format.
      *
      * \param filename name of output file
      */
      void writeCBasis(const std::string & filename) const;
   
      /**
      * Write concentration fields in real space grid (r-grid) format.
      *
      * \param filename name of output file
      */
      void writeCRGrid(const std::string & filename) const;
   
      /**
      * Convert a field from symmetry-adapted basis to r-grid format.
      *
      * \param inFileName name of input file
      * \param outFileName name of output file
      */
      void basisToRGrid(const std::string & inFileName, 
                        const std::string & outFileName) const;
   
      /**
      * Convert a field from real-space grid to symmetrized basis format.
      *
      * \param inFileName name of input file
      * \param outFileName name of output file
      */
      void rGridToBasis(const std::string & inFileName,
                        const std::string & outFileName) const;
   
      /**
      * Convert fields from Fourier (k-grid) to real-space (r-grid) format.
      *
      * \param inFileName name of input file
      * \param outFileName name of output file
      */
      void kGridToRGrid(const std::string& inFileName, 
                        const std::string& outFileName) const;
   
      /**
      * Convert fields from real-space (r-grid) to Fourier (k-grid) format.
      *
      * \param inFileName name of input file
      * \param outFileName name of output file
      */
      void rGridToKGrid(const std::string & inFileName, 
                        const std::string & outFileName) const;
  
      /** 
      * Check if r-grid fields have the declared space group symmetry.
      *
      * \param inFileName name of input file
      * \return true if fields all have symmetry, false otherwise
      */ 
      bool checkRGridFieldSymmetry(const std::string & inFileName) const;

      /**
      * Construct trial w-fields from c-fields.
      *
      * This function reads concentration fields in symmetrized basis format
      * from a file named inFileName, and constructs an initial guess for 
      * corresponding chemical potential fields by setting the Lagrange 
      * multiplier field xi to zero. The resulting guess is stored in the 
      * internal System wFields array and is also output to a file named 
      * outFileName. Upon return, hasWFields is set true and hasCFields 
      * is set false. 
      *
      * \param inFileName name of input file
      * \param outFileName name of output file
      */
      void rhoToOmega(const std::string& inFileName, 
                      const std::string& outFileName);
   
      /**
      * Output information about stars and symmetrized basis functions.
      *
      * This function opens a file with the specified filename and then
      * calls Basis::outputStars.
      *
      * \param outFileName name of output file
      */
      void outputStars(const std::string & outFileName) const;
   
      /**
      * Output information about waves.
      *
      * This function opens a file with the specified filename and then
      * calls Basis::outputWaves.
      *
      * \param outFileName name of output file
      */
      void outputWaves(const std::string & outFileName) const;

      //@}

   private:

      // Private member variables

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
      * Pointer to Interaction (free energy model).
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
      * Pointer to a Sweep object
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
      DArray< DArray<double> > wFields_;

      /**
      * Array of chemical potential fields for monomer types.
      *
      * Indexed by monomer typeId, size = nMonomer.
      */
      DArray<WField> wFieldsRGrid_;

      /**
      * Array of concentration fields for monomer types.
      *
      * Indexed by monomer typeId, size = nMonomer.
      */
      DArray< DArray<double> > cFields_;

      /**
      * Array of concentration fields on real space grid.
      *
      * Indexed by monomer typeId, size = nMonomer.
      */
      DArray<CField> cFieldsRGrid_;

      /**
      * Work array of field coefficients for all monomer types.
      *
      * Indexed by monomer typeId, size = nMonomer.
      */
      mutable DArray< DArray<double> > tmpFields_;

      /**
      * Work array of fields on real space grid.
      *
      * Indexed by monomer typeId, size = nMonomer.
      */
      mutable DArray<CField> tmpFieldsRGrid_;

      /**
      * Work array of fields on Fourier grid (k-grid).
      *
      * Indexed by monomer typeId, size = nMonomer.
      */
      mutable DArray<RFieldDft<D> > tmpFieldsKGrid_;

      /**
      * Work array (size = # of grid points).
      */
      mutable DArray<double> f_;

      /**
      * Work array (size = # of monomer types).
      */
      mutable DArray<double> c_;

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
      * True iff both wFields_ and wFieldsRGrid_ are set and consistent.
      */
      bool hasWFields_;

      /**
      * Have C fields been computed by solving MDEs for current w fields?
      *
      * Set true when c fields are computed, set false when w fields or
      * unit cell are reset.
      */
      bool hasCFields_;

      /**
      * Does this system have a Sweep object?
      */
      bool hasSweep_;
      
      /**
      * Does this system have an iterator object?
      */
      // bool hasIterator_;

      // Private member functions

      /**
      * Allocate memory for fields (private)
      */
      void allocate();

      /**
      * Initialize Homogeneous::Mixture object.
      */
      void initHomogeneous();

      /**
      * Reader header of field file (fortran pscf format)
      *
      * \param in input stream (i.e., input file)
      */
      void readFieldHeader(std::istream& in);

      /**
      * Write header for field file (fortran pscf format)
      *
      * \param out output stream (i.e., output file)
      */
      void writeFieldHeader(std::ostream& out) const;

      /**
      * Read a filename string and echo to log file.
      *
      * Used to read filenames in readCommands.
      *
      * \param in  input stream (i.e., input file)
      * \param string  string to read and echo
      */
      void readEcho(std::istream& in, std::string& string) const;

   };

   // Inline member functions

   // Get the associated UnitCell<D> object by const reference.
   template <int D>
   inline UnitCell<D> const & System<D>::unitCell() const
   { return domain_.unitCell(); }

   // Get the Mesh<D> object.
   template <int D>
   inline Mesh<D> const & System<D>::mesh() const
   { return domain_.mesh(); }

   // Get the associated Mixture object.
   template <int D>
   inline Mixture<D>& System<D>::mixture()
   { return mixture_; }

   // Get the Basis<D> object.
   template <int D>
   inline Basis<D> const & System<D>::basis() const
   {  return domain_.basis(); }

   // Get the FFT<D> object.
   template <int D>
   inline FFT<D> const & System<D>::fft() const
   {  return domain_.fft(); }

   // Get the FieldIo<D> object.
   template <int D>
   inline FieldIo<D> const & System<D>::fieldIo() const
   {  return domain_.fieldIo(); }

   // Get the groupName string.
   template <int D>
   inline std::string System<D>::groupName() const
   { return domain_.groupName(); }

   // Get the FileMaster.
   template <int D>
   inline FileMaster& System<D>::fileMaster()
   {  return fileMaster_; }

   // Get the Homogeneous::Mixture object.
   template <int D>
   inline Homogeneous::Mixture& System<D>::homogeneous()
   {  return homogeneous_; }

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

   template <int D>
   inline
   DArray< DArray<double> > const & System<D>::wFields() const
   {  return wFields_; }

   template <int D>
   inline
   DArray<double> const & System<D>::wField(int id) const
   {  return wFields_[id]; }

   // Get an array of monomer chemical potential fields on r-space grids.
   template <int D>
   inline 
   DArray< typename System<D>::WField > const & System<D>::wFieldsRGrid() 
   const
   {  return wFieldsRGrid_; }

   // Get a single monomer chemical potential field on an r-space grid.
   template <int D>
   inline 
   typename System<D>::WField const & System<D>::wFieldRGrid(int id) const
   {  return wFieldsRGrid_[id]; }

   template <int D>
   inline
   DArray< DArray<double> > const & System<D>::cFields() const
   { return cFields_; }

   template <int D>
   inline
   DArray<double> const & System<D>::cField(int id) const
   { return cFields_[id]; }

   // Get array of all monomer concentration fields on grids.
   template <int D>
   inline
   DArray< typename System<D>::CField > const & System<D>::cFieldsRGrid() 
   const
   {  return cFieldsRGrid_; }

   // Get a single monomer concentration field on an r-space grid.
   template <int D>
   inline typename System<D>::CField const & System<D>::cFieldRGrid(int id)
   const
   {  return cFieldsRGrid_[id]; }

   // Have the w fields been set?
   template <int D>
   inline bool System<D>::hasWFields() const
   {  return hasWFields_; }

   // Have the c fields been computed for the current w fields?
   template <int D>
   inline bool System<D>::hasCFields() const
   {  return hasCFields_; }

   // Get the precomputed Helmoltz free energy per monomer / kT.
   template <int D>
   inline double System<D>::fHelmholtz() const
   {  return fHelmholtz_; }

   // Get the precomputed pressure (units of kT / monomer volume).
   template <int D>
   inline double System<D>::pressure() const
   {  return pressure_; }

   #ifndef PSPC_SYSTEM_TPP
   // Suppress implicit instantiation
   extern template class System<1>;
   extern template class System<2>;
   extern template class System<3>;
   #endif

} // namespace Pspc
} // namespace Pscf
#endif

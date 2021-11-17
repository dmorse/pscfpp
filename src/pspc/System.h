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
   * Main class in SCFT simulation of one system.
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
      /// \name Chemical Potential Field (w-Field) Accessor Functions
      //@{

      /**
      * Get array of all chemical potential fields expanded in a basis.
      *
      * The array capacity is equal to the number of monomer types.
      */
      DArray< DArray<double> >& wFields();

      /**
      * Get chemical potential field for a specific monomer type.
      *
      * \param monomerId integer monomer type index
      */
      DArray<double>& wField(int monomerId);
      
      /**
      * Get array of all chemical potential fields on an r-space grid.
      *
      * The array capacity is equal to the number of monomer types.
      */
      DArray<WField>& wFieldsRGrid();

      /**
      * Get chemical potential field for one monomer type on r-space grid.
      *
      * \param monomerId integer monomer type index
      */
      WField& wFieldRGrid(int monomerId);

      /**
      * Get array of all chemical potential fields in k-space.
      * 
      * The array capacity is equal to the number of monomer types.
      */
      DArray<RFieldDft<D> >& wFieldsKGrid();

      /**
      * Get chemical potential field for one monomer type in k-space.
      *
      * \param monomerId integer monomer type index
      */
      RFieldDft<D>& wFieldKGrid(int monomerId);

      //@}
      /// \name Concentration Field (c-Field) Accessor Functions
      //@{

      /**
      * Get array of all concentration fields expanded in a basis.
      *
      * The array capacity is equal to the number of monomer types.
      */
      DArray< DArray<double> >& cFields();

      /**
      * Get concentration field for one monomer type expanded in a basis.
      *
      * \param monomerId integer monomer type index
      */
      DArray<double>& cField(int monomerId);

      /**
      * Get array of all concentration fields on r-space grid.
      *
      * The array capacity is equal to the number of monomer types.
      */
      DArray<CField>& cFieldsRGrid();

      /**
      * Get concentration field for one monomer type on r-space grid.
      *
      * \param monomerId integer monomer type index
      */
      CField& cFieldRGrid(int monomerId);

      /**
      * Get array of all concentration fields in k-space.
      *
      * The array capacity is equal to the number of monomer types.
      */
      DArray<RFieldDft<D> >& cFieldsKGrid();

      /**
      * Get concentration field for one monomer type on k-space grid.
      *
      * \param monomerId integer monomer type index
      */
      RFieldDft<D>& cFieldKGrid(int monomerId);

      //@}
      /// \name Miscellaneous Accessors 
      //@{

      /**
      * Get Mixture by reference.
      */
      Mixture<D>& mixture();

      /**
      * Get spatial discretization mesh by reference.
      */
      Mesh<D>& mesh();

      /**
      * Get UnitCell (i.e., lattice type and parameters) by reference.
      */
      UnitCell<D>& unitCell();

      /**
      * Get Interaction (i.e., excess free energy model) by reference.
      */
      ChiInteraction& interaction();

      /**
      * Get the Iterator by reference.
      */
      Iterator<D>& iterator();

      /**
      * Get associated Basis object by reference.
      */
      Basis<D>& basis();

      /**
      * Get associated FFT object.
      */
      FFT<D>& fft();

      /**
      * Get associated FieldIo object.
      */
      FieldIo<D>& fieldIo();

      /**
      * Get homogeneous mixture (for reference calculations).
      */
      Homogeneous::Mixture& homogeneous();

      /**
      * Get FileMaster by reference.
      */
      FileMaster& fileMaster();

      /** 
      * Get group name.
      */  
      std::string groupName() const;

      /** 
      * Have monomer chemical potential fields (w fields) been set?
      *
      * A true value is returned if and only if consistent values have 
      * been set for both components in a symmetrized basis (wFields) and 
      * for values on a regular real space grid (wFieldsRGrid). Commands 
      * that read w fields from file in either of these formats must 
      * immediately convert to the other.
      */
      bool hasWFields() const;

      /** 
      * Have monomer concentration fields (c fields) been computed?
      *
      * A true value is returned if and only if monomer concentration fields
      * have been computed by solving the modified diffusion equation for 
      * the current w fields, and consistent values have been set for both 
      * values on an r-grid (wFieldsRGrid) and for coefficients in a basis 
      * (cFields).  To satisfy this requirement, solution of the MDE on a 
      * r-grid should always be immediately followed by conversion of c 
      * fields to a basis.
      */  
      bool hasCFields() const;

      //@}
      /// \name Commands (one-to-one correspondence with command file commands)
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
      * Solve the modified diffusion equation once, without iteration.
      *
      * This function calls the mixture().compute() function to solve
      * the statistical mechanics problem for a non-interacting system
      * subjected to the currrent chemical potential fields (wFields
      * and wFieldRGrid). This requires solution of the modified 
      * diffusion equation for all polymers, computation of Boltzmann
      * weights for all solvents, computation and setting of molecular 
      * partition functions for all species, and computation and 
      * setting of block and setting of concentration fields for all
      * blocks and solvents, and computation of overall concentrations
      * for all monomer types. This function does not compute the 
      * canonical (Helmholtz) free energy or grand-canonical free
      * energy (i.e., pressure). Upon eturn, the flag hasCFields is 
      * set true.
      */
      void compute();
   
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
      void writeWRGrid(const std::string & filename);
   
      /**
      * Write concentrations in symmetry-adapted basis format.
      *
      * \param filename name of output file
      */
      void writeCBasis(const std::string & filename);
   
      /**
      * Write concentration fields in real space grid format.
      *
      * \param filename name of output file
      */
      void writeCRGrid(const std::string & filename);
   
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
   
      /**
      * Convert fields from Fourier (k-grid) to real-space (r-grid) format.
      *
      * This function corrupts monomer concentration field arrays, and so
      * cFields is set false on return.
      *
      * \param inFileName name of input file
      * \param outFileName name of output file
      */
      void kGridToRGrid(const std::string& inFileName, 
                        const std::string& outFileName);
   
      /**
      * Convert fields from real-space (r-grid) to Fourier (k-grid) format.
      *
      * This function corrupts monomer concentration field arrays, and so
      * cFields is set false on return.
      *
      * \param inFileName name of input file
      * \param outFileName name of output file
      */
      void rGridToKGrid(const std::string & inFileName, 
                        const std::string & outFileName);
  
      /** 
      * Check if r-grid fields have the declared space group symmetry.
      *
      * \param inFileName name of input file
      * \return true if fields all have symmetry, false otherwise
      */ 
      bool checkRGridFieldSymmetry(const std::string & inFileName);

      /**
      * Construct proposed chemical potential fields from concentration fields.
      *
      * This function reads concentration fields in symmetrized basis format
      * from a file named inFileName, constructs an initial guess for 
      * corresponding chemical potential fields by setting the Lagrange 
      * multiplier field to zero, stores the result in the system wFields 
      * array and also outputs the result to file named outFileName. Upon
      * return, hasWFields is set true and hasCFields is set false. 
      *
      * \param inFileName name of input file
      * \param outFileName name of output file
      */
      void rhoToOmega(const std::string& inFileName, 
                      const std::string& outFileName);
   
      /**
      * Output information about stars and symmetry-adapted basis functions.
      *
      * \param outFileName name of output file
      */
      void outputStars(const std::string & outFileName);
   
      /**
      * Output information about waves.
      *
      * \param outFileName name of output file
      */
      void outputWaves(const std::string & outFileName);

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

      #if 0
      /*
      * Crystallographic unit cell (crystal system and cell parameters).
      */
      UnitCell<D> unitCell_;

      /**
      * Spatial discretization mesh.
      */
      Mesh<D> mesh_;

      /**
      * FFT object to be used by iterator
      */
      FFT<D> fft_;

      /**
      * Pointer to a Basis object
      */
      Basis<D> basis_;

      /**
      * FieldIo object for field input/output operations
      */
      FieldIo<D> fieldIo_;

      /**
      * Group name.
      */
      std::string groupName_;

      #endif

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

      #if 0
      /**
      * Pointer to an Sweep object
      */
      Sweep<D>* sweepPtr_;

      /**
      * Pointer to SweepFactory object
      */
      SweepFactory<D>* sweepFactoryPtr_;
      #endif

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
      * Work space for chemical potential fields
      */
      DArray<RFieldDft<D> > wFieldsKGrid_;

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
      * Array of concentration fields on Fourier grid (k-grid).
      *
      * Indexed by monomer typeId, size = nMonomer.
      */
      DArray<RFieldDft<D> > cFieldsKGrid_;

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
      DArray<CField> tmpFieldsRGrid_;

      /**
      * Work array of fields on Fourier grid (k-grid).
      *
      * Indexed by monomer typeId, size = nMonomer.
      */
      DArray<RFieldDft<D> > tmpFieldsKGrid_;

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

      #if 0
      /**
      * Has the UnitCell been initialized?
      */
      bool hasUnitCell_;

      /**
      * Has the Mesh been initialized?
      */
      bool hasMesh_;
      #endif

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
      * Have C fields been computed by solving the MDE ?
      *
      * True iff both cFields_ and cFieldsRGrid_ are set and consistent.
      */
      bool hasCFields_;

      #if 0
      /**
      * Does this system have a Sweep object?
      */
      bool hasSweep_;
      #endif

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

   // Get the associated Mixture object.
   template <int D>
   inline Mixture<D>& System<D>::mixture()
   { return mixture_; }

   // Get the associated UnitCell<D> object.
   template <int D>
   inline UnitCell<D>& System<D>::unitCell()
   { return domain_.unitCell(); }

   // Get the Mesh<D> object.
   template <int D>
   inline Mesh<D>& System<D>::mesh()
   { return domain_.mesh(); }

   // Get the Basis<D> object.
   template <int D>
   inline Basis<D>& System<D>::basis()
   {  return domain_.basis(); }

   // Get the FFT<D> object.
   template <int D>
   inline FFT<D>& System<D>::fft()
   {  return domain_.fft(); }

   // Get the FieldIo<D> object.
   template <int D>
   inline FieldIo<D>& System<D>::fieldIo()
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
   DArray< DArray<double> >& System<D>::wFields()
   {  return wFields_; }

   template <int D>
   inline
   DArray<double>& System<D>::wField(int id)
   {  return wFields_[id]; }

   // Get an array of monomer chemical potential fields on r-space grids.
   template <int D>
   inline 
   DArray< typename System<D>::WField >& System<D>::wFieldsRGrid()
   {  return wFieldsRGrid_; }

   // Get a single monomer chemical potential field on an r-space grid.
   template <int D>
   inline 
   typename System<D>::WField& System<D>::wFieldRGrid(int id)
   {  return wFieldsRGrid_[id]; }

   template <int D>
   inline
   DArray<RFieldDft<D> >& System<D>::wFieldsKGrid()
   { return wFieldsKGrid_; }

   template <int D>
   inline
   RFieldDft<D>& System<D>::wFieldKGrid(int id)
   { return wFieldsKGrid_[id]; }

   template <int D>
   inline
   DArray< DArray<double> >& System<D>::cFields()
   { return cFields_; }

   template <int D>
   inline
   DArray<double>& System<D>::cField(int id)
   { return cFields_[id]; }

   template <int D>
   inline
   DArray<RFieldDft<D> >& System<D>::cFieldsKGrid()
   { return cFieldsKGrid_; }

   template <int D>
   inline
   RFieldDft<D>& System<D>::cFieldKGrid(int id)
   { return cFieldsKGrid_[id]; }

   // Get array of all monomer concentration fields on grids.
   template <int D>
   inline
   DArray< typename System<D>::CField >& System<D>::cFieldsRGrid()
   {  return cFieldsRGrid_; }

   // Get a single monomer concentration field on an r-space grid.
   template <int D>
   inline typename System<D>::CField& System<D>::cFieldRGrid(int id)
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

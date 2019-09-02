#ifndef PSSP_SYSTEM_H
#define PSSP_SYSTEM_H

/*
* PSCF - Polymer Self-Consistent Field Theory
*
* Copyright 2016, The Regents of the University of Minnesota
* Distributed under the terms of the GNU General Public License.
*/

#include <util/param/ParamComposite.h>     // base class
#include <pssp/solvers/Mixture.h>          // member
#include <pssp/field/FFT.h>                // member
#include <pssp/basis/Basis.h>              // member
#include <pssp/field/FieldIo.h>            // member
#include <pscf/mesh/Mesh.h>                // member
#include <pssp/iterator/AmIterator.h>
#include <pscf/crystal/UnitCell.h>         // member
#include <pscf/homogeneous/Mixture.h>      // member
#include <util/misc/FileMaster.h>          // member
#include <util/containers/DArray.h>        // member template
#include <util/containers/Array.h>         // function parameter
#include <pssp/field/RField.h>             // typedef

namespace Pscf { class ChiInteraction; }

namespace Pscf {
namespace Pssp
{
   template <int D> class AmIterator;
   class Sweep;
   class SweepFactory;

   using namespace Util;

   /**
   * Main class in SCFT simulation of one system.
   *
   * \ingroup Pscf_Pssp_Module
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
      /// \name Field File I/O
      //@{

      /**
      * Read concentration or chemical potential field components from file.
      *
      * This function reads components in a symmetry adapted basis from 
      * file in.
      *
      * The capacity of DArray fields is equal to nMonomer, and element
      * fields[i] is a DArray containing components of the field 
      * associated with monomer type i.
      *
      * \param in input stream (i.e., input file)
      */
      void readFields(std::istream& in, DArray< DArray <double> >& fields);

      /**
      * Read concentration or chemical potential field components from file.
      *
      * This function opens an input file with the specified filename, 
      * reads components in symmetry-adapted form from that file, and 
      * closes the file.
      *
      * \param filename name of input file
      */
      void readFields(std::string filename, 
                      DArray< DArray <double> >& fields);

      /**
      * Write concentration or chemical potential field components to file.
      *
      * This function writes components in a symmetry adapted basis.
      *
      * \param out output stream (i.e., output file)
      * \param fields array of fields for different species
      */
      void writeFields(std::ostream& out, 
                       DArray< DArray <double> > const & fields);

      /**
      * Write concentration or chemical potential field components to file.
      *
      * This function opens an output file with the specified filename, 
      * writes components in symmetry-adapted form to that file, and then
      * closes the file. 
      *
      * \param filename name of input file
      */
      void writeFields(std::string filename, 
                       DArray< DArray <double> > const & fields);

      /**
      * Read array of RField objects (fields on an r-space grid) from file.
      *
      * The capacity of array fields is equal to nMonomer, and element
      * fields[i] is the RField<D> associated with monomer type i.
      * 
      * \param in input stream (i.e., input file)
      */
      void readRFields(std::istream& in, DArray< RField<D> >& fields);

      /**
      * Write array of RField objects (fields on an r-space grid) to file.
      *
      * \param out output stream (i.e., output file)
      * \param fields array of fields for different species
      */
      void writeRFields(std::ostream& out, 
                        DArray< RField<D> > const& fields);

      /**
      * Read array of RFieldDft objects (k-space fields) from file.
      *
      * The capacity of the array is equal to nMonomer, and element
      * fields[i] is the discrete Fourier transform of the field for 
      * monomer type i.
      * 
      * \param in input stream (i.e., input file)
      * \param fields array of fields for different species
      */
      void readKFields(std::istream& in, DArray< RFieldDft<D> >& fields);

      /**
      * Write array of RFieldDft objects (k-space fields) to file.
      *
      * \param out output stream (i.e., output file)
      * \param fields array of fields for different species
      */
      void writeKFields(std::ostream& out, 
                        DArray< RFieldDft<D> > const& fields);

      //@}
      /// \name Field Accessor Functions
      //@{

      /**
      * Get array of all chemical potential fields expanded in a basis.
      *
      * The array capacity is equal to the number of monomer types.
      */
      DArray<DArray <double> >& wFields();

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
      DArray<WField>& wFieldGrids();

      /**
      * Get chemical potential field for one monomer type on r-space grid.
      *
      * \param monomerId integer monomer type index
      */
      WField& wFieldGrid(int monomerId);

      /**
      * Get array of all chemical potential fields in k-space.
      * 
      * The array capacity is equal to the number of monomer types.
      */
      DArray<RFieldDft<D> >& wFieldDfts();

      /**
      * Get chemical potential field for one monomer type in k-space.
      *
      * \param monomerId integer monomer type index
      */
      RFieldDft<D>& wFieldDft(int monomerId);

      /**
      * Get array of all concentration fields expanded in a basis.
      *
      * The array capacity is equal to the number of monomer types.
      */
      DArray<DArray <double> >& cFields();

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
      DArray<CField>& cFieldGrids();

      /**
      * Get concentration field for one monomer type on r-space grid.
      *
      * \param monomerId integer monomer type index
      */
      CField& cFieldGrid(int monomerId);

      /**
      * Get array of all concentration fields in k-space.
      *
      * The array capacity is equal to the number of monomer types.
      */
      DArray<RFieldDft<D> >& cFieldDfts();

      /**
      * Get concentration field for one monomer type on k-space grid.
      *
      * \param monomerId integer monomer type index
      */
      RFieldDft<D>& cFieldDft(int monomerId);

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
      //temporarily changed to allow testing on member functions
      AmIterator<D>& iterator();

      /**
      * Get associated Basis object by reference.
      */
      Basis<D>& basis();

      /**
      * Get associated FFT object.
      */
      FFT<D>& fft();

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
      std::string groupName();

      //@}

   private:

      /**
      * Mixture object (solves MDE for all species).
      */
      Mixture<D> mixture_;

      /**
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
      * Group name.
      */
      std::string groupName_;

      /**
      * Pointer to a Basis object
      */
      Basis<D> basis_;

      /**
      * Filemaster (holds paths to associated I/O files).
      */
      FileMaster fileMaster_;

      /**
      * FieldIo object for field input/output operations
      */
      FieldIo<D> fieldIo_;

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
      AmIterator<D>* iteratorPtr_;

      #if 0
      /**
      * Pointer to an Sweep object
      */
      Sweep* sweepPtr_;

      /**
      * Pointer to SweepFactory object
      */
      SweepFactory* sweepFactoryPtr_;
      #endif

      /**
      * Array of chemical potential fields for monomer types.
      *
      * Indexed by monomer typeId, size = nMonomer.
      */
      DArray<DArray <double> > wFields_;

      /**
      * Array of chemical potential fields for monomer types.
      *
      * Indexed by monomer typeId, size = nMonomer.
      */
      DArray<WField> wFieldGrids_;

      /**
      * Work space for chemical potential fields
      */
      DArray<RFieldDft<D> > wFieldDfts_;

      /**
      * Array of concentration fields for monomer types.
      *
      * Indexed by monomer typeId, size = nMonomer.
      */
      DArray<DArray <double> > cFields_;

      /**
      * Array of concentration fields on real space grid.
      *
      * Indexed by monomer typeId, size = nMonomer.
      */
      DArray<CField> cFieldGrids_;

      /**
      * Array of concentration fields on Fourier space grid.
      *
      * Indexed by monomer typeId, size = nMonomer.
      */
      DArray<RFieldDft<D> > cFieldDfts_;

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
      * Has the UnitCell been initialized?
      */
      bool hasUnitCell_;

      /**
      * Has the Mesh been initialized?
      */
      bool hasMesh_;

      /**
      * Have initial chemical potential fields been read from file?
      */
      bool hasFields_;

      /**
      * Does this system have a Sweep object?
      */
      bool hasSweep_;

      /**
      * Allocate memory for fields (private)
      */
      void allocateFields();

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

   };

   // Inline member functions

   /*
   * Get the associated Mixture object.
   */
   template <int D>
   inline Mixture<D>& System<D>::mixture()
   { return mixture_; }

   /*
   * Get the UnitCell<D>.
   */
   template <int D>
   inline UnitCell<D>& System<D>::unitCell()
   { return unitCell_; }

   /*
   * Get the mesh.
   */
   template <int D>
   inline Mesh<D>& System<D>::mesh()
   { return mesh_; }

   template <int D>
   inline FFT<D>& System<D>::fft()
   {  return fft_; }

   /*
   * Get group name.
   */
   template <int D>
   inline std::string System<D>::groupName()
   { return groupName_; }

   template <int D>
   inline Basis<D>& System<D>::basis()
   {  return basis_; }

   /*
   * Get the FileMaster.
   */
   template <int D>
   inline FileMaster& System<D>::fileMaster()
   {  return fileMaster_; }

   /*
   * Get the Homogeneous::Mixture object.
   */
   template <int D>
   inline 
   Homogeneous::Mixture& System<D>::homogeneous()
   {  return homogeneous_; }

   /*
   * Get the Interaction (excess free energy model).
   */
   template <int D>
   inline ChiInteraction& System<D>::interaction()
   {
      UTIL_ASSERT(interactionPtr_);
      return *interactionPtr_;
   }

   /*
   * Get the Iterator.
   */
   template <int D>
   inline AmIterator<D>& System<D>::iterator()
   {
      UTIL_ASSERT(iteratorPtr_);
      return *iteratorPtr_;
   }

   template <int D>
   inline
   DArray<DArray <double> >& System<D>::wFields()
   {  return wFields_; }

   template <int D>
   inline
   DArray<double>& System<D>::wField(int id)
   {  return wFields_[id]; }

   /*
   * Get an array of all monomer chemical potential fields.
   */
   template <int D>
   inline 
   DArray< typename System<D>::WField >& System<D>::wFieldGrids()
   {  return wFieldGrids_; }

   /*
   * Get a single monomer chemical potential field.
   */
   template <int D>
   inline 
   typename System<D>::WField& System<D>::wFieldGrid(int id)
   {  return wFieldGrids_[id]; }

   template <int D>
   inline
   DArray<RFieldDft<D> >& System<D>::wFieldDfts()
   { return wFieldDfts_; }

   template <int D>
   inline
   RFieldDft<D>& System<D>::wFieldDft(int id)
   { return wFieldDfts_[id]; }

   template <int D>
   inline
   DArray<DArray <double> >& System<D>::cFields()
   { return cFields_; }

   template <int D>
   inline
   DArray<double>& System<D>::cField(int id)
   { return cFields_[id]; }

   template <int D>
   inline
   DArray<RFieldDft<D> >& System<D>::cFieldDfts()
   { return cFieldDfts_; }

   template <int D>
   inline
   RFieldDft<D>& System<D>::cFieldDft(int id)
   { return cFieldDfts_[id]; }

   /*
   * Get array of all monomer concentration fields.
   */
   template <int D>
   inline
   DArray< typename System<D>::CField >& System<D>::cFieldGrids()
   {  return cFieldGrids_; }

   /*
   * Get a single monomer concentration field.
   */
   template <int D>
   inline typename System<D>::CField& System<D>::cFieldGrid(int id)
   {  return cFieldGrids_[id]; }

   /*
   * Get precomputed Helmoltz free energy per monomer / kT.
   */
   template <int D>
   inline double System<D>::fHelmholtz() const
   {  return fHelmholtz_; }

   /*
   * Get precomputed pressure (units of kT / monomer volume).
   */
   template <int D>
   inline double System<D>::pressure() const
   {  return pressure_; }

} // namespace Pssp
} // namespace Pscf
#include "System.tpp"
#endif

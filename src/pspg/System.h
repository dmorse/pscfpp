#ifndef PSPG_SYSTEM_H
#define PSPG_SYSTEM_H

/*
* PSCF - Polymer Self-Consistent Field Theory
*
* Copyright 2016 - 2019, The Regents of the University of Minnesota
* Distributed under the terms of the GNU General Public License.
*/

#include <pspg/iterator/AmIterator.h>
#include <pspg/field/FieldIo.h>            // member
#include <pspg/solvers/Mixture.h>          // member
#include <pspg/wavelist/WaveList.h>        // member
#include <pspg/field/RDField.h>            // typedef
#include <pspg/field/RDFieldDft.h>         // typedef

#include <pscf/crystal/Basis.h>            // member
#include <pscf/mesh/Mesh.h>                // member
#include <pscf/crystal/UnitCell.h>         // member
#include <pscf/homogeneous/Mixture.h>      // member
#include <pscf/inter/ChiInteraction.h>     // member

#include <util/param/ParamComposite.h>     // base class
#include <util/misc/FileMaster.h>          // member
#include <util/containers/DArray.h>        // member template
#include <util/containers/Array.h>         // function parameter

namespace Pscf {
namespace Pspg
{
   class Sweep;
   class SweepFactory;

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

      /// Base class for WField and CField
      typedef RDField<D> Field;

      /// Monomer chemical potential field type.
      typedef typename Propagator<D>::WField WField;

      /// Monomer concentration / volume fraction field type.
      typedef typename Propagator<D>::CField CField;

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
      * Read body of input parameters block (without opening and closing lines).
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
      */
      void readCommands();

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

      //@}
      /// \name Chemical Potential Fields (W Fields)
      //@{

      /**
      * Get an array of chemical potential fields, in a basis.
      *
      * This function returns an array in which each element is an
      * array containing the coefficients of the chemical potential
      * field (w field) in a symmetry-adapted basis for one monomer 
      * type. The array capacity is the number of monomer types.
      */
      DArray<RDField <D> >& wFields();

      /**
      * Get chemical potential field for one monomer type, in a basis.
      *
      * This function returns an array containing coefficients of 
      * the chemical potential field (w field) in a symmetry-adapted
      * basis for a specified monomer type.
      *
      * \param monomerId integer monomer type index
      */
      RDField<D>& wField(int monomerId);

      /**
      * Get array of chemical potential fields, on an r-space grid.
      *
      * This function returns an array in which each element is a
      * WField object containing values of the chemical potential field 
      * (w field) on a regular grid for one monomer type. The array 
      * capacity is the number of monomer types.
      */
      DArray<WField>& wFieldsRGrid();

      /**
      * Get the chemical potential field for one monomer type, on a grid.
      *
      * \param monomerId integer monomer type index
      */
      WField& wFieldRGrid(int monomerId);

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
      DArray<RDField <D> >& cFields();

      /**
      * Get the concentration field for one monomer type, in a basis.
      *
      * This function returns an array containing the coefficients of 
      * the monomer concentration / volume fraction field (c field) 
      * for a specific monomer type. 
      *
      * \param monomerId integer monomer type index
      */
      RDField<D>& cField(int monomerId);

      /**
      * Get array of all concentration fields (c fields), on a grid.
      *
      * This function returns an array in which each element is the
      * monomer concentration field for one monomer type on a regular 
      * grid (an r-grid). 
      */
      DArray<CField>& cFieldsRGrid();

      /**
      * Get the concentration (c field) for one monomer type, on a grid.
      *
      * \param monomerId integer monomer type index
      */
      CField& cFieldRGrid(int monomerId);

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
      * Get spatial discretization mesh by reference.
      */
      Mesh<D>& mesh();

      /**
      * Get crystal unitCell (i.e., lattice type and parameters) by reference.
      */
      UnitCell<D>& unitCell();

      /**
      * Get interaction (i.e., excess free energy model) by reference.
      */
      ChiInteraction& interaction();
	  
      /**
      * Get the Iterator by reference.
      */
      //temporarily changed to allow testing on member functions
      AmIterator<D>& iterator();

      /**
      * Get basis object by reference.
      */
      Basis<D>& basis();

      /**
      * Get container for wavevector data.
      */
      WaveList<D>& wavelist();

      /**
      * Get associated FieldIo object.
      */
      FieldIo<D>& fieldIo();

      /**
      * Get associated FFT objecti by reference.
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

      //@}
      /// \name Accessors (return values)
      //@{
      
      /** 
      * Get the group name string.
      */  
      std::string groupName();

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

      //@}

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
      * Spatial discretization mesh.
      */
      Mesh<D> mesh_;

      /**
      * group name.
      */
      std::string groupName_;

      /**
      * Crystallographic unit cell (type and dimensions).
      */
      UnitCell<D> unitCell_;

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
      AmIterator<D>* iteratorPtr_;

      /**
      * Pointer to a Basis object
      */
      Basis<D>* basisPtr_;
      
     
      /**
      * Container for wavevector data.   
      */ 
      WaveList<D>* wavelistPtr_;

      /** 
      * FieldIo object for field input/output operations
      */
      FieldIo<D> fieldIo_;

      /**
      * FFT object to be used by iterator
      */
      FFT<D> fft_;

      /**
      * Pointer to an Sweep object
      */
      Sweep* sweepPtr_;

      /**
      * Pointer to SweepFactory object
      */
      SweepFactory* sweepFactoryPtr_;

      /**
      * Array of chemical potential fields for monomer types.
      *
      * Indexed by monomer typeId, size = nMonomer.
      */
      DArray<RDField<D> > wFields_;

      /**
      * Array of chemical potential fields for monomer types.
      *
      * Indexed by monomer typeId, size = nMonomer.
      */
      DArray<WField> wFieldsRGrid_;

      /**
      * work space for chemical potential fields
      *
      */
      DArray<RDFieldDft<D> > wFieldsKGrid_;

      /**
      * Array of concentration fields for monomer types.
      *
      * Indexed by monomer typeId, size = nMonomer.
      */
      DArray< RDField<D> > cFields_;

      DArray<CField> cFieldsRGrid_;

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
      * Has the Mesh been initialized?
      */
      bool hasMesh_;

      /**
      * Has the UnitCell been initialized?
      */
      bool hasUnitCell_;

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
      * Does this system have a Sweep object?
      */
      bool hasSweep_;

      IntVec<D> kMeshDimensions_;

      RDField<D> workArray;

      cufftReal* d_kernelWorkSpace_;

      cufftReal* kernelWorkSpace_;

      /**
      * Allocate memory for fields (private)
      */
      void allocate();

      /**
      * Initialize Homogeneous::Mixture object (private).
      */
      void initHomogeneous();

      /**
      * Compute inner product of two RDField fields (private, on GPU).
      */
      cufftReal innerProduct(const RDField<D>& a, const RDField<D>& b, int size);

      /**
      * Compute reduction of an RDField (private, on GPU).
      */
      cufftReal reductionH(const RDField<D>& a, int size);

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

   // Get the associated Mixture<D> object.
   template <int D>
   inline Mixture<D>& System<D>::mixture()
   { return mixture_; }

   // Get the Mesh<D>
   template <int D>
   inline Mesh<D>& System<D>::mesh()
   { return mesh_; }

   // Get the UnitCell<D>.
   template <int D>
   inline UnitCell<D>& System<D>::unitCell()
   { return unitCell_; }

   // Get the Interaction (excess free energy model).
   template <int D>
   inline ChiInteraction& System<D>::interaction()
   {
      UTIL_ASSERT(interactionPtr_);
      return *interactionPtr_;
   }

   // Get the FFT<D> object.
   template <int D>
   inline FFT<D>& System<D>::fft()
   {  return fft_; }

   // Get the Basis<D> object.
   template <int D>
   inline Basis<D>& System<D>::basis()
   {
      UTIL_ASSERT(basisPtr_);
      return *basisPtr_;
   }

   // Get the FieldIo<D> object.
   template <int D>
   inline FieldIo<D>& System<D>::fieldIo()
   {  return fieldIo_; }

   // Get the Iterator.
   template <int D>
   inline AmIterator<D>& System<D>::iterator()
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

   // Get all monomer chemical potential (w field), in a basis.
   template <int D>
   inline
   DArray<RDField<D> >& System<D>::wFields()
   {  return wFields_; }

   // Get a single monomer chemical potential (w field), in a basis.
   template <int D>
   inline
   RDField<D>& System<D>::wField(int id)
   {  return wFields_[id]; }

   // Get all monomer excess chemical potential fields, on a grid.
   template <int D>
   inline 
   DArray< typename System<D>::WField >& System<D>::wFieldsRGrid()
   {  return wFieldsRGrid_; }

   // Get a single monomer hemical potential field, on a grid.
   template <int D>
   inline 
   typename System<D>::WField& System<D>::wFieldRGrid(int id)
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
   DArray<RDField<D> >& System<D>::cFields()
   {  return cFields_; }

   // Get a single monomer concentration field, in a basis.
   template <int D>
   inline
   RDField<D>& System<D>::cField(int id)
   {  return cFields_[id]; }

   // Get all monomer concentration fields, on a grid.
   template <int D>
   inline
   DArray< typename System<D>::CField >& System<D>::cFieldsRGrid()
   {  return cFieldsRGrid_; }

   // Get a single monomer concentration field, on a grid.
   template <int D>
   inline typename System<D>::CField& System<D>::cFieldRGrid(int id)
   {  return cFieldsRGrid_[id]; }

   template <int D>
   inline
   DArray<RDFieldDft<D> >& System<D>::cFieldsKGrid()
   {  return cFieldsKGrid_; }

   template <int D>
   inline
   RDFieldDft<D>& System<D>::cFieldKGrid(int id)
   {  return cFieldsKGrid_[id]; }

   // Get group name string.
   template <int D>
   inline std::string System<D>::groupName()
   { return groupName_; }

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

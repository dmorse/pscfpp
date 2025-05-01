#ifndef PRDC_DOMAIN_REAL_H
#define PRDC_DOMAIN_REAL_H

/*
* PSCF - Polymer Self-Consistent Field Theory
*
* Copyright 2015 - 2025, The Regents of the University of Minnesota
* Distributed under the terms of the GNU General Public License.
*/

#include <util/param/ParamComposite.h>    // base class

#include <prdc/cpu/WaveList.h>            // member
#include <prdc/cpu/FFT.h>                 // member
#include <prdc/crystal/Basis.h>           // member
#include <prdc/crystal/SpaceGroup.h>      // member
#include <prdc/crystal/UnitCell.h>        // member

#include <pscf/mesh/Mesh.h>               // member

#include <string>

namespace Pscf {
namespace Prdc {

   using namespace Util;
   using namespace Pscf::Prdc;
   using namespace Pscf::Prdc::Cpu;

   /**
   * Spatial domain and spatial discretization for a periodic structure.
   *
   * A DomainReal has (among other components):
   *
   *    - a Pscf::Mesh spatial discretization mesh
   *    - a Prdc::UnitCell crystallographic unit cell
   *    - a Prdc::SpaceGroup crystallographic space group
   *    - a Prdc::Basis symmetry-adapated Fourier basis
   *    - a Prdc::Cpu::FFT Fast Fourier Transform
   *    - a Prdc::Cpu::WaveList container for wavevector properties
   *    - a lattice system enum (type Prdc::UnitCell\<D\>::LatticeSystem)
   *    - a groupName string
   *
   * Note: Class names Pscf::Mesh, Prdc::UnitCell, etc. mentioned above are
   * actually all class templates with a template parameter D. Actual class
   * names are Mesh \<D\>, Prdc::UnitCell \<D\>, etc. with D=1, 2, or 3.
   *
   * \ingroup Prdc_Field_Module
   */
   template <int D>
   class DomainReal : public ParamComposite
   {

   public:


      /**
      * Constructor.
      */
      DomainReal();

      /**
      * Destructor.
      */
      ~DomainReal();

      /// \name Initialization and Mutators
      ///@{

      /**
      * Read body of parameter block (without opening and closing lines).
      *
      * \param in input parameter stream
      */
      virtual void readParameters(std::istream& in);

      /**
      * Read initialization data from header of an r-grid field file.
      *
      * This is an alternative to reading the parameter file that is only 
      * used for unit testing.
      *
      * \param in  input parameter stream
      * \param nMonomer  number of monomers in field file (output)
      */
      void readRGridFieldHeader(std::istream& in, int& nMonomer);

      /**
      * Set unit cell by copying another UnitCell<D> object.
      *
      * The lattice system in the unitCell must match any value that was
      * read from the parameter file. This function initializes the basis
      * if needed.
      *
      * \param unitCell new unit cell
      */
      void setUnitCell(UnitCell<D> const & unitCell);

      /**
      * Set unit cell state, given the lattice system and parameters.
      *
      * The "lattice" enumeration value must match any value that was
      * read from the parameter file. This function initializes the basis 
      * if needed. 
      *
      * \param lattice  lattice system
      * \param parameters array of unit cell parameters
      */
      void setUnitCell(typename UnitCell<D>::LatticeSystem lattice,
                       FSArray<double, 6> const & parameters);

      /**
      * Set unit cell parameters.
      *
      * The lattice system must be set to non-null value on entry. The
      * size of the parameters array must match the number expected for
      * the lattice type. The basis is initialized if needed. 
      *
      * \param parameters array of unit cell parameters
      */
      void setUnitCell(FSArray<double, 6> const & parameters);

      /**
      * Construct basis if not done already.
      */
      void makeBasis();

      ///@}
      /// \name Accessors (return objects by reference)
      ///@{

      /**
      * Get the UnitCell by non-const reference.
      */
      UnitCell<D>& unitCell();

      /**
      * Get the UnitCell by const reference.
      */
      UnitCell<D> const & unitCell() const;

      /**
      * Get the Mesh by non-const reference.
      */
      Mesh<D>& mesh();

      /**
      * Get the Mesh by const reference.
      */
      Mesh<D> const & mesh() const;

      /**
      * Get the SpaceGroup by const reference.
      */
      SpaceGroup<D> const & group() const ;

      /**
      * Get the Basis object by non-const reference.
      */
      Basis<D>& basis();

      /**
      * Get the Basis by const reference.
      */
      Basis<D> const & basis() const ;

      /**
      * Get the FFT by non-const reference.
      */
      FFT<D>& fft();

      /**
      * Get the FFT object by non-const reference.
      */
      FFT<D> const & fft() const;

      /**
      * Get the WaveList by non-const reference.
      */
      WaveList<D>& waveList();

      /**
      * Get the WaveList by const reference.
      */
      WaveList<D> const & waveList() const;

      ///@}
      /// \name Accessors (return by value)
      ///@{

      /** 
      * Get the lattice system (enumeration value).
      */  
      typename UnitCell<D>::LatticeSystem lattice() const;

      /** 
      * Get the group name.
      */  
      std::string groupName() const;

      /**
      * Has a space group been identified?
      */ 
      bool hasGroup() const;

      /**
      * Has a symmetry-adapted Fourier basis been initialized?
      */ 
      bool hasBasis() const;

      ///@}

   private:

      // Private member variables

      /**
      * Crystallographic unit cell (crystal system and cell parameters).
      */
      UnitCell<D> unitCell_;

      /**
      * Spatial discretization mesh.
      */
      Mesh<D> mesh_;

      /**
      * SpaceGroup object
      */
      SpaceGroup<D> group_;

      /**
      * Basis object.
      */
      Basis<D> basis_;

      /**
      * FFT object to be used by solvers.
      */
      FFT<D> fft_;

      /**
      * WaveList object.
      */
      WaveList<D> waveList_;

      /**
      * Lattice system (enumeration value).
      */
      typename UnitCell<D>::LatticeSystem lattice_;

      /**
      * Group name.
      */
      std::string groupName_;

      /**
      * Has a space group been indentified?
      */
      bool hasGroup_;

      /**
      * Has this DomainReal object been initialized?
      */
      bool isInitialized_;

   };

   // Inline member functions

   // Get the UnitCell by non-const reference.
   template <int D>
   inline UnitCell<D>& DomainReal<D>::unitCell()
   {  return unitCell_; }

   // Get the UnitCell by const reference.
   template <int D>
   inline UnitCell<D> const & DomainReal<D>::unitCell() const
   {  return unitCell_; }

   // Get the Mesh by non-const reference.
   template <int D>
   inline Mesh<D>& DomainReal<D>::mesh()
   {  return mesh_; }

   // Get the Mesh by const reference.
   template <int D>
   inline Mesh<D> const & DomainReal<D>::mesh() const
   {  return mesh_; }

   // Get the SpaceGroup by const reference.
   template <int D>
   inline SpaceGroup<D> const & DomainReal<D>::group() const
   {  return group_; }

   // Get the Basis by non-const reference.
   template <int D>
   inline Basis<D>& DomainReal<D>::basis()
   {  return basis_; }

   // Get the Basis by const reference.
   template <int D>
   inline Basis<D> const & DomainReal<D>::basis() const
   {  return basis_; }

   // Get the FFT by non-const reference.
   template <int D>
   inline FFT<D>& DomainReal<D>::fft()
   {  return fft_; }

   // Get the FFT by const reference.
   template <int D>
   inline FFT<D> const & DomainReal<D>::fft() const
   {  return fft_; }

   // Get the WaveList by non-const reference.
   template <int D>
   inline WaveList<D>& DomainReal<D>::waveList()
   {  return waveList_; }

   // Get the WaveList by const reference.
   template <int D>
   inline WaveList<D> const & DomainReal<D>::waveList() const
   {  return waveList_; }

   // Get the lattice system enumeration value
   template <int D>
   inline 
   typename UnitCell<D>::LatticeSystem DomainReal<D>::lattice() 
   const
   {  return lattice_; }

   // Get the groupName string.
   template <int D>
   inline std::string DomainReal<D>::groupName() const
   {  return groupName_; }

   // Has a space group been identified?
   template <int D>
   inline bool DomainReal<D>::hasGroup() const
   {  return hasGroup_; }

   // Has a symmetry-adapted Fourier basis been initialized ?
   template <int D>
   inline bool DomainReal<D>::hasBasis() const
   {  return basis_.isInitialized(); }

   #ifndef PRDC_DOMAIN_REAL_TPP
   // Suppress implicit instantiation
   extern template class DomainReal<1>;
   extern template class DomainReal<2>;
   extern template class DomainReal<3>;
   #endif

} // namespace Prdc
} // namespace Pscf
#endif

#ifndef PRDC_DOMAIN_REAL_H
#define PRDC_DOMAIN_REAL_H

/*
* PSCF - Polymer Self-Consistent Field Theory
*
* Copyright 2015 - 2025, The Regents of the University of Minnesota
* Distributed under the terms of the GNU General Public License.
*/

#include <util/param/ParamComposite.h>    // base class

#include <prdc/crystal/Basis.h>           // member
#include <prdc/crystal/SpaceGroup.h>      // member
#include <prdc/crystal/UnitCell.h>        // member
#include <pscf/mesh/Mesh.h>               // member
#include <string>                         // member (groupName)

// Forward declaration
namespace Util{
   class FileMaster;
}

namespace Pscf {
namespace Prdc {

   using namespace Util;
   using namespace Pscf::Prdc;

   /**
   * Spatial domain for a periodic structure with real fields.
   *
   * Partial specializations of the DomainReal class template are used
   * as base classes for Rpc::Domain<int D> and Rpg::Domain<int D>.
   *
   * Template Parameters:
   *
   *   - D    : integer dimension of space (D=1, 2, or 3)
   *   - FFT  : Fast Fourier transform calculator type, e.g., FFT<D>
   *   - WLT  : WaveList container type, e.g., WaveList<D>
   *   - FIT  : FieldIo class for field operations, e.g., FieldIo<D>
   *
   * A DomainReal template instance has:
   *
   *  - a Mesh spatial discretization mesh
   *  - a UnitCell crystallographic unit cell
   *  - a SpaceGroup crystallographic space group
   *  - a Basis symmetry-adapated Fourier basis
   *  - a FFT Fast Fourier Transform calculator (class FFT)
   *  - a WaveList container for wavevector properties (class WLT)
   *  - a FieldIo object for field IO & conversion operations (class FIT)
   *  - a lattice system enum (type Prdc::UnitCell\<D\>::LatticeSystem)
   *  - a groupName string
   *
   * Note: Class names Pscf::Mesh, Prdc::UnitCell, etc. mentioned above are
   * actually class templates with a template parameter D. Actual class 
   * names are Mesh \<D\>, Prdc::UnitCell \<D\>, etc. with D=1, 2, or 3.
   *
   * \ingroup Prdc_Field_Module
   */
   template <int D, class FFT, class WLT, class FIT>
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

      /// \name Initialization
      ///@{

      /**
      * Create association with a FileMaster, needed by FieldIo.
      *
      * \param fileMaster associated FileMaster object.
      */
      void setFileMaster(FileMaster& fileMaster);

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
      * Construct basis if not done already.
      */
      void makeBasis();

      ///@}
      /// \name Unit Cell Modifiers
      ///@{

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

      ///@}
      /// \name Accessors (return component objects by reference)
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
      Basis<D> const & basis() const;

      /**
      * Get the FFT by non-const reference.
      */
      FFT& fft();

      /**
      * Get the FFT object by non-const reference.
      */
      FFT const & fft() const;

      /**
      * Get the WaveList by non-const reference.
      */
      WLT& waveList();

      /**
      * Get the WaveList by const reference.
      */
      WLT const & waveList() const;

      /**
      * Get the FieldIo by non-const reference.
      */
      FIT& fieldIo();

      /**
      * Get the FieldIo by const reference.
      */
      FIT const & fieldIo() const;

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
      /// \name Crystallography Information Output
      ///@{

      /**
      * Output information about waves.
      *
      * This function opens a file with the specified filename, calls
      * Basis<D>::outputWaves, and closes the file before returning.
      *
      * \param filename name of output file
      */
      void writeWaves(std::string const & filename) const;

      /**
      * Output information about stars and symmetrized basis functions.
      *
      * This function opens a file with the specified filename, calls
      * Basis<D>::outputStars, and closes the file before returning.
      *
      * \param filename name of output file
      */
      void writeStars(std::string const & filename) const;

      /**
      * Output all elements of the space group.
      *
      * \param filename name of output file
      */
      void writeGroup(std::string const & filename) const;

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
      FFT fft_;

      /**
      * WaveList object.
      */
      WLT waveList_;

      /**
      * FieldIo object for field input/output operations
      */
      FIT fieldIo_;

      /**
      * Lattice system (enumeration value).
      */
      typename UnitCell<D>::LatticeSystem lattice_;

      /**
      * Group name.
      */
      std::string groupName_;

      /**
      * Pointer to associated FileMaster.
      */
      FileMaster* fileMasterPtr_;

      /**
      * Has a space group been indentified?
      */
      bool hasGroup_;

      /**
      * Has this DomainReal object been initialized?
      */
      bool isInitialized_;

      /*
      * Get FileMaster as const reference.
      */
      FileMaster const & fileMaster() const
      {
         UTIL_CHECK(fileMasterPtr_);
         return * fileMasterPtr_;
      }

      // Members of parent class with non-dependent names
      using ParamComposite::read;
      using ParamComposite::readOptional;

   };

   // Inline member functions

   // Get the UnitCell by non-const reference.
   template <int D, class FFT, class WLT, class FIT>
   inline UnitCell<D>& DomainReal<D,FFT,WLT,FIT>::unitCell()
   {  return unitCell_; }

   // Get the UnitCell by const reference.
   template <int D, class FFT, class WLT, class FIT>
   inline UnitCell<D> const & DomainReal<D,FFT,WLT,FIT>::unitCell() const
   {  return unitCell_; }

   // Get the Mesh by non-const reference.
   template <int D, class FFT, class WLT, class FIT>
   inline Mesh<D>& DomainReal<D,FFT,WLT,FIT>::mesh()
   {  return mesh_; }

   // Get the Mesh by const reference.
   template <int D, class FFT, class WLT, class FIT>
   inline Mesh<D> const & DomainReal<D,FFT,WLT,FIT>::mesh() const
   {  return mesh_; }

   // Get the SpaceGroup by const reference.
   template <int D, class FFT, class WLT, class FIT>
   inline SpaceGroup<D> const & DomainReal<D,FFT,WLT,FIT>::group() const
   {  return group_; }

   // Get the Basis by non-const reference.
   template <int D, class FFT, class WLT, class FIT>
   inline Basis<D>& DomainReal<D,FFT,WLT,FIT>::basis()
   {  return basis_; }

   // Get the Basis by const reference.
   template <int D, class FFT, class WLT, class FIT>
   inline Basis<D> const & DomainReal<D,FFT,WLT,FIT>::basis() const
   {  return basis_; }

   // Get the FFT by non-const reference.
   template <int D, class FFT, class WLT, class FIT>
   inline FFT& DomainReal<D,FFT,WLT,FIT>::fft()
   {  return fft_; }

   // Get the FFT by const reference.
   template <int D, class FFT, class WLT, class FIT>
   inline FFT const & DomainReal<D,FFT,WLT,FIT>::fft() const
   {  return fft_; }

   // Get the WaveList by non-const reference.
   template <int D, class FFT, class WLT, class FIT>
   inline WLT& DomainReal<D,FFT,WLT,FIT>::waveList()
   {  return waveList_; }

   // Get the WaveList by const reference.
   template <int D, class FFT, class WLT, class FIT>
   inline WLT const & DomainReal<D,FFT,WLT,FIT>::waveList() const
   {  return waveList_; }

   // Get the FieldIo by const reference.
   template <int D, class FFT, class WLT, class FIT>
   inline FIT& DomainReal<D,FFT,WLT,FIT>::fieldIo()
   {  return fieldIo_; }

   // Get the FieldIo by const reference.
   template <int D, class FFT, class WLT, class FIT>
   inline FIT const & DomainReal<D,FFT,WLT,FIT>::fieldIo() const
   {  return fieldIo_; }

   // Get the lattice system enumeration value
   template <int D, class FFT, class WLT, class FIT>
   inline 
   typename UnitCell<D>::LatticeSystem DomainReal<D,FFT,WLT,FIT>::lattice() 
   const
   {  return lattice_; }

   // Get the groupName string.
   template <int D, class FFT, class WLT, class FIT>
   inline std::string DomainReal<D,FFT,WLT,FIT>::groupName() const
   {  return groupName_; }

   // Has a space group been identified?
   template <int D, class FFT, class WLT, class FIT>
   inline bool DomainReal<D,FFT,WLT,FIT>::hasGroup() const
   {  return hasGroup_; }

   // Has a symmetry-adapted Fourier basis been initialized ?
   template <int D, class FFT, class WLT, class FIT>
   inline bool DomainReal<D,FFT,WLT,FIT>::hasBasis() const
   {  return basis_.isInitialized(); }

} // namespace Prdc
} // namespace Pscf
#endif

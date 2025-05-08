#ifndef RPG_DOMAIN_H
#define RPG_DOMAIN_H

/*
* PSCF - Polymer Self-Consistent Field Theory
*
* Copyright 2015 - 2025, The Regents of the University of Minnesota
* Distributed under the terms of the GNU General Public License.
*/

#include <util/param/ParamComposite.h>     // base class

#include <rpg/field/FieldIo.h>             // member

#include <prdc/cuda/WaveList.h>            // member
#include <prdc/cuda/FFT.h>                 // member
#include <prdc/crystal/Basis.h>            // member
#include <prdc/crystal/SpaceGroup.h>       // member
#include <prdc/crystal/UnitCell.h>         // member

#include <pscf/mesh/Mesh.h>                // member

#include <string>

namespace Pscf {
namespace Rpg {

   using namespace Util;
   using namespace Pscf::Prdc;
   using namespace Pscf::Prdc::Cuda;

   /**
   * Spatial domain and spatial discretization for a periodic structure.
   *
   * A Domain has (among other components):
   *
   *    - a Mesh
   *    - a UnitCell
   *    - a SpaceGroup
   *    - a Basis
   *    - an Fft 
   *    - a WaveList
   *    - a FieldIo
   *    - a lattice system enumeration value
   *    - a groupName string
   *
   * \ingroup Rpg_Field_Module
   */
   template <int D>
   class Domain : public ParamComposite
   {

   public:

      /**
      * Constructor.
      */
      Domain();

      /**
      * Destructor.
      */
      ~Domain();

      /// \name Initialization and Mutators
      //@{

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
      * Set the unit cell by copying another UnitCell<D> object.
      *
      * The lattice system in the unitCell must match any value that was
      * read from the parameter file. This function initialized the basis
      * if needed.
      *
      * \param unitCell new unit cell
      */
      void setUnitCell(UnitCell<D> const & unitCell);

      /**
      * Set the unit cell state, given the lattice system and parameters.
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
      * The lattice system must be set to a non-null value on entry. The
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

      //@}
      /// \name Accessors (return objects by reference)
      //@{

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
      SpaceGroup<D> const & group() const;

      /**
      * Get the Basis by non-const reference.
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
      * Get the FFT by const reference.
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

      /**
      * Get the FieldIo by non-const reference.
      */
      FieldIo<D>& fieldIo();

      /**
      * Get the FieldIo by const reference.
      */
      FieldIo<D> const & fieldIo() const;

      //@}
      /// \name Accessors (return by value)
      //@{

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

      //@}
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
      * SpaceGroup object.
      */
      SpaceGroup<D> group_;

      /**
      * Basis object
      */
      Basis<D> basis_;

      /**
      * FFT object.
      */
      FFT<D> fft_;

      /**
      * WaveList object.
      */
      WaveList<D> waveList_;

      /**
      * FieldIo object for field input/output operations
      */
      FieldIo<D> fieldIo_;

      /**
      * Lattice system (enumeration value).
      */
      typename UnitCell<D>::LatticeSystem lattice_;

      /**
      * Group name.
      */
      std::string groupName_;

      /**
      * Pointer to associated FileMaster
      */ 
      FileMaster* fileMasterPtr_;

      /**
      * Has a space group been indentified?
      */
      bool hasGroup_;

      /**
      * Has this Domain object been initialized?
      */
      bool isInitialized_;

      /**
      * Get the FileMaster by reference (private).
      */
      FileMaster const & fileMaster() const
      {
         UTIL_CHECK(fileMasterPtr_);
         return *fileMasterPtr_; 
      }

   };

   // Inline member functions

   // Get the UnitCell by non-const reference.
   template <int D>
   inline UnitCell<D>& Domain<D>::unitCell()
   {  return unitCell_; }

   // Get the UnitCell by const reference.
   template <int D>
   inline UnitCell<D> const & Domain<D>::unitCell() const
   {  return unitCell_; }

   // Get the Mesh by non-const reference.
   template <int D>
   inline Mesh<D>& Domain<D>::mesh()
   {  return mesh_; }

   // Get the Mesh by const reference.
   template <int D>
   inline Mesh<D> const & Domain<D>::mesh() const
   {  return mesh_; }

   // Get the SpaceGroup by const reference.
   template <int D>
   inline SpaceGroup<D> const & Domain<D>::group() const
   {  return group_; }

   // Get the Basis by non-const reference.
   template <int D>
   inline Basis<D>& Domain<D>::basis()
   {  return basis_; }

   // Get the Basis by const reference.
   template <int D>
   inline Basis<D> const & Domain<D>::basis() const
   {  return basis_; }

   // Get the FFT by non-const reference.
   template <int D>
   inline FFT<D>& Domain<D>::fft()
   {  return fft_; }

   // Get the FFT by const reference.
   template <int D>
   inline FFT<D> const & Domain<D>::fft() const
   {  return fft_; }

   // Get the WaveList by non-const reference.
   template <int D>
   inline WaveList<D>& Domain<D>::waveList()
   {  return waveList_; }

   // Get the WaveList by const reference.
   template <int D>
   inline WaveList<D> const & Domain<D>::waveList() const
   {  return waveList_; }

   // Get the FieldIo by non-const reference.
   template <int D>
   inline FieldIo<D>& Domain<D>::fieldIo()
   {  return fieldIo_; }

   // Get the FieldIo by const reference.
   template <int D>
   inline FieldIo<D> const & Domain<D>::fieldIo() const
   {  return fieldIo_; }

   // Get the lattice system enumeration value.
   template <int D>
   inline 
   typename UnitCell<D>::LatticeSystem Domain<D>::lattice() const
   {  return lattice_; }

   // Get the groupName string.
   template <int D>
   inline std::string Domain<D>::groupName() const
   {  return groupName_; }

   // Has a space group been identified?
   template <int D>
   inline bool Domain<D>::hasGroup() const
   {  return hasGroup_; }

   // Has a symmetry-adapted Fourier basis been initialized ?
   template <int D>
   inline bool Domain<D>::hasBasis() const
   {  return basis_.isInitialized(); }

   #ifndef RPG_DOMAIN_TPP
   // Suppress implicit instantiation
   extern template class Domain<1>;
   extern template class Domain<2>;
   extern template class Domain<3>;
   #endif

} // namespace Rpg
} // namespace Pscf
#endif

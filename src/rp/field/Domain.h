#ifndef RP_DOMAIN_H
#define RP_DOMAIN_H

/*
* PSCF - Polymer Self-Consistent Field
*
* Copyright 2015 - 2025, The Regents of the University of Minnesota
* Distributed under the terms of the GNU General Public License.
*/

#include <util/param/ParamComposite.h>    // base class
#include <prdc/crystal/UnitCell.h>        // member
#include <pscf/mesh/Mesh.h>               // member
#include <string>                         // member (groupName)

// Forward declarations
namespace Util {
   class FileMaster;
   template <typename T> class Signal;
   template <> class Signal<void>;
}
namespace Pscf {
   namespace Prdc {
      template <int D> class SpaceGroup;
      template <int D> class Basis;
      template <int D, class T> class FFT;
      template <int D, class T> class WaveList;
   }
   namespace Rp {
      template <int D, class T> class FieldIo;
   }
}

namespace Pscf {
namespace Rp {

   using namespace Util;
   using namespace Pscf::Prdc;

   /**
   * Spatial domain for a periodic structure with real fields.
   *
   * A Domain instance has:
   *
   *  - a Mesh spatial discretization mesh
   *  - a UnitCell crystallographic unit cell
   *  - a SpaceGroup crystallographic space group
   *  - a Basis symmetry-adapated Fourier basis
   *  - a FFT Fast Fourier Transform calculator (T::FFT)
   *  - a WaveList container for wavevector properties (T::WaveList)
   *  - a FieldIo object for field IO & conversion (T::FieldIo)
   *  - a lattice system enum (type Prdc::UnitCell\<D\>::LatticeSystem)
   *  - a groupName string
   *
   * Note: Class names Pscf::Mesh, Prdc::UnitCell, etc. mentioned above
   * refer to class templates with an integer template parameter D. Actual
   * class names are Mesh \<D\>, Prdc::UnitCell \<D\>, etc. with D=1, 2,
   * or 3.
   *
   * <b> Template Parameters </b>:
   *
   *   - D  : integer dimension of space (D=1, 2, or 3)
   *   - T  : Types class (CppTp<D> or CudaTp<D>)
   *
   * <b> Subclasses </b>: Partial specializations of the Domain class
   * template are used as base classes for classes Rpc::Domain \<D\> and
   * Rpg::Domain \<D\>.
   *
   * \ingroup Rp_Field_Module
   */
   template <int D, class T>
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
      /// \name Accessors (return component objects by reference)
      ///@{

      /**
      * Get the Mesh by non-const reference.
      */
      Mesh<D>& mesh();

      /**
      * Get the Mesh by const reference.
      */
      Mesh<D> const & mesh() const;

      /**
      * Get the UnitCell by non-const reference.
      */
      UnitCell<D>& unitCell();

      /**
      * Get the UnitCell by const reference.
      */
      UnitCell<D> const & unitCell() const;

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
      FFT<D,T>& fft();

      /**
      * Get the FFT object by non-const reference.
      */
      FFT<D,T> const & fft() const;

      /**
      * Get the WaveList by non-const reference.
      */
      WaveList<D,T>& waveList();

      /**
      * Get the WaveList by const reference.
      */
      WaveList<D,T> const & waveList() const;

      /**
      * Get the FieldIo by non-const reference.
      */
      FieldIo<D,T>& fieldIo();

      /**
      * Get the FieldIo by const reference.
      */
      FieldIo<D,T> const & fieldIo() const;

      ///@}
      /// \name Accessors (return by value)
      ///@{

      /**
      * Get the lattice type (enumeration value).
      */
      typename UnitCell<D>::LatticeSystem lattice() const;

      /**
      * Get the group name.
      */
      std::string groupName() const;

      /**
      * Has a space group been declared?
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
      * Spatial discretization mesh.
      */
      Mesh<D> mesh_;

      /**
      * Crystallographic unit cell (crystal system and cell parameters).
      */
      UnitCell<D> unitCell_;

      /**
      * Lattice system (enumeration value).
      */
      typename UnitCell<D>::LatticeSystem lattice_;

      /**
      * Group name.
      */
      std::string groupName_;

      // Pointers to owned and associated objects.

      /**
      * Pointer to a SpaceGroup object (owned).
      */
      SpaceGroup<D>* groupPtr_;

      /**
      * Pointer to a Basis object (owned).
      */
      Basis<D>* basisPtr_;

      /**
      * Pointer to a FFT (Fast Fourier Transform) object (owned).
      */
      FFT<D,T>* fftPtr_;

      /**
      * Pointer to a WaveList object (owned).
      */
      WaveList<D,T>* waveListPtr_;

      /**
      * Pointer to a FieldIo object (owned).
      */
      FieldIo<D,T>* fieldIoPtr_;

      /**
      * Pointer to a Signal object (owned).
      */
      Signal<void>* signalPtr_;

      /**
      * Pointer to associated FileMaster.
      */
      FileMaster* fileMasterPtr_;

      // Boolean flags

      /**
      * Has a space group been indentified?
      */
      bool hasGroup_;

      /**
      * Has this Domain object been initialized?
      */
      bool isInitialized_;

      // Private member function

      /*
      * Get FileMaster as const reference.
      */
      FileMaster const & fileMaster() const;

   };

   // Public inline member functions

   // Get the UnitCell by non-const reference.
   template <int D, class T> inline 
   UnitCell<D>& Domain<D,T>::unitCell()
   {  return unitCell_; }

   // Get the UnitCell by const reference.
   template <int D, class T> inline 
   UnitCell<D> const & Domain<D,T>::unitCell() const
   {  return unitCell_; }

   // Get the Mesh by non-const reference.
   template <int D, class T> inline 
   Mesh<D>& Domain<D,T>::mesh()
   {  return mesh_; }

   // Get the Mesh by const reference.
   template <int D, class T> inline 
   Mesh<D> const & Domain<D,T>::mesh() const
   {  return mesh_; }

   // Get the SpaceGroup by const reference.
   template <int D, class T> inline 
   SpaceGroup<D> const & Domain<D,T>::group() const
   {  return *groupPtr_; }

   // Get the Basis by non-const reference.
   template <int D, class T> inline 
   Basis<D>& Domain<D,T>::basis()
   {  return *basisPtr_; }

   // Get the Basis by const reference.
   template <int D, class T> inline 
   Basis<D> const & Domain<D,T>::basis() const
   {  return *basisPtr_; }

   // Get the FFT by non-const reference.
   template <int D, class T> inline 
   FFT<D,T>& Domain<D,T>::fft()
   {  return *fftPtr_; }

   // Get the FFT by const reference.
   template <int D, class T> inline 
   FFT<D,T> const & Domain<D,T>::fft() const
   {  return *fftPtr_; }

   // Get the WaveList by non-const reference.
   template <int D, class T> inline 
   WaveList<D,T>& Domain<D,T>::waveList()
   {  return *waveListPtr_; }

   // Get the WaveList by const reference.
   template <int D, class T> inline 
   WaveList<D,T> const & Domain<D,T>::waveList() const
   {  return *waveListPtr_; }

   // Get the FieldIo by const reference.
   template <int D, class T> inline
   FieldIo<D,T>& Domain<D,T>::fieldIo()
   {  return *fieldIoPtr_; }

   // Get the FieldIo by const reference.
   template <int D, class T> inline
   FieldIo<D,T> const & Domain<D,T>::fieldIo() const
   {  return *fieldIoPtr_; }

   // Get the lattice type enumeration value.
   template <int D, class T> inline
   typename UnitCell<D>::LatticeSystem Domain<D,T>::lattice()
   const
   {  return lattice_; }

   // Get the groupName string identifier.
   template <int D, class T> inline 
   std::string Domain<D,T>::groupName() const
   {  return groupName_; }

   // Has a space group been identified?
   template <int D, class T> inline 
   bool Domain<D,T>::hasGroup() const
   {  return hasGroup_; }

   // Private inline member function

   // Get FileMaster as const reference.
   template <int D, class T> inline
   FileMaster const & Domain<D,T>::fileMaster() const
   {
      UTIL_CHECK(fileMasterPtr_);
      return * fileMasterPtr_;
   }

} // namespace Rp
} // namespace Pscf
#endif

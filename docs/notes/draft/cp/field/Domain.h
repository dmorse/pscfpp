#ifndef PRDC_CL_DOMAIN_H
#define PRDC_CL_DOMAIN_H

/*
* PSCF - Polymer Self-Consistent Field
*
* Copyright 2015 - 2025, The Regents of the University of Minnesota
* Distributed under the terms of the GNU General Public License.
*/

#include <util/param/ParamComposite.h>    // base class
#include <prdc/crystal/UnitCell.h>        // member
#include <pscf/mesh/Mesh.h>               // member

// Forward declarations
namespace Util {
   class FileMaster;
   template <typename T> class Signal;
   template <> class Signal<void>;
}
namespace Pscf {
   namespace Prdc {
      template <int D, class T> class WaveList;
   }
}

namespace Pscf {
namespace Cp {

   using namespace Util;
   using namespace Prdc;

   /**
   * Spatial domain for a periodic structure with complex fields.
   *
   * A Domain template instance has:
   *
   *  - a Mesh spatial discretization mesh
   *  - a UnitCell crystallographic unit cell
   *  - a FFT Fast Fourier Transform calculator (class FFT)
   *  - a WaveList container for wavevector properties (class WLT)
   *  - a FieldIo object for field IO & conversion operations (class FIT)
   *  - a lattice system enum (type Prdc::UnitCell\<D\>::LatticeSystem)
   *
   * Note: Class names Pscf::Mesh, Prdc::UnitCell, etc. mentioned above
   * refer to class templates with an integer template parameter D. Actual
   * class names are Mesh \<D\>, Prdc::UnitCell \<D\>, etc. with D=1, 2,
   * or 3.
   *
   * <b> Template Parameters </b>:
   *
   *   - D    : integer dimension of space (D=1, 2, or 3)
   *   - FFT  : Fast Fourier transform calculator type, e.g., FFT<D>
   *   - WLT  : WaveList container type, e.g., WaveList<D>
   *   - FIT  : FieldIo class for field operations, e.g., FieldIo<D>
   *
   * <b> Subclasses </b>: Partial specializations of the Domain class
   * template are used as base classes for classes Cpc::Domain \<D\> and
   * Cpg::Domain \<D\>.
   *
   * \ingroup Cp_Field_Module
   */
   template <int D, class FFT, class WLT, class FIT>
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
      void readFieldHeader(std::istream& in, int& nMonomer);

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

      /**
      * Get the lattice system (enumeration value).
      */
      typename UnitCell<D>::LatticeSystem lattice() const;

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
      * Pointer to a FFT (Fast Fourier Transform) object (owned).
      */
      FFT* fftPtr_;

      /**
      * Pointer to a WaveList object (owned).
      */
      WLT* waveListPtr_;

      /**
      * Pointer to a FieldIo object (owned).
      */
      FIT* fieldIoPtr_;

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
      * Has this Domain object been initialized?
      */
      bool isInitialized_;

      // Private member function

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
   inline UnitCell<D>& Domain<D,FFT,WLT,FIT>::unitCell()
   {  return unitCell_; }

   // Get the UnitCell by const reference.
   template <int D, class FFT, class WLT, class FIT>
   inline UnitCell<D> const & Domain<D,FFT,WLT,FIT>::unitCell() const
   {  return unitCell_; }

   // Get the Mesh by non-const reference.
   template <int D, class FFT, class WLT, class FIT>
   inline Mesh<D>& Domain<D,FFT,WLT,FIT>::mesh()
   {  return mesh_; }

   // Get the Mesh by const reference.
   template <int D, class FFT, class WLT, class FIT>
   inline Mesh<D> const & Domain<D,FFT,WLT,FIT>::mesh() const
   {  return mesh_; }

   // Get the FFT by non-const reference.
   template <int D, class FFT, class WLT, class FIT>
   inline FFT& Domain<D,FFT,WLT,FIT>::fft()
   {  return *fftPtr_; }

   // Get the FFT by const reference.
   template <int D, class FFT, class WLT, class FIT>
   inline FFT const & Domain<D,FFT,WLT,FIT>::fft() const
   {  return *fftPtr_; }

   // Get the WaveList by non-const reference.
   template <int D, class FFT, class WLT, class FIT>
   inline WLT& Domain<D,FFT,WLT,FIT>::waveList()
   {  return *waveListPtr_; }

   // Get the WaveList by const reference.
   template <int D, class FFT, class WLT, class FIT>
   inline WLT const & Domain<D,FFT,WLT,FIT>::waveList() const
   {  return *waveListPtr_; }

   // Get the FieldIo by const reference.
   template <int D, class FFT, class WLT, class FIT>
   inline FIT& Domain<D,FFT,WLT,FIT>::fieldIo()
   {  return *fieldIoPtr_; }

   // Get the FieldIo by const reference.
   template <int D, class FFT, class WLT, class FIT>
   inline FIT const & Domain<D,FFT,WLT,FIT>::fieldIo() const
   {  return *fieldIoPtr_; }

   // Get the lattice system enumeration value.
   template <int D, class FFT, class WLT, class FIT>
   inline
   typename UnitCell<D>::LatticeSystem Domain<D,FFT,WLT,FIT>::lattice()
   const
   {  return lattice_; }

} // namespace Cp
} // namespace Pscf
#endif

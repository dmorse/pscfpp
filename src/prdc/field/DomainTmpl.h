#ifndef PRDC_DOMAIN_TMPL_H
#define PRDC_DOMAIN_TMPL_H

/*
* PSCF - Polymer Self-Consistent Field
*
* Copyright 2015 - 2025, The Regents of the University of Minnesota
* Distributed under the terms of the GNU General Public License.
*/

#include <util/param/ParamComposite.h>    // base class

#include <prdc/crystal/Basis.h>           // member
#include <prdc/crystal/UnitCell.h>        // member
#include <pscf/mesh/Mesh.h>               // member
#include <string>                         // member (groupName)

// Forward declaration
namespace Util {
   class FileMaster;
   template <typename T> class Signal;
   template <> class Signal<void>;
}
namespace Pscf {
   namespace Prdc {
      template <int D> class SpaceGroup;
   }
}

namespace Pscf {
namespace Prdc {

   using namespace Util;
   using namespace Pscf::Prdc;

   /**
   * Spatial domain for a periodic structure with real fields.
   *
   * A DomainTmpl template instance has:
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
   * actually class templates with an integer template parameter D. Actual 
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
   * <b> Subclasses </b>: Partial specializations of the DomainTmpl class 
   * template are used as base classes for classes Rpc::Domain \<D\> and
   * Rpg::Domain \<D\>.
   *
   * \ingroup Prdc_Field_Module
   */
   template <int D, class FFT, class WLT, class FIT>
   class DomainTmpl : public ParamComposite
   {

   public:

      /**
      * Constructor.
      */
      DomainTmpl();

      /**
      * Destructor.
      */
      ~DomainTmpl();

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

      #if 0
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
      #endif

      /**
      * Lattice system (enumeration value).
      */
      typename UnitCell<D>::LatticeSystem lattice_;

      /**
      * Group name.
      */
      std::string groupName_;

      // Pointers to associated objects

      /**
      * Pointer to a SpaceGroup object
      */
      SpaceGroup<D>* groupPtr_;

      /**
      * Pointer to a Basis object (symmetry-adapted Fourier basis).
      */
      Basis<D>* basisPtr_;

      /**
      * Pointer to a FFT (Fast Fourier Transform) object.
      */
      FFT* fftPtr_;

      /**
      * Pointer to a FieldIo object for field input/output operations.
      */
      WLT* waveListPtr_;

      /**
      * Pointer to a FieldIo object for field input/output operations.
      */
      FIT* fieldIoPtr_;

      /**
      * Pointer to a Signal owned by this DomainTmpl.
      */
      Signal<void>* signalPtr_;

      /**
      * Pointer to associated FileMaster.
      */
      FileMaster* fileMasterPtr_;

      /**
      * Has a space group been indentified?
      */
      bool hasGroup_;

      /**
      * Has this DomainTmpl object been initialized?
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
   inline UnitCell<D>& DomainTmpl<D,FFT,WLT,FIT>::unitCell()
   {  return unitCell_; }

   // Get the UnitCell by const reference.
   template <int D, class FFT, class WLT, class FIT>
   inline UnitCell<D> const & DomainTmpl<D,FFT,WLT,FIT>::unitCell() const
   {  return unitCell_; }

   // Get the Mesh by non-const reference.
   template <int D, class FFT, class WLT, class FIT>
   inline Mesh<D>& DomainTmpl<D,FFT,WLT,FIT>::mesh()
   {  return mesh_; }

   // Get the Mesh by const reference.
   template <int D, class FFT, class WLT, class FIT>
   inline Mesh<D> const & DomainTmpl<D,FFT,WLT,FIT>::mesh() const
   {  return mesh_; }

   // Get the SpaceGroup by const reference.
   template <int D, class FFT, class WLT, class FIT>
   inline SpaceGroup<D> const & DomainTmpl<D,FFT,WLT,FIT>::group() const
   {  return *groupPtr_; }

   // Get the Basis by non-const reference.
   template <int D, class FFT, class WLT, class FIT>
   inline Basis<D>& DomainTmpl<D,FFT,WLT,FIT>::basis()
   {  return *basisPtr_; }

   // Get the Basis by const reference.
   template <int D, class FFT, class WLT, class FIT>
   inline Basis<D> const & DomainTmpl<D,FFT,WLT,FIT>::basis() const
   {  return *basisPtr_; }

   // Get the FFT by non-const reference.
   template <int D, class FFT, class WLT, class FIT>
   inline FFT& DomainTmpl<D,FFT,WLT,FIT>::fft()
   {  return *fftPtr_; }

   // Get the FFT by const reference.
   template <int D, class FFT, class WLT, class FIT>
   inline FFT const & DomainTmpl<D,FFT,WLT,FIT>::fft() const
   {  return *fftPtr_; }

   // Get the WaveList by non-const reference.
   template <int D, class FFT, class WLT, class FIT>
   inline WLT& DomainTmpl<D,FFT,WLT,FIT>::waveList()
   {  return *waveListPtr_; }

   // Get the WaveList by const reference.
   template <int D, class FFT, class WLT, class FIT>
   inline WLT const & DomainTmpl<D,FFT,WLT,FIT>::waveList() const
   {  return *waveListPtr_; }

   // Get the FieldIo by const reference.
   template <int D, class FFT, class WLT, class FIT>
   inline FIT& DomainTmpl<D,FFT,WLT,FIT>::fieldIo()
   {  return *fieldIoPtr_; }

   // Get the FieldIo by const reference.
   template <int D, class FFT, class WLT, class FIT>
   inline FIT const & DomainTmpl<D,FFT,WLT,FIT>::fieldIo() const
   {  return *fieldIoPtr_; }

   // Get the lattice system enumeration value
   template <int D, class FFT, class WLT, class FIT>
   inline 
   typename UnitCell<D>::LatticeSystem DomainTmpl<D,FFT,WLT,FIT>::lattice() 
   const
   {  return lattice_; }

   // Get the groupName string.
   template <int D, class FFT, class WLT, class FIT>
   inline std::string DomainTmpl<D,FFT,WLT,FIT>::groupName() const
   {  return groupName_; }

   // Has a space group been identified?
   template <int D, class FFT, class WLT, class FIT>
   inline bool DomainTmpl<D,FFT,WLT,FIT>::hasGroup() const
   {  return hasGroup_; }

   // Has a symmetry-adapted Fourier basis been initialized ?
   template <int D, class FFT, class WLT, class FIT>
   inline bool DomainTmpl<D,FFT,WLT,FIT>::hasBasis() const
   {  return basis().isInitialized(); }

} // namespace Prdc
} // namespace Pscf
#endif

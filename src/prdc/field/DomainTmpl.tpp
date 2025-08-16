#ifndef PRDC_DOMAIN_TMPL_TPP
#define PRDC_DOMAIN_TMPL_TPP

/*
* PSCF - Polymer Self-Consistent Field
*
* Copyright 2015 - 2025, The Regents of the University of Minnesota
* Distributed under the terms of the GNU General Public License.
*/

#include "DomainTmpl.h"
#include <prdc/field/fieldIoUtil.h>
#include <util/signal/Signal.h>
#include <util/misc/FileMaster.h>

namespace Pscf {
namespace Prdc {

   using namespace Util;

   /*
   * Constructor.
   */
   template <int D, class FFT, class WLT, class FIT>
   DomainTmpl<D,FFT,WLT,FIT>::DomainTmpl()
    : mesh_(),
      unitCell_(),
      group_(),
      basis_(),
      fft_(),
      waveList_(),
      fieldIo_(),
      lattice_(UnitCell<D>::Null),
      groupName_(""),
      signalPtr_(nullptr),
      fileMasterPtr_(nullptr),
      hasGroup_(false),
      isInitialized_(false)
   {
      setClassName("DomainTmpl");
      fieldIo_.associate(mesh_, fft_, lattice_,
                         hasGroup_, groupName_, group_, basis_);

      // Create Signal used by UnitCell, triggered by modification
      signalPtr_ = new Signal<void>();
      unitCell_.setSignal(*signalPtr_);
   }

   /*
   * Destructor.
   */
   template <int D, class FFT, class WLT, class FIT>
   DomainTmpl<D,FFT,WLT,FIT>::~DomainTmpl()
   {
      delete signalPtr_;
   }

   /*
   * Create association with a FileMaster.
   */
   template <int D, class FFT, class WLT, class FIT>
   void DomainTmpl<D,FFT,WLT,FIT>::setFileMaster(FileMaster& fileMaster)
   {
      fileMasterPtr_ = &fileMaster;
      fieldIo_.setFileMaster(fileMaster);
   }

   /*
   * Read parameters and initialize.
   */
   template <int D, class FFT, class WLT, class FIT>
   void DomainTmpl<D,FFT,WLT,FIT>::readParameters(std::istream& in)
   {
      // Preconditions
      UTIL_CHECK(!isInitialized_);
      UTIL_CHECK(fileMasterPtr_);

      // Read computational mesh dimensions (required)
      read(in, "mesh", mesh_);
      UTIL_CHECK(mesh_.size() > 0);
      fft_.setup(mesh_.dimensions());

      // Read lattice system enumeration value (required)
      read(in, "lattice", lattice_);
      unitCell_.set(lattice_);
      UTIL_CHECK(unitCell_.lattice() != UnitCell<D>::Null);
      UTIL_CHECK(unitCell_.nParameter() > 0);

      // Allocate memory for WaveList
      waveList_.allocate(mesh_, unitCell_);

      // Optionally read groupName_ (string identifier for space group)
      hasGroup_ = false;
      bool hasGroupName = false;
      hasGroupName = readOptional(in, "groupName", groupName_).isActive();

      // If groupName_ exists, read and construct group_ (space group)
      if (hasGroupName) {
         // Read group symmetry operations from file
         // An Exception is thrown if groupName_ string is not recognized
         readGroup(groupName_, group_);
         hasGroup_ = true;
      }

      isInitialized_ = true;
   }

   /*
   * Read header of r-grid field to initialize the DomainTmpl.
   *
   * Alternative to parameter file, used only for unit testing.
   */
   template <int D, class FFT, class WLT, class FIT>
   void
   DomainTmpl<D,FFT,WLT,FIT>::readRGridFieldHeader(std::istream& in,
                                                   int& nMonomer)
   {
      // Preconditions - confirm that nothing is initialized
      UTIL_CHECK(!isInitialized_);
      UTIL_CHECK(lattice_ == UnitCell<D>::Null);
      UTIL_CHECK(!unitCell_.isInitialized());
      UTIL_CHECK(!hasGroup_);
      UTIL_CHECK(groupName_ == "");

      // Read common section of standard field header
      int ver1, ver2;
      Pscf::Prdc::readFieldHeader(in, ver1, ver2,
                                  unitCell_, groupName_, nMonomer);

      // Set lattice_ (lattice system identifier)
      lattice_ = unitCell_.lattice();
      UTIL_CHECK(lattice_ != UnitCell<D>::Null);
      UTIL_CHECK(unitCell_.isInitialized());

      // Read grid dimensions
      std::string label;
      in >> label;
      if (label != "mesh" && label != "ngrid") {
         std::string msg =  "\n";
         msg += "Error reading field file:\n";
         msg += "Expected mesh or ngrid, but found [";
         msg += label;
         msg += "]";
         UTIL_THROW(msg.c_str());
      }
      IntVec<D> nGrid;
      in >> nGrid;

      // Initialize mesh and fft
      if (mesh_.size() == 0) {
         mesh_.setDimensions(nGrid);
         fft_.setup(mesh_.dimensions());
      }

      // Allocate waveList
      if (!waveList_.isAllocated()) {
         waveList_.allocate(mesh_, unitCell_);
      }

      // If groupName is present, construct group and basis
      if (groupName_ != "") {
         readGroup(groupName_, group_);
         hasGroup_ = true;
         basis_.makeBasis(mesh_, unitCell_, group_);
      }

      isInitialized_ = true;
   }

   /*
   * Make basis if needed.
   */
   template <int D, class FFT, class WLT, class FIT>
   void DomainTmpl<D,FFT,WLT,FIT>::makeBasis()
   {
      UTIL_CHECK(mesh_.size() > 0);
      UTIL_CHECK(unitCell_.lattice() != UnitCell<D>::Null);
      UTIL_CHECK(unitCell_.isInitialized());
      UTIL_CHECK(hasGroup_);

      // Check basis, construct if not initialized
      if (!basis_.isInitialized()) {
         basis_.makeBasis(mesh_, unitCell_, group_);
      }
      UTIL_CHECK(basis_.isInitialized());
   }

   // Crystallographic Data Output

   /*
   * Write description of symmetry-adapted stars and basis to file.
   */
   template <int D, class FFT, class WLT, class FIT>
   void
   DomainTmpl<D,FFT,WLT,FIT>::writeStars(std::string const & filename)
   const
   {
      UTIL_CHECK(hasGroup());
      UTIL_CHECK(basis_.isInitialized());
      std::ofstream file;
      fileMaster().openOutputFile(filename, file);
      bool isSymmetric = true;
      int nMonomer = 0;
      fieldIo_.writeFieldHeader(file, nMonomer, unitCell_, isSymmetric);
      basis_.outputStars(file);
      file.close();
   }

   /*
   * Write a list of waves and associated stars to file.
   */
   template <int D, class FFT, class WLT, class FIT>
   void
   DomainTmpl<D,FFT,WLT,FIT>::writeWaves(std::string const & filename)
   const
   {
      UTIL_CHECK(hasGroup());
      UTIL_CHECK(basis_.isInitialized());
      std::ofstream file;
      fileMaster().openOutputFile(filename, file);
      bool isSymmetric = true;
      int nMonomer = 0;
      fieldIo_.writeFieldHeader(file, nMonomer, unitCell_, isSymmetric);
      basis_.outputWaves(file);
      file.close();
   }

   /*
   * Write all elements of the space group to a file.
   */
   template <int D, class FFT, class WLT, class FIT>
   void
   DomainTmpl<D,FFT,WLT,FIT>::writeGroup(std::string const & filename)
   const
   {
      UTIL_CHECK(hasGroup());
      std::ofstream file;
      fileMaster().openOutputFile(filename, file);
      file << group_;
      file.close();
   }

} // namespace Prdc
} // namespace Pscf
#endif

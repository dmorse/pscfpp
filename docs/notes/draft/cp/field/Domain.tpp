#ifndef PRDC_CL_DOMAIN_TPP
#define PRDC_CL_DOMAIN_TPP

/*
* PSCF - Polymer Self-Consistent Field
*
* Copyright 2015 - 2025, The Regents of the University of Minnesota
* Distributed under the terms of the GNU General Public License.
*/

#include "Domain.h"
#include <prdc/fieldIo/fieldHeader.h>
#include <util/signal/Signal.h>
#include <util/misc/FileMaster.h>

namespace Pscf {
namespace Cp {

   using namespace Util;
   using namespace Prdc;

   /*
   * Constructor.
   */
   template <int D, class FFT, class WLT, class FIT>
   Domain<D,FFT,WLT,FIT>::Domain()
    : mesh_(),
      unitCell_(),
      lattice_(UnitCell<D>::Null),
      fftPtr_(nullptr),
      waveListPtr_(nullptr),
      fieldIoPtr_(nullptr),
      signalPtr_(nullptr),
      fileMasterPtr_(nullptr),
      isInitialized_(false)
   {
      setClassName("Domain");

      // Construct associated objects
      fftPtr_ = new FFT();
      bool isRealField = false;
      waveListPtr_ = new WLT(isRealField);
      fieldIoPtr_ = new FIT();
      signalPtr_ = new Signal<void>();

      // Note: Passing the WLT (i.e. WaveList<D>) constructor
      // an argument isRealField = false configures it to use
      // a full k-space grid appropriate for a complex field,
      // rather than the half-space grid used for the DFT of
      // a real-valued field. 

      // Create associations between objects
      fieldIo().associate(mesh_, fft(), lattice_);
      unitCell_.setSignal(*signalPtr_);
   }

   /*
   * Destructor.
   */
   template <int D, class FFT, class WLT, class FIT>
   Domain<D,FFT,WLT,FIT>::~Domain()
   {
      delete fftPtr_;
      delete waveListPtr_;
      delete fieldIoPtr_;
      delete signalPtr_;
   }

   /*
   * Create association with a FileMaster.
   */
   template <int D, class FFT, class WLT, class FIT>
   void Domain<D,FFT,WLT,FIT>::setFileMaster(FileMaster& fileMaster)
   {
      fileMasterPtr_ = &fileMaster;
      fieldIo().setFileMaster(fileMaster);
   }

   /*
   * Read parameters and initialize.
   */
   template <int D, class FFT, class WLT, class FIT>
   void Domain<D,FFT,WLT,FIT>::readParameters(std::istream& in)
   {
      // Preconditions
      UTIL_CHECK(!isInitialized_);
      UTIL_CHECK(fileMasterPtr_);

      // Read computational mesh dimensions (required)
      read(in, "mesh", mesh_);
      UTIL_CHECK(mesh_.size() > 0);
      fft().setup(mesh_.dimensions());

      // Read lattice system enumeration value (required)
      read(in, "lattice", lattice_);
      unitCell_.set(lattice_);
      UTIL_CHECK(unitCell_.lattice() != UnitCell<D>::Null);
      UTIL_CHECK(unitCell_.nParameter() > 0);

      // Allocate memory for WaveList
      waveList().allocate(mesh_, unitCell_);

      // Optionally read groupName, and discard information if found
      std::string groupName;
      readOptional(in, "groupName", groupName);

      isInitialized_ = true;
   }

   /*
   * Read header of r-grid field to initialize the Domain.
   *
   * Alternative to parameter file, used only for unit testing.
   */
   template <int D, class FFT, class WLT, class FIT>
   void
   Domain<D,FFT,WLT,FIT>::readFieldHeader(std::istream& in,
                                          int& nMonomer)
   {
      // Preconditions - confirm that nothing is initialized
      UTIL_CHECK(!isInitialized_);
      UTIL_CHECK(lattice_ == UnitCell<D>::Null);
      UTIL_CHECK(!unitCell_.isInitialized());

      // Read common section of standard field header
      int ver1, ver2;
      std::string groupName;
      Pscf::Prdc::readFieldHeader(in, ver1, ver2,
                                  unitCell_, groupName, nMonomer);

      // Set lattice_ (lattice system identifier)
      lattice_ = unitCell_.lattice();
      UTIL_CHECK(lattice_ != UnitCell<D>::Null);
      UTIL_CHECK(unitCell_.isInitialized());

      // Read mesh dimensions
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
         fft().setup(mesh_.dimensions());
      }

      // Allocate waveList
      if (!waveList().isAllocated()) {
         waveList().allocate(mesh_, unitCell_);
      }

      isInitialized_ = true;
   }

} // namespace Cp
} // namespace Pscf
#endif

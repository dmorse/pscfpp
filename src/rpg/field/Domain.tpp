#ifndef RPG_DOMAIN_TPP
#define RPG_DOMAIN_TPP

/*
* PSCF - Polymer Self-Consistent Field Theory
*
* Copyright 2015 - 2025, The Regents of the University of Minnesota
* Distributed under the terms of the GNU General Public License.
*/

#include "Domain.h"
#include <prdc/crystal/fieldHeader.h>

namespace Pscf {
namespace Rpg {

   using namespace Util;
   using namespace Pscf::Prdc;
   using namespace Pscf::Prdc::Cuda;

   /*
   * Constructor.
   */
   template <int D>
   Domain<D>::Domain()
    : unitCell_(),
      mesh_(),
      group_(),
      basis_(),
      fft_(),
      waveList_(),
      fieldIo_(),
      lattice_(UnitCell<D>::Null),
      groupName_(""),
      hasGroup_(false),
      hasFileMaster_(false),
      isInitialized_(false)
   {
      setClassName("Domain");
      fieldIo_.associate(mesh_, fft_, lattice_,
                         hasGroup_, groupName_, group_, basis_);
   }

   /*
   * Destructor.
   */
   template <int D>
   Domain<D>::~Domain()
   {}

   template <int D>
   void Domain<D>::setFileMaster(FileMaster& fileMaster)
   {
      fieldIo_.setFileMaster(fileMaster);
      hasFileMaster_ = true;
   }

   /*
   * Read parameters and initialize.
   */
   template <int D>
   void Domain<D>::readParameters(std::istream& in)
   {
      // Preconditions
      UTIL_CHECK(!isInitialized_);
      UTIL_CHECK(hasFileMaster_);

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

      // If groupName_ exists, construct group_ (space group)
      if (hasGroupName) {
         // Read group symmetry operations from file
         // An Exception is thrown if groupName_ string is not recognized
         readGroup(groupName_, group_);
         hasGroup_ = true;
      }

      isInitialized_ = true;
   }

   /*
   * Read header of r-grid field to initialize the Domain.
   *
   * Alternative to parameter file, used only for unit testing.
   */
   template <int D>
   void Domain<D>::readRGridFieldHeader(std::istream& in, int& nMonomer)
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
   * Set the unit cell by copying a UnitCell<D>, make basis if needed.
   */
   template <int D>
   void Domain<D>::setUnitCell(UnitCell<D> const & unitCell)
   {
      if (lattice_ == UnitCell<D>::Null) {
         lattice_ = unitCell.lattice();
      } else {
         UTIL_CHECK(lattice_ == unitCell.lattice());
      }
      unitCell_ = unitCell;

      UTIL_CHECK(waveList_.isAllocated());
      waveList_.clearUnitCellData();

      if (hasGroup_ && !basis_.isInitialized()) {
         makeBasis();
      }
   }

   /*
   * Set the unit cell, make basis if needed.
   */
   template <int D>
   void Domain<D>::setUnitCell(typename UnitCell<D>::LatticeSystem lattice,
                               FSArray<double, 6> const & parameters)
   {
      if (lattice_ == UnitCell<D>::Null) {
         lattice_ = lattice;
      } else {
         UTIL_CHECK(lattice_ == lattice);
      }
      unitCell_.set(lattice, parameters);

      UTIL_CHECK(waveList_.isAllocated());
      waveList_.clearUnitCellData();

      if (hasGroup_ && !basis_.isInitialized()) {
         makeBasis();
      }
   }

   /*
   * Set unit cell parameters, make basis if needed.
   */
   template <int D>
   void Domain<D>::setUnitCell(FSArray<double, 6> const & parameters)
   {
      UTIL_CHECK(unitCell_.lattice() != UnitCell<D>::Null);
      UTIL_CHECK(unitCell_.nParameter() == parameters.size());
      unitCell_.setParameters(parameters);

      UTIL_CHECK(waveList_.isAllocated());
      waveList_.clearUnitCellData();

      if (hasGroup_ && !basis_.isInitialized()) {
         makeBasis();
      }
   }

   /*
   * Make basis if needed.
   */
   template <int D>
   void Domain<D>::makeBasis()
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

} // namespace Rpg
} // namespace Pscf
#endif

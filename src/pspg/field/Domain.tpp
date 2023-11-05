#ifndef PSPG_DOMAIN_TPP
#define PSPG_DOMAIN_TPP

/*
* PSCF - Polymer Self-Consistent Field Theory
*
* Copyright 2016 - 2022, The Regents of the University of Minnesota
* Distributed under the terms of the GNU General Public License.
*/

#include "Domain.h"

namespace Pscf {
namespace Pspg
{

   using namespace Util;

   /*
   * Constructor.
   */
   template <int D>
   Domain<D>::Domain()
    : unitCell_(),
      mesh_(),
      basis_(),
      fft_(),
      fieldIo_(),
      lattice_(UnitCell<D>::Null),
      groupName_(""),
      hasGroup_(false),
      hasFileMaster_(false),
      isInitialized_(false)
   {  setClassName("Domain"); }

   /*
   * Destructor.
   */
   template <int D>
   Domain<D>::~Domain()
   {}

   template <int D>
   void Domain<D>::setFileMaster(FileMaster& fileMaster)
   {
      fieldIo_.associate(mesh_, fft_, 
                         lattice_, hasGroup_, groupName_, group_, basis_, 
                         fileMaster);
      hasFileMaster_ = true;
   }

   /*
   * Read parameters and initialize.
   */
   template <int D>
   void Domain<D>::readParameters(std::istream& in)
   {
      UTIL_CHECK(hasFileMaster_);

      read(in, "mesh", mesh_);
      UTIL_CHECK(mesh().size() > 0);
      fft_.setup(mesh_.dimensions());

      // Read lattice system identifier (enumeration)
      read(in, "lattice", lattice_);
      unitCell_.set(lattice_);
      UTIL_CHECK(unitCell().lattice() != UnitCell<D>::Null);
      UTIL_CHECK(unitCell().nParameter() > 0);

      // Allocate memory for WaveList
      waveList().allocate(mesh(), unitCell());

      // Optionally read group name 
      hasGroup_ = false;
      bool hasGroupName = false;
      hasGroupName = readOptional(in, "groupName", groupName_).isActive();

      // If group name is found, construct the group
      if (hasGroupName) {
         // Read group symmetry operations from file
         readGroup(groupName_, group_);
         hasGroup_ = true;
      }

      isInitialized_ = true;
   }
 
   template <int D> 
   void Domain<D>::readRGridFieldHeader(std::istream& in, int& nMonomer)
   {
      // Read common section of standard field header
      int ver1, ver2;
      Pscf::Prdc::readFieldHeader(in, ver1, ver2, 
                                  unitCell_, groupName_, nMonomer);

      lattice_ = unitCell_.lattice();
 
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

      // Initialize mesh, fft and basis
      if (mesh_.size() == 0) {
         // Initialize mesh, fft
         mesh_.setDimensions(nGrid);
         fft_.setup(mesh_.dimensions());
      }

      // Initialize group and basis
      if (groupName_ != "") {
         readGroup(groupName_, group_);
         hasGroup_ = true;
         basis_.makeBasis(mesh_, unitCell_, group_);
      }
      
      isInitialized_ = true;
   }

   /*
   * Set the unit cell, make basis if needed.
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
      if (hasGroup_) {
         if (!basis_.isInitialized()) {
            makeBasis();
         }
      }
      waveList_.computeKSq(unitCell_);
      waveList_.computedKSq(unitCell_);
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
      if (!basis_.isInitialized()) {
         makeBasis();
      }
      waveList_.computeKSq(unitCell_);
      waveList_.computedKSq(unitCell_);
   }

   /*
   * Set parameters of the associated unit cell, make basis if needed.
   */
   template <int D>
   void Domain<D>::setUnitCell(FSArray<double, 6> const & parameters)
   {
      UTIL_CHECK(unitCell_.lattice() != UnitCell<D>::Null);
      UTIL_CHECK(unitCell_.nParameter() == parameters.size());
      unitCell_.setParameters(parameters);
      if (!basis_.isInitialized()) {
         makeBasis();
      }
      waveList_.computeKSq(unitCell_);
      waveList_.computedKSq(unitCell_);
   }

   template <int D>
   void Domain<D>::makeBasis() 
   {
      UTIL_CHECK(mesh_.size() > 0);
      UTIL_CHECK(unitCell_.lattice() != UnitCell<D>::Null);
      UTIL_CHECK(unitCell_.nParameter() > 0);
      UTIL_CHECK(unitCell_.isInitialized());
      UTIL_CHECK(hasGroup_);

      // Check basis, construct if not initialized
      if (!basis().isInitialized()) {
         basis_.makeBasis(mesh_, unitCell_, group_);
      }
      UTIL_CHECK(basis().isInitialized());

      // Compute minimum images in WaveList
      UTIL_CHECK(waveList().isAllocated());
      if (!waveList().hasMinimumImages()) {
         waveList().computeMinimumImages(mesh(), unitCell());
      }

   }

} // namespace Pspg
} // namespace Pscf
#endif

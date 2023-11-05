#ifndef PSPC_DOMAIN_TPP
#define PSPC_DOMAIN_TPP

/*
* PSCF - Polymer Self-Consistent Field Theory
*
* Copyright 2016 - 2022, The Regents of the University of Minnesota
* Distributed under the terms of the GNU General Public License.
*/

#include "Domain.h"
#include <prdc/crystal/fieldHeader.h>

namespace Pscf {
namespace Pspc
{

   using namespace Util;
   using namespace Pscf::Prdc;
   using namespace Pscf::Prdc::Cpu;

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
      fieldIo_.associate(mesh_, fft_, lattice_, 
                         hasGroup_, groupName_, group_, basis_, 
                         fileMaster);
      hasFileMaster_ = true;
   }

   /*
   * Read parameters and initialize.
   */
   template <int D>
   void Domain<D>::readParameters(std::istream& in)
   {
      // Preconditins
      UTIL_CHECK(!isInitialized_);
      UTIL_CHECK(hasFileMaster_);
      bool hasUnitCell = false; 

      read(in, "mesh", mesh_);
      fft_.setup(mesh_.dimensions());

      read(in, "lattice", lattice_);
      unitCell_.set(lattice_);

      // Optionally read group name 
      bool hasGroupName;
      hasGroupName = readOptional(in, "groupName", groupName_).isActive();

      hasGroup_ = false;
      if (hasGroupName) {

         // Read the space group from file.
         readGroup(groupName_, group_);
         hasGroup_ = true;

         // If unit cell parameters are known, construct a basis
         if (hasUnitCell) { 
            basis().makeBasis(mesh(), unitCell(), group_);
         }
      }

      isInitialized_ = true;
   }

   /*
   * Read header of r-grid field in order to initialize the Domain.
   *
   * Used only for unit testing.
   */ 
   template <int D> 
   void Domain<D>::readRGridFieldHeader(std::istream& in, int& nMonomer)
   {
      // Preconditions - confirm that nothing is initialized
      UTIL_CHECK(!isInitialized_)
      UTIL_CHECK(lattice_ == UnitCell<D>::Null);
      UTIL_CHECK(!unitCell_.isInitialized());
      UTIL_CHECK(!hasGroup_);
      UTIL_CHECK(groupName_ == "");

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

      if (mesh_.size() == 0) {
         // Initialize mesh, fft 
         mesh_.setDimensions(nGrid);
         fft_.setup(mesh_.dimensions());
      }

      if (groupName_ != "") {
         readGroup(groupName_, group_);
         hasGroup_ = true;
         basis_.makeBasis(mesh_, unitCell_, group_);
      }
      
      isInitialized_ = true;
   }

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
   }

   /*
   * Set parameters of the associated unit cell.
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
      if (hasGroup_) {
         if (!basis_.isInitialized()) {
            makeBasis();
         }
      }
   }

   /*
   * Set parameters of the associated unit cell.
   */
   template <int D>
   void Domain<D>::setUnitCell(FSArray<double, 6> const & parameters)
   {
      UTIL_CHECK(unitCell_.lattice() != UnitCell<D>::Null);
      UTIL_CHECK(unitCell_.nParameter() == parameters.size());
      unitCell_.setParameters(parameters);
      if (hasGroup_) {
         if (!basis_.isInitialized()) {
            makeBasis();
         }
      }
   }

   template <int D>
   void Domain<D>::makeBasis() 
   {
      UTIL_CHECK(mesh_.size() > 0);
      UTIL_CHECK(unitCell_.lattice() != UnitCell<D>::Null);
      UTIL_CHECK(hasGroup_);

      #if 0
      // Check group, read from file if necessary
      if (!hasGroup_ && groupName_ != "") {
         readGroup(groupName_, group_);
         hasGroup_ = true;
      }
      #endif

      // Check basis, construct if not initialized
      if (!basis().isInitialized()) {
         basis_.makeBasis(mesh_, unitCell_, group_);
      }
      UTIL_CHECK(basis().isInitialized());
   }

} // namespace Pspc
} // namespace Pscf
#endif

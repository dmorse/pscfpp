#ifndef PSPC_DOMAIN_TPP
#define PSPC_DOMAIN_TPP

/*
* PSCF - Polymer Self-Consistent Field Theory
*
* Copyright 2016 - 2022, The Regents of the University of Minnesota
* Distributed under the terms of the GNU General Public License.
*/

#include "Domain.h"

namespace Pscf {
namespace Pspc
{

   using namespace Util;
   using namespace Pscf::Prdc;

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
      groupName_(),
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
                         lattice_, groupName_, group_, basis_, 
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
      // Uncomment for backwards compatibility with old format (< v1.0)
      #if 0
      readOptional(in, "unitCell", unitCell_);
      if (unitCell_.lattice() != UnitCell<D>::Null) { 
         lattice_ = unitCell_.lattice();
         hasUnitCell = true; 
      }
      #endif

      read(in, "mesh", mesh_);
      fft_.setup(mesh_.dimensions());

      // If no unit cell was read, read lattice system
      if (lattice_ == UnitCell<D>::Null) { 
         read(in, "lattice", lattice_);
         unitCell_.set(lattice_);
      }

      // Read group name and initialized space group
      read(in, "groupName", groupName_);
      readGroup(groupName_, group_);

      // Initialize unit cell if unit cell parameters are known
      if (hasUnitCell) { 
         basis().makeBasis(mesh(), unitCell(), group_);
      }
      isInitialized_ = true;
   }

   /*
   * Read header of r-grid field in order to initialize the Domain.
   *
   * Useful for unit testing.
   */ 
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

      // Initialize mesh, fft 
      mesh_.setDimensions(nGrid);
      fft_.setup(mesh_.dimensions());

      // Initialize group and basis
      readGroup(groupName_, group_);
      basis_.makeBasis(mesh_, unitCell_, group_);
      
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
      if (!basis_.isInitialized()) {
         makeBasis();
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
      if (!basis_.isInitialized()) {
         makeBasis();
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
      if (!basis_.isInitialized()) {
         makeBasis();
      }
   }

   template <int D>
   void Domain<D>::makeBasis() 
   {
      UTIL_CHECK(mesh_.size() > 0);
      UTIL_CHECK(unitCell_.lattice() != UnitCell<D>::Null);

      // Check group, read from file if necessary
      if (group_.size() == 1) {
         if (groupName_ != "I") {
            UTIL_CHECK(groupName_ != "");
            readGroup(groupName_, group_);
         }
      }

      // Check basis, construct if not initialized
      if (!basis().isInitialized()) {
         basis_.makeBasis(mesh_, unitCell_, group_);
      }
      UTIL_CHECK(basis().isInitialized());
   }

} // namespace Pspc
} // namespace Pscf
#endif

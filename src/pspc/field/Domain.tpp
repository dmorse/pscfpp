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
      UTIL_CHECK(hasFileMaster_);

      readOptional(in, "unitCell", unitCell_);
      bool hasUnitCell = false; 
      if (unitCell_.lattice() != UnitCell<D>::Null) { 
         lattice_ = unitCell_.lattice();
         hasUnitCell = true; 
      }

      read(in, "mesh", mesh_);
      fft_.setup(mesh_.dimensions());

      read(in, "groupName", groupName_);
      readGroup(groupName_, group_);

      if (hasUnitCell) { 
         basis().makeBasis(mesh(), unitCell(), group_);
      }
      isInitialized_ = true;
   }
   
 
   template <int D> 
   void Domain<D>::readFieldHeader(std::istream& in, int& nMonomer)
   {
      // Read common section of standard field header
      int ver1, ver2;
      Pscf::readFieldHeader(in, ver1, ver2, 
                           unitCell_, groupName_, nMonomer);

      lattice_ = unitCell_.lattice();
 
      // Read grid dimensions
      std::string label;
      in >> label;
      UTIL_CHECK(label == "ngrid");
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
      if (lattice_ != UnitCell<D>::Null) {
         lattice_ = unitCell.lattice();
      } else {
         UTIL_CHECK(lattice_ = unitCell.lattice());
      }
      unitCell_ = unitCell;
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

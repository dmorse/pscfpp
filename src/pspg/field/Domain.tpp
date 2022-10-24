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

      // Optionally read unit cell
      read(in, "unitCell", unitCell_);
      bool hasUnitCell = false;
      if (unitCell_.lattice() != UnitCell<D>::Null) {
         lattice_ = unitCell_.lattice();
         hasUnitCell = true;
      }

      read(in, "mesh", mesh_);
      fft_.setup(mesh_.dimensions());

      // If no unit cell was read, read lattice system 
      if (!hasUnitCell) {
         read(in, "lattice", lattice_);
         unitCell_.set(lattice_);
      }

      // Read group name and initialize space group
      read(in, "groupName", groupName_);
      readGroup(groupName_, group_);

      // Make symmetry-adapted basis
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
 
      // Read grid dimensions
      std::string label;
      in >> label;
      UTIL_CHECK(label == "ngrid");
      IntVec<D> nGrid;
      in >> nGrid;

      // Initialize mesh, fft and basis
      mesh_.setDimensions(nGrid);
      fft_.setup(mesh_.dimensions());
      basis_.makeBasis(mesh_, unitCell_, groupName_);
      
      isInitialized_ = true;
   }

} // namespace Pspg
} // namespace Pscf
#endif

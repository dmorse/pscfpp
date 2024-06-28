#ifndef RPC_FILM_ITERATOR_BASE_TPP
#define RPC_FILM_ITERATOR_BASE_TPP

/*
* PSCF - Polymer Self-Consistent Field Theory
*
* Copyright 2016 - 2021, The Regents of the University of Minnesota
* Distributed under the terms of the GNU General Public License.
*/

#include "FilmIteratorBase.h"
#include "rpc/System.h"
#include "rpc/field/FieldIo.h"

#include "prdc/cpu/RField.h"
#include "prdc/crystal/UnitCell.h"
#include "prdc/crystal/groupFile.h"
#include "prdc/crystal/SpaceGroup.h"

#include "pscf/math/RealVec.h"

#include "util/containers/FArray.h"
#include "util/format/Dbl.h"

#include <cmath>
#include <iostream>

namespace Pscf {
namespace Rpc
{

   using namespace Util;
   using namespace Pscf::Prdc;
   using namespace Pscf::Prdc::Cpu;

   /*
   * Constructor.
   */
   template <int D, typename IteratorType>
   FilmIteratorBase<D, IteratorType>::FilmIteratorBase(System<D>& system)
    : Iterator<D>(system),
      iterator_(system),
      parameters_(),
      normalVecId_(-1),
      interfaceThickness_(-1.0),
      excludedThickness_(-1.0),
      chiBottom_(),
      chiTop_(),
      chiBottomCurrent_(),
      chiTopCurrent_(),
      ungenerated_(true)
   {  
      isSymmetric_ = true;
      setClassName(iterator_.className().append("FilmBase").c_str());
      system.mask().setFieldIo(system.fieldIo());
      system.h().setFieldIo(system.fieldIo());
   }

   /*
   * Destructor.
   */
   template <int D, typename IteratorType>
   FilmIteratorBase<D, IteratorType>::~FilmIteratorBase()
   {}

   /*
   * Read iterator and wall data from parameter file
   */
   template <int D, typename IteratorType>
   void FilmIteratorBase<D, IteratorType>::readParameters(std::istream& in)
   {
      // Make iterator_ a child paramComponent of this object, so that it 
      // will be read/written correctly to/from param file with correct 
      // indentation
      setParent(iterator_,false);
      addComponent(iterator_,false);
      
      // Read the Iterator parameters
      iterator_.readParameters(in);

      // Read required data defining the walls
      read(in, "normalVecId", normalVecId_);
      read(in, "interfaceThickness", interfaceThickness_);
      read(in, "excludedThickness", excludedThickness_);

      // Make sure inputs are valid
      if (normalVecId_ > D || normalVecId_ < 0) {
         UTIL_THROW("bad value for normalVecId, must be in [0,D)");
      }
      if (interfaceThickness_ > excludedThickness_) {
         UTIL_THROW("excludedThickness must be larger than interfaceThickness");
      }
      if ((excludedThickness_ <= 0) || (interfaceThickness_ <= 0)) {
         UTIL_THROW("excludedThickness and interfaceThickness must be >0");
      }

      // Allocate chiBottom_ and chiTop_ and set to zero before 
      // reading them in
      int nm = system().mixture().nMonomer();
      chiBottom_.allocate(nm);
      chiTop_.allocate(nm);
      for (int i = 0; i < nm; i++) {
         chiBottom_[i] = 0.0;
         chiTop_[i] = 0.0;
      }
      readDArray(in, "chiBottom", chiBottom_, nm);
      readDArray(in, "chiTop", chiTop_, nm);

      // If lattice parameters are flexible, determine which parameters
      // are allowed to vary, store them in this object, and pass them
      // into iterator_. The flexibleParams_ member of the iterator_
      // should always be matched to that of this class.
      setFlexibleParams();
   }

   /*
   * Allocate required memory, perform necessary checks to ensure user 
   * input is compatible with a film constraint, and create the mask / 
   * external fields that will be used to represent the walls during 
   * iteration.
   */
   template <int D, typename IteratorType>
   void FilmIteratorBase<D, IteratorType>::setup()
   {
      UTIL_CHECK(system().basis().isInitialized());
      UTIL_CHECK(system().unitCell().isInitialized());

      // Allocate the mask and external field containers if needed
      if (!system().mask().isAllocated()) {
         system().mask().allocate(system().basis().nBasis(), 
                                 system().mesh().dimensions());
      }
      if (!isAthermal()) {
         if (!system().h().isAllocatedRGrid()) {
            system().h().allocateRGrid(system().mesh().dimensions());
         }
         if (!system().h().isAllocatedBasis()) {
            system().h().allocateBasis(system().basis().nBasis());
         }
      }

      // Ensure that space group symmetry is compatible with the wall
      checkSpaceGroup();

      if (ungenerated_) {
         // Generate the field representation of the walls (and external
         // fields if needed) and provide them to iterator_ to use during 
         // iteration. 
         generateWallFields();
         ungenerated_ = false;
      } else {
         // Update the concentration field of the walls based on the
         // current state of Domain<D> and provide it to iterator_, as
         // well as the external fields for each species if needed
         updateWallFields(); 
      }
   }

   /*
   * Iterate to a solution
   */
   template <int D, typename IteratorType>
   int FilmIteratorBase<D, IteratorType>::solve(bool isContinuation)
   {
      // Setup the FilmIteratorBase object, set external fields and mask
      setup();

      // solve
      return iterator_.solve(isContinuation);
   }

   /*
   * Generate the concentration field for the walls and pass it into 
   * iterator. Adjust phi values of species in the system to accommodate 
   * the volume occupied by the wall. Finally, generate the external 
   * fields that are created by the wall, if needed (if the wall is not 
   * athermal), and provide them to the iterator.
   */
   template <int D, typename IteratorType>
   void FilmIteratorBase<D, IteratorType>::generateWallFields() 
   {
      UTIL_CHECK(interfaceThickness() > 0);
      UTIL_CHECK(excludedThickness() > interfaceThickness());
      UTIL_CHECK(system().mask().isAllocated());

      if (ungenerated_) ungenerated_ = false;

      // Ensure that unit cell is compatible with wall
      checkLatticeVectors();

      // Get the length L of the lattice basis vector normal to the walls
      RealVec<D> a;
      a = system().domain().unitCell().rBasis(normalVecId_);
      double norm_sqd(0.0); // norm squared
      for (int i = 0; i < D; i++) {
         norm_sqd += a[i]*a[i];
      }
      double L(sqrt(norm_sqd));

      // Create a 3 element vector 'dim' that contains the grid dimensions.
      // If system is 2D (1D), then the z (y & z) dimensions are set to 1.
      IntVec<3> dim;
      for (int ind = 0; ind < 3; ind++) {
         if (ind < D) {
            dim[ind] = system().domain().mesh().dimensions()[ind];
         } else {
            dim[ind] = 1;
         }
      }

      // Generate an r-grid representation of the walls
      RField<D> rGrid;
      rGrid.allocate(system().domain().mesh().dimensions());
      int x, y, z;
      int counter = 0;
      FArray<int,3> coords;
      double d, rho_w;

      for (x = 0; x < dim[0]; x++) {
         coords[0] = x;
         for (y = 0; y < dim[1]; y++) {
            coords[1] = y;
            for (z = 0; z < dim[2]; z++) {
               coords[2] = z;

               // Get the distance 'd' traveled along the lattice basis 
               // vector that is normal to the walls
               d = coords[normalVecId_] * L / dim[normalVecId_];

               // Calculate wall volume fraction (rho_w) at gridpoint (x,y,z)
               rho_w = 0.5 * (1 + tanh(4 * (((.5 * (excludedThickness_-L)) + 
                                  fabs(d - (L/2))) / interfaceThickness_)));
               rGrid[counter++] = 1-rho_w;
            }
         }
      }

      // Store this mask in System
      system().mask().setRGrid(rGrid,true);

      // Store lattice parameters associated with this maskBasis
      parameters_ = system().domain().unitCell().parameters();

      // Generate external fields if needed
      generateExternalFields();
   }

   /*
   * Update the mask in the iterator (the region in which the polymers are
   * confined) if the lattice parameters have been updated since mask was 
   * last generated. Also update external fields.
   * 
   * If the lattice parameters have not changed but the wall/polymer chi 
   * parameters have been updated since external fields were last generated,
   * update the external fields (but not the mask).
   */
   template <int D, typename IteratorType>
   void FilmIteratorBase<D, IteratorType>::updateWallFields() 
   {
      // Check if system lattice parameters are different than parameters_
      FSArray<double, 6> sysParams = 
         system().domain().unitCell().parameters();
      UTIL_CHECK(sysParams.size() == parameters_.size());
      bool identical = true; // are the two arrays identical?
      for (int i = 0; i < parameters_.size(); i++) {
         if (fabs(sysParams[i] - parameters_[i]) > 1e-10) {
            identical = false;
            break;
         }
      }

      // Regenerate wall fields if the lattice parameters have changed
      if (!identical) {
         generateWallFields();
      } else {
         // If we did not regenerate wall fields, check if we need to
         // regenerate external field (if chi array has changed)
         bool newExternalFields = false; // Do we need new external fields?
         for (int i = 0; i < chiBottom_.capacity(); i++) {
            if ((chiBottom_[i] != chiBottomCurrent_[i]) || 
                (chiTop_[i] != chiTopCurrent_[i])) {
               newExternalFields = true;
               break;
            }
         }
         if (newExternalFields) {
            generateExternalFields();
         }
      }
   }

   /*
   * Generate the external fields that are imposed as a result of chemical
   * interactions between the walls and the monomer species.
   */
   template <int D, typename IteratorType>
   void FilmIteratorBase<D, IteratorType>::generateExternalFields() 
   {
      // Set chiBottomCurrent_ and chiTopCurrent_ equal to the current chi 
      // arrays associated with the external fields that are about to be set
      chiBottomCurrent_ = chiBottom_;
      chiTopCurrent_ = chiTop_;

      // If walls are athermal then there is no external field needed.
      // If an external field already exists in the System, we need to
      // overwrite it with a field of all zeros, otherwise do nothing
      if ((isAthermal()) && (!system().h().hasData())) { 
         return; 
      }

      // If this point is reached, external field must be generated
      UTIL_CHECK(system().h().isAllocatedRGrid());
      UTIL_CHECK(system().h().isAllocatedBasis());
      int nm = system().mixture().nMonomer();

      // Get length L of the lattice basis vector normal to the walls
      RealVec<D> a;
      a = system().domain().unitCell().rBasis(normalVecId_);
      double norm_sqd(0.0); // norm squared
      for (int i = 0; i < D; i++) {
         norm_sqd += a[i]*a[i];
      }
      double L(sqrt(norm_sqd));

      // Create a 3 element vector 'dim' that contains the grid 
      // dimensions. If system is 2D (1D), then the z (y and z) 
      // dimensions are set to 1.
      IntVec<3> dim;
      for (int ind = 0; ind < 3; ind++) {
         if (ind < D) {
            dim[ind] = system().domain().mesh().dimensions()[ind];
         } else {
            dim[ind] = 1;
         }
      }

      // Generate an r-grid representation of the external fields
      DArray< RField<D> > hRGrid;
      hRGrid.allocate(nm);
      for (int i = 0; i < nm; i++) {
         hRGrid[i].allocate(system().domain().mesh().dimensions());
      }

      int i, x, y, z;
      int counter = 0;
      FArray<int,3> coords;
      double d, rho_w;

      for (i = 0; i < nm; i++) {
         for (x = 0; x < dim[0]; x++) {
            coords[0] = x;
            for (y = 0; y < dim[1]; y++) {
               coords[1] = y;
               for (z = 0; z < dim[2]; z++) {
                  coords[2] = z;

                  // Get the distance 'd' traveled along the lattice 
                  // basis vector that is orthogonal to the walls
                  d = coords[normalVecId_] * L / dim[normalVecId_];

                  // Calculate wall volume fraction (rho_w) at gridpoint 
                  // (x,y,z)
                  rho_w = 0.5 * (1 + tanh(4 * (((.5 * (excludedThickness_-L)) + 
                                     fabs(d - (L/2))) / interfaceThickness_)));
                  if (d < (L/2)) {
                     hRGrid[i][counter++] = rho_w * chiBottom_[i];
                  } else {
                     hRGrid[i][counter++] = rho_w * chiTop_[i];
                  }
               }
            }
         }
         counter = 0;
      } 

      // Pass h into the System
      system().h().setRGrid(hRGrid,true);
   }

   /*
   * Check that user-defined space group is compatible 
   * with the thin film constraint
   */
   template <int D, typename IteratorType>
   void FilmIteratorBase<D, IteratorType>::checkSpaceGroup() const
   {
      // Setup
      std::string groupName = system().groupName();
      SpaceGroup<D> group;
      std::ifstream in;

      // Open and read file containing space group's symmetry operations
      if (groupName == "I") {
         // Create identity group by default
         group.makeCompleteGroup();
      } else {
         bool foundFile = false;
         {
            // Search first in this directory
            in.open(groupName);
            if (in.is_open()) {
               in >> group;
               UTIL_CHECK(group.isValid());
               foundFile = true;
            }
         }
         if (!foundFile) {
            // Search in the data directory containing standard space groups
            std::string fileName = makeGroupFileName(D, groupName);
            in.open(fileName);
            if (in.is_open()) {
               in >> group;
               UTIL_CHECK(group.isValid());
            } else {
               Log::file() << "\nFailed to open group file: " 
                           << fileName << "\n";
               Log::file() << "\n Error: Unknown space group\n";
               UTIL_THROW("Unknown space group");
            }
         } 
      }

      // Make sure all symmetry operations are allowed
      int nv = normalVecId();
      bool symmetric = hasSymmetricWalls();
      std::string msg = "Space group contains forbidden symmetry operations";
      for (int i = 0; i < group.size(); i++) {
         for (int j = 0; j < D; j++) {
            int r = group[i].R(nv,j);
            if (j == nv) {
               if (r != 1) {
                  if ((r != -1) || (!symmetric)) {
                     UTIL_THROW(msg.c_str());
                  }
               }
            } else { // j != nv
               if (r != 0) {
                  UTIL_THROW(msg.c_str());
               }
            }
         }
         if (group[i].t(nv) != 0) {
            UTIL_THROW(msg.c_str());
         }
      }
   }

   /*
   * Check whether or not the two walls are chemically identical using
   * the chi array.
   */
   template <int D, typename IteratorType>
   bool FilmIteratorBase<D, IteratorType>::hasSymmetricWalls() const 
   {
      int nm = system().mixture().nMonomer();

      UTIL_CHECK(nm > 0);
      UTIL_CHECK(chiBottom_.capacity() == nm);
      UTIL_CHECK(chiTop_.capacity() == nm);

      for (int i = 0; i < nm; i++) {
         if (fabs(chiBottom_[i]-chiTop_[i]) > 1e-7) {
            return false;
         }
      }
      return true;
   }

   /*
   * Check whether or not the walls are athermal (only true if all values
   * in chi array are zero)
   */
   template <int D, typename IteratorType>
   bool FilmIteratorBase<D, IteratorType>::isAthermal() const 
   {
      int nm = system().mixture().nMonomer();

      UTIL_CHECK(nm > 0);
      UTIL_CHECK(chiBottom_.capacity() == nm);
      UTIL_CHECK(chiTop_.capacity() == nm);

      for (int i = 0; i < nm; i++) {
         if ((fabs(chiBottom_[i]) >= 1e-7) || (fabs(chiTop_[i]) >= 1e-7)) {
            return false;
         }
      }
      return true;
   }

   /*
   * Add specialized sweep parameter types chi_top and chi_bottom to 
   * the Sweep object.
   */
   template <int D, typename IteratorType>
   void FilmIteratorBase<D, IteratorType>::addParameterTypes(Sweep<D>& sweep)
   {
      sweep.addParameterType("chi_top", 1, *this);
      sweep.addParameterType("chi_bottom", 1, *this);
   }

   /*
   * Set the value of a specialized sweep parameter (chi_top or chi_bottom).
   */
   template <int D, typename IteratorType>
   void FilmIteratorBase<D, IteratorType>::setParameter(std::string name, 
                                            DArray<int> ids, double value)
   {
      if (name == "chi_top") {
         chiTop_[ids[0]] = value;
      } else if (name == "chi_bottom") {
         chiBottom_[ids[0]] = value;
      } else {
         UTIL_THROW(("Parameter name " + name + " not recognized.").c_str());
      }
   }

   /*
   * Get the value of a specialized sweep parameter (chi_top or chi_bottom).
   */
   template <int D, typename IteratorType>
   double FilmIteratorBase<D, IteratorType>::getParameter(std::string name, 
                                                      DArray<int> ids) const
   {
      if (name == "chi_top") {
         return chiTop_[ids[0]];
      } else if (name == "chi_bottom") {
         return chiBottom_[ids[0]];
      } else {
         UTIL_THROW(("Parameter name " + name + " not recognized.").c_str());
         return 0.0;
      }
   }

} // namespace Rpc
} // namespace Pscf
#endif

#ifndef PSPC_FILM_ITERATOR_BASE_TPP
#define PSPC_FILM_ITERATOR_BASE_TPP

/*
* PSCF - Polymer Self-Consistent Field Theory
*
* Copyright 2016 - 2021, The Regents of the University of Minnesota
* Distributed under the terms of the GNU General Public License.
*/

#include "FilmIteratorBase.h"
#include "pscf/crystal/UnitCell.h"
#include "pscf/crystal/groupFile.h"
#include "pscf/math/RealVec.h"
#include "pscf/crystal/SpaceGroup.h"
#include "pspc/System.h"
#include "pspc/field/RField.h"
#include "pspc/field/FieldIo.h"
#include "util/containers/FArray.h"
#include <cmath>
#include <iostream>

namespace Pscf {
namespace Pspc
{

   using namespace Util;

   /*
   * Constructor.
   */
   template <int D, typename IteratorType>
   FilmIteratorBase<D, IteratorType>::FilmIteratorBase(System<D>& system)
    : Iterator<D>(system),
      iterator_(system),
      wallCField_(),
      parameters_(),
      normalVec_(-1),
      t_(-1.0),
      T_(-1.0),
      chi_(),
      writeBasis_(),
      writeRGrid_()
   {  setClassName(iterator_.className().append("FilmBase").c_str()); }

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
      // First, read the Iterator parameters
      iterator_.readParameters(in);

      // Then, read required data defining the walls
      read(in, "normalVec", normalVec_);
      read(in, "interfaceThickness", t_);
      read(in, "wallThickness", T_);

      // Make sure inputs are valid
      if (normalVec_ > D || normalVec_ < 0) {
         UTIL_THROW("bad value for normalVec, must be in [0,D)");
      }
      if (t_ > T_) {
         UTIL_THROW("wallThickness must be larger than interfaceThickness");
      }
      if ((T_ <= 0) || (t_ <= 0)) {
         UTIL_THROW("wallThickness and interfaceThickness must be >0");
      }

      // Allocate chi_ and set it to zero before optionally reading it in
      int nm = system().mixture().nMonomer();
      chi_.allocate(nm,2);
      for (int i = 0; i < nm; i++) {
         for (int j = 0; j < nm; j++) {
            chi_(i,j) = 0;
         }
      }
      readOptionalDMatrix(in, "chi", chi_, nm, 2);

      // Read optional filenames if user wants to write wall basis / rgrid files
      readOptional(in, "writeWallBasis", writeBasis_);
      readOptional(in, "writeWallRGrid", writeRGrid_);
   }

   /*
   * Allocate required memory, perform necessary checks to ensure user input is
   * compatible with a film constraint, and create the mask / external fields
   * that will be used to represent the walls during iteration.
   */
   template <int D, typename IteratorType>
   void FilmIteratorBase<D, IteratorType>::setup()
   {
      // Set up the iterator
      iterator_.setup();

      // Ensure that the space group symmetry is compatible with the wall
      checkSpaceGroup();
      
      // Generate the field representation of the walls (and external fields
      // if needed) and provide them to iterator_ to use during iteration
      generateWallFields();

      // Write wall field to a file if requested in param file
      if (writeBasis_.size() > 0) {
         writeWallField(writeBasis_);
      }
      if (writeRGrid_.size() > 0) {
         writeWallFieldRGrid(writeRGrid_);
      }

      // If lattice parameters are flexible, determine which parameters are 
      // allowed to vary and pass them into iterator_.
      if (isFlexible()) {
         iterator_.setFlexibleParams(flexibleParams());
      }
      
      // Adjust phi values for all species to account for 
      // the volume that is occupied by the wall
      adjustPhiVals();
   }

   /*
   * Iterate to a solution
   */
   template <int D, typename IteratorType>
   int FilmIteratorBase<D, IteratorType>::solve()
   {
      // Update the concentration field of the walls based on the current
      // state of Domain<D> and provide it to iterator_, as well as the 
      // external fields for each monomer species, if needed
      updateWallFields(); 

      // solve
      return iterator_.solve();
   }

   /*
   * Generate the concentration field for the walls and store it in wallCField_
   */
   template <int D, typename IteratorType>
   void FilmIteratorBase<D, IteratorType>::generateWallFields() 
   {
      UTIL_CHECK(interfaceThickness() > 0);
      UTIL_CHECK(wallThickness() > interfaceThickness());

      // Ensure that unit cell is compatible with wall
      checkLatticeVectors();

      int nb = system().domain().basis().nBasis();

      if (!wallCField_.isAllocated()) {
         wallCField_.allocate(nb);
      } else {
         UTIL_CHECK(wallCField_.capacity() == nb);
      }

      // Get the length L of the lattice basis vector orthogonal to the walls
      RealVec<D> a;
      a = system().domain().unitCell().rBasis(normalVec_);
      double norm_sqd; // norm squared
      for (int i = 0; i < D; i++) {
         norm_sqd = a[i]*a[i];
      }
      double L(sqrt(norm_sqd));

      // Create a 3 element vector 'dim' that contains the grid dimensions. If 
      // system is 2D (1D), then the z (y and z) dimensions are set to 1.
      IntVec<3> dim;
      for (int ind = 0; ind < 3; ind++) {
         if (ind < D) {
            dim[ind] = system().domain().mesh().dimensions()[ind];
         } else {
            dim[ind] = 1;
         }
      }

      // Generate an r-grid representation of the walls
      DArray< RField<D> > rGrid;
      rGrid.allocate(1);
      rGrid[0].allocate(system().domain().mesh().dimensions());
      int x, y, z;
      int counter = 0;
      FArray<int,3> coords;
      double d, rho;

      for (x = 0; x < dim[0]; x++) {
         coords[0] = x;
         for (y = 0; y < dim[1]; y++) {
            coords[1] = y;
            for (z = 0; z < dim[2]; z++) {
               coords[2] = z;

               // Get the distance 'd' traveled along the lattice basis vector 
               // that is orthogonal to the walls
               d = coords[normalVec_] * L / dim[normalVec_];

               // Calculate volume fraction of wall (rho) at gridpoint (x,y,z)
               rho = 0.5*(1+tanh(4*(((.5*(T_-L))+fabs(d-(L/2)))/t_)));
               rGrid[0][counter++] = rho;
            }
         }
      }

      // Convert r-grid representation of walls to basis format, stored in wallCField_
      DArray< DArray <double> > cField;
      cField.allocate(1);
      cField[0].allocate(nb);
      system().domain().fieldIo().convertRGridToBasis(rGrid, cField);

      wallCField_ = cField[0];
      iterator_.setMask(wallCField_);

      // Store lattice parameters associated with this wallCField_
      parameters_ = system().domain().unitCell().parameters();
   }

   /*
   * Update the concentration field for the walls in wallCField_ if the lattice
   * parameters have been updated since we last generated wallCField_
   */
   template <int D, typename IteratorType>
   void FilmIteratorBase<D, IteratorType>::updateWallFields() 
   {
      // Check if system lattice parameters are different than parameters_
      FSArray<double, 6> sysParams = system().domain().unitCell().parameters();
      UTIL_CHECK(sysParams.size() == parameters_.size()); // ensure same size
      bool identical = true; // are the two arrays identical?
      for (int i = 0; i < parameters_.size(); i++) {
         if (sysParams[i] != parameters_[i]) {
            identical = false;
            break;
         }
      }

      // Regenerate wall fields if the lattice parameters have changed
      if (!identical) {
         generateWallFields();
      }
   }

   template <int D, typename IteratorType>
   void FilmIteratorBase<D, IteratorType>::adjustPhiVals()
   {
      int np = system().mixture().nPolymer(); 
      int ns = system().mixture().nSolvent(); 
      UTIL_CHECK(np + ns > 0);

      // Get the total volume fraction currently occupied by polymers/solvents
      double frac = 0;
      for (int i = 0; i < np; i++) {
         frac += system().mixture().polymer(i).phi();
      }
      for (int j = 0; j < ns; j++) {
         frac += system().mixture().solvent(j).phi();
      }

      // Adjust phi for each species based on phi_walls and frac, if needed
      double multiplier = (1 - phi()) / frac;
      if (multiplier - 1 > 1e-7) { // if multiplier != 1
         Log::file() << "NOTE: phi values of all species are being scaled by a factor of"
                     << std::endl << multiplier
                     << " due to the volume occupied by the wall." 
                     << std::endl;

         double phi;
         for (int i = 0; i < np; i++) {
            phi = system().mixture().polymer(i).phi();
            system().mixture().polymer(i).setPhi(phi * multiplier);
         }
         for (int j = 0; j < ns; j++) {
            phi = system().mixture().solvent(j).phi();
            system().mixture().solvent(j).setPhi(phi * multiplier);
         }
      }
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
      std::string fileName = makeGroupFileName(D, groupName);
      std::ifstream in;

      // Open and read file containing the space group's symmetry operations
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

      // Make sure all symmetry operations are allowed
      int nv = normalVec();
      bool symmetric = isSymmetric();
      for (int i = 0; i < group.size(); i++) {
         for (int j = 0; j < D; j++) {
            int r = group[i].R(nv,j);
            if (j == nv) {
               if (r != 1) {
                  if ((r != -1) || (!symmetric)) {
                     UTIL_THROW("Space group contains forbidden symmetry operations");
                  }
               }
            } else { // j != nv
               if (r != 0) {
                  UTIL_THROW("Space group contains forbidden symmetry operations");
               }
            }
         }
         if (group[i].t(nv) != 0) {
            UTIL_THROW("Space group contains forbidden symmetry operations");
         }
      }
   }

   /*
   * Check whether or not the two walls are chemically identical using
   * the chi array.
   */
   template <int D, typename IteratorType>
   bool FilmIteratorBase<D, IteratorType>::isSymmetric() const 
   {
      int nm = system().mixture().nMonomer();

      UTIL_CHECK(nm > 0);
      UTIL_CHECK(chi().capacity1() == nm);
      UTIL_CHECK(chi().capacity2() == 2);

      for (int i = 0; i < nm; i++) {
         if (chi(i,0) != chi(i,1)) {
            return false;
         }
      }
      return true;
   }

   /*
   * Check whether or not the walls are athermal (only true if all values in
   * chi array are zero)
   */
   template <int D, typename IteratorType>
   bool FilmIteratorBase<D, IteratorType>::isAthermal() const 
   {
      int nm = system().mixture().nMonomer();

      UTIL_CHECK(nm > 0);
      UTIL_CHECK(chi().capacity1() == nm);
      UTIL_CHECK(chi().capacity2() == 2);

      for (int i = 0; i < nm; i++) {
         if ((chi(i,0) != 0) || (chi(i,1) != 0)) {
            return false;
         }
      }
      return true;
   }

   /*
   * Return the value of phi, the overall volume fraction of unit cell that is
   * occupied by the walls.
   */
   template <int D, typename IteratorType>
   double const FilmIteratorBase<D, IteratorType>::phi() const
   {
      // Must call generateWallFields() before calling phi()
      UTIL_CHECK(wallCField().isAllocated());
      UTIL_CHECK(wallCField()[0] > 0);

      return wallCField()[0];
   }

   /*
   * Write the concentration field of the walls to an output stream, in basis format.
   */
   template <int D, typename IteratorType>
   void FilmIteratorBase<D, IteratorType>::writeWallField(std::ostream &out) const 
   {
      UTIL_CHECK(wallCField_.isAllocated());

      DArray< DArray <double> > cField;
      cField.allocate(1);
      cField[0] = wallCField_;
      system().domain().fieldIo().writeFieldsBasis(out, cField, system().domain().unitCell());
   }

   /*
   * Write the concentration field of the walls to a file specified by 'filename', 
   * in basis format.
   */
   template <int D, typename IteratorType>
   void FilmIteratorBase<D, IteratorType>::writeWallField(std::string filename) const 
   {
      UTIL_CHECK(wallCField_.isAllocated());

      DArray< DArray <double> > cField;
      cField.allocate(1);
      cField[0] = wallCField_;
      system().domain().fieldIo().writeFieldsBasis(filename, cField, system().domain().unitCell());
   }

   /*
   * Write the concentration field of the walls to an output stream, in r-grid format.
   */
   template <int D, typename IteratorType>
   void FilmIteratorBase<D, IteratorType>::writeWallFieldRGrid(std::ostream &out) const 
   {
      UTIL_CHECK(wallCField_.isAllocated());

      DArray< DArray <double> > cField;
      cField.allocate(1);
      cField[0] = wallCField_;

      DArray< RField<D> > rGrid;
      rGrid.allocate(1);
      rGrid[0].allocate(system().domain().mesh().dimensions());

      system().domain().fieldIo().convertBasisToRGrid(cField, rGrid);
      system().domain().fieldIo().writeFieldsRGrid(out, rGrid, system().domain().unitCell());
   }

   /*
   * Write the concentration field of the walls to a file specified by 'filename', 
   * in r-grid format.
   */
   template <int D, typename IteratorType>
   void FilmIteratorBase<D, IteratorType>::writeWallFieldRGrid(std::string filename) const 
   {
      UTIL_CHECK(wallCField_.isAllocated());

      DArray< DArray <double> > cField;
      cField.allocate(1);
      cField[0] = wallCField_;

      DArray< RField<D> > rGrid;
      rGrid.allocate(1);
      rGrid[0].allocate(system().domain().mesh().dimensions());

      system().domain().fieldIo().convertBasisToRGrid(cField, rGrid);
      system().domain().fieldIo().writeFieldsRGrid(filename, rGrid, system().domain().unitCell());
   }

} // namespace Pspc
} // namespace Pscf
#endif
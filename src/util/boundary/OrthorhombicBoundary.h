#ifndef UTIL_ORTHORHOMBIC_BOUNDARY_H
#define UTIL_ORTHORHOMBIC_BOUNDARY_H

/*
* Simpatico - Simulation Package for Polymeric and Molecular Liquids
*
* Copyright 2010 - 2014, The Regents of the University of Minnesota
* Distributed under the terms of the GNU General Public License.
*/

#include "OrthoRegion.h"                // base class
#include <util/crystal/LatticeSystem.h> // member
#include <util/containers/FSArray.h>    // member template
#include <util/space/Vector.h>          // member template argument
#include <util/space/IntVector.h>       // inline functions
#include <util/space/Dimension.h>       // member template argument
#include <util/global.h>                // asserts in inline functions

#include <cmath>
#include <iostream>

class OrthorhombicBoundaryTest;

namespace Util
{

   class Random;

   /**
   * An orthorhombic periodic unit cell.
   *
   * \ingroup Boundary_Module
   */
   class OrthorhombicBoundary : private OrthoRegion
   {

   public:

      /**
      * Constructor.
      */
      OrthorhombicBoundary();

      // Use default copy constructor.

      /**
      * Set unit cell dimensions for orthorhombic boundary.
      *
      * Also sets all related lengths and volume.
      *
      * \param lengths  Vector of unit cell lengths
      */
      void setOrthorhombic(const Vector &lengths);

      /**
      * Set unit cell dimensions for tetragonal boundary.
      *
      * \param a  unit cell dimensions in unique x-direction
      * \param bc unit cell length in y and z directions
      */
      void setTetragonal(double a, double bc);

      /**
      * Set unit cell dimensions for a cubic boundary.
      *
      * \param a unit cell length in x, y, and z directions.
      */
      void setCubic(double a);

      /**
      * Serialize to/from an archive.
      * 
      * \param ar       saving or loading archive
      * \param version  archive version id
      */
      template <class Archive>
      void serialize(Archive& ar, const unsigned int version);

      ///\name Periodic Boundary Conditions - Primitive Cell
      //@{

      /**
      * Shift Cartesian Vector r to its primary image.
      *
      * Upon completion, each component r[i] of Vector r is shifted by
      * a multiple of length[i] so as to lie within the allowed range 
      * minima_[i] <= r[i] < maxima_[i].
      *
      * Precondition: The algorithm assumes that on input, for each i,
      * minima_[i] - lengths_[i] < r[i] < maxima_[i] + lengths_[i]
      *
      * \param r Vector of Cartesian coordinates
      */
      void shift(Vector &r) const;

      /**
      * Shift Cartesian Vector r to its primary image.
      *
      * This function maps an atomic position to its primary image, and
      * also increments the atomic shift IntVector:
      *
      * If    r[i]     ->  r[i] - t*length_[i], 
      * then  shift[i] ->  shift[i] + t.
      *
      * \sa Atom:shift()
      *
      * \param r     Vector of Cartesian coordinates
      * \param shift integer shifts required to obtain "true" coordinates.
      */
      void shift(Vector &r, IntVector& shift) const;

      /**
      * Shift Cartesian Vector r by multiple t of a Bravais lattice vector.
      *
      * This function shifts the Vector r by a specified amount:
      *
      *   r ->  r + t*a'[i]
      *
      * where a[i] is Bravais lattice vector number i.
      *
      * \param r Cartesian position Vector 
      * \param i direction index
      * \param t multiple of Bravais lattice vector i
      */
      void applyShift(Vector &r, int i, int t) const;

      /**
      * Shift generalized Vector r to its primary image.
      *
      * One output, each coordinate r[i] is shifted by an integer, so as to
      * lie within the range 0 < r[i] < 1.0
      *
      * Precondition: The algorithm assumes that on input, for each i=0,..,2,
      * -1.0 < r[i] < 2.0
      *
      * \param r Vector of scaled / generalized coordinates
      */
      void shiftGen(Vector &r) const;

      /**
      * Shift generalized Vector r to its image within the primary unit cell.
      *
      * This function maps a vector of generalized coordinates to lie in the 
      * primary cell, and also increments the atomic shift IntVector:
      *
      * If r[i] -> r[i] - t, then shift[i] -> shift[i] + t.
      *
      * \sa Atom:shift()
      *
      * \param r     Vector of generalized coordinates  (in/out)
      * \param shift integer shifts, modified on output (in/out)
      */
      void shiftGen(Vector &r, IntVector& shift) const;

      //@}
      ///\name Minimum Image Pair Separations for Cartesian Vectors
      //@{

      /**
      * Return square distance between positions r1 and r2.
      *
      * This function returns the square distance between two Cartesian 
      * vectors, using the nearest image convention for the separation Vector.
      *
      * \param r1  first Cartesian position Vector
      * \param r2  second Cartesian position Vector
      * \return square of distance between r1 and r2, using nearest image.
      */
      double distanceSq(const Vector &r1, const Vector &r2) const;

      /**
      * Return square distance between positions r1 and r2.
      *
      * This function returns the square distance between two Cartesian 
      * vectors, using the nearest image convention for the separation Vector.
      * On return, the IntVector shift contains the coefficients of the
      * Bravais lattice vectors that were added to r1 to create a nearest
      * image of r2.
      *
      * \param r1  first Cartesian position Vector
      * \param r2  second Cartesian position Vector
      * \param shift  shift added to r1 to create nearest image of r2.
      * \return square of distance between r1 and r2, using nearest image.
      */
      double distanceSq(const Vector &r1, const Vector &r2, IntVector& shift) 
      const;

      /**
      * Compute distance and separation between r1 and r2.
      *
      * This function returns the square distance between two Cartesian 
      * vectors, using the nearest image convention for the separation Vector.
      * Upon return, Vector dr contains the separation r1 - r2 computed using
      * the nearest image convention, and the return value is the square of dr.
      *
      * \param r1 first Cartesian position Vector
      * \param r2 second Cartesian position Vector
      * \param dr separation Vector (upon return)
      * \return square of separation Vector dr
      */
      double distanceSq(const Vector &r1, const Vector &r2, Vector &dr) const;

      //@}
      ///\name Minimum Image Vector Tests
      //@{
     
      /**
      * Is a generalized separation vector a minimimum image of itself?
      *
      * This function returns true if the scaled separation dr has a Cartesian
      * norm less than or equal to that of any vector that can be produced by
      * adding a Bravais vector to the corresponding Cartesian vector, or false 
      * otherwise.
      *
      * \param dr separation vector in generalized coordinates.
      * \return true if minimum image, false otherwise
      */
      bool isMinImageGen(const Vector &dr);

      /**
      * Is a Cartesian separation vector a minimimum image of itself?
      *
      * This function returns true if the scaled separation dr has a Cartesian
      * norm less than or equal to that of any vector that can be produced by
      * adding a Bravais vector to dr, false otherwise.
      *
      * \param dr separation vector in Cartesian coordinates.
      * \return true if minimum image, false otherwise
      */
      bool isMinImageCart(const Vector &dr);

      //@}
      ///\name Coordinate Transformations
      //@{

      /**
      * Transform Cartesian Vector to scaled / generalized coordinates.
      *
      * Generalized coordinates range from 0.0 < Rg[i] < 1.0 within the
      * primitive cell, for i=0,..,2.
      *
      * \param Rc Vector of Cartesian coordinates (input)
      * \param Rg Vector of generalized coordinates (output)
      */
      void transformCartToGen(const Vector& Rc, Vector& Rg) const;

      /**
      * Transform Vector of generalized coordinates to Cartesian Vector.
      *
      * \param Rg Vector of generalized coordinates (input)
      * \param Rc Vector of Cartesian coordinates (output)
      */
      void transformGenToCart(const Vector& Rg, Vector& Rc) const;

      //@}
      ///\name Accessors
      //@{

      /**
      * Return actual lattice system.
      *
      * Value can be Cubic, Tetragonal, or Orthorhombic.
      */
      LatticeSystem latticeSystem();

      /**
      * Get Vector of unit cell lengths by const reference.
      */
      const Vector& lengths() const;

      /**
      * Get length in Cartesian direction i.
      *
      * \param i index of Cartesian direction, 0 <= i < Dimension
      */
      double length(int i) const;

      /**
      * Get minimum length across primitive unit cell.
      */
      double minLength() const;

      /**
      * Return unit cell volume.
      */
      double volume() const;

      /**
      * Return Bravais lattice vector i
      *  
      * \param i basis Vector index.
      */ 
      const Vector& bravaisBasisVector(int i) const;

      /**
      * Return reciprocal lattice basis vector i
      *  
      * \param i basis Vector index.
      */ 
      const Vector& reciprocalBasisVector(int i) const;

      /**
      * Generate random position within the primary unit cell.
      *
      * Generates Vector r[i], i=1,..,3 with minima_[i] < r[i] < maxima_[i].
      *
      * \param random random number generator object
      * \param r      Vector of random coordinates (upon return)
      */
      void randomPosition(Random &random, Vector &r) const;

      /**
      * Return true if valid, or throw Exception.
      */
      bool isValid();

      //@}

   private:

      /**
      * Array of Bravais lattice vectors.
      */
      FSArray<Vector, Dimension>  bravaisBasisVectors_;

      /**
      * Array of Reciprocal lattice vectors.
      */
      FSArray<Vector, Dimension>  reciprocalBasisVectors_;

      /**
      * Vector of inverse box dimensions.
      */
      Vector invLengths_;

      /**
      * Minimum distance across the unit cell.
      */
      double minLength_;

      /**
      * Actual lattice system (Orthorhombic, Tetragonal, or Cubic)
      */
      LatticeSystem lattice_;

   // friends:

      /// Unit test
      friend class ::OrthorhombicBoundaryTest;

      /// istream extractor
      friend std::istream& operator >> (std::istream& in, 
                                        OrthorhombicBoundary& boundary);

      /// ostream inserter
      friend std::ostream& operator << (std::ostream& out, 
                                        const OrthorhombicBoundary& boundary);

      /**
      * Reset all quantities that depend upon lengths.
      */ 
      void reset();

   };

   // Friend function declarations

   /**
   * istream extractor for a OrthorhombicBoundary.
   *
   * \param  in        input stream
   * \param  boundary  OrthorhombicBoundary to be read from stream
   * \return modified input stream
   */
   std::istream& operator >> (std::istream& in, OrthorhombicBoundary& boundary);

   /**
   * ostream inserter for an OrthorhombicBoundary.
   *
   * \param  out      output stream
   * \param  boundary OrthorhombicBoundary to be written to stream
   * \return modified output stream
   */
   std::ostream& operator << (std::ostream& out, 
                              const OrthorhombicBoundary& boundary);

   // Inline member function definitions

   /* 
   * Return Vector of lengths by const reference.
   */
   inline const Vector& OrthorhombicBoundary::lengths() const 
   {  return lengths_; }

   /* 
   * Get length = maximum - minimum in direction i.
   */
   inline double OrthorhombicBoundary::length(int i) const 
   {  return lengths_[i]; }

   /* 
   * Return the maximum validity range of the distances.
   */
   inline double OrthorhombicBoundary::minLength() const
   {  return minLength_; }

   /* 
   * Return region volume.
   */
   inline double OrthorhombicBoundary::volume() const
   {  return volume_; }

   /* 
   * Return Bravais lattice basis vector number i.
   */
   inline const Vector& OrthorhombicBoundary::bravaisBasisVector(int i) const
   {  return bravaisBasisVectors_[i]; }

   /* 
   * Return reciprocal lattice basis vector number i.
   */
   inline const Vector& OrthorhombicBoundary::reciprocalBasisVector(int i) const
   {  return reciprocalBasisVectors_[i]; }

   /* 
   * Shift Cartesian Vector r to primitive unit cell.
   */
   inline void OrthorhombicBoundary::shift(Vector& r) const
   {
      for (int i = 0; i < Dimension; ++i) {
         if( r[i] >= maxima_[i] ) {
            r[i] = r[i] - lengths_[i];
            assert(r[i] < maxima_[i]);
         } else
         if ( r[i] <  minima_[i] ) {
            r[i] = r[i] + lengths_[i];
            assert(r[i] >= minima_[i]);
         }
      }
   }

   /* 
   * Shift Vector r to periodic cell, R[axis] < r[axis] < maxima_[axis].
   */
   inline void OrthorhombicBoundary::shift(Vector& r, IntVector& shift) const
   {
      for (int i = 0; i < Dimension; ++i) {
         if (r[i] >= maxima_[i]) {
            r[i] = r[i] - lengths_[i];
            ++(shift[i]);
            assert(r[i] < maxima_[i]);
         } else
         if (r[i] <  minima_[i]) {
            r[i] = r[i] + lengths_[i];
            --(shift[i]);
            assert(r[i] >= minima_[i]);
         }
      }
   }

   /*
   * Shift Cartesian Vector r by multiple t of a Bravais lattice vector.
   */
   inline void OrthorhombicBoundary::applyShift(Vector &r, int i, int t) const
   {  r[i] += t*lengths_[i]; }

   /* 
   * Shift generalized Vector r to primitive unit cell.
   */
   inline void OrthorhombicBoundary::shiftGen(Vector& r) const
   {
      for (int i = 0; i < Dimension; ++i) {
         if( r[i] >= 1.0 ) {
            r[i] = r[i] - 1.0;
            assert(r[i] < 1.0);
         } else
         if ( r[i] <  0.0 ) {
            r[i] = r[i] + 1.0;
            assert(r[i] >= 0.0);
         }
      }
   }

   /* 
   * Shift generalized Vector r to primitive cell, 0 < r[axis] < 1.0.
   */
   inline 
   void OrthorhombicBoundary::shiftGen(Vector& r, IntVector& shift) const
   {
      for (int i = 0; i < Dimension; ++i) {
         if (r[i] >= 1.0) {
            r[i] = r[i] - 1.0;
            ++(shift[i]);
            assert(r[i] < 1.0);
         } else
         if (r[i] <  0.0) {
            r[i] = r[i] + 1.0;
            --(shift[i]);
            assert(r[i] >= 0.0);
         }
      }
   }

   /* 
   * Calculate squared distance by minimum image convention.
   */
   inline 
   double OrthorhombicBoundary::distanceSq(const Vector &r1, const Vector &r2, 
                               IntVector& translate) const
   {
      double dr;
      double norm = 0.0;
      for (int i = 0; i < Dimension; ++i) {
         dr = r1[i] - r2[i];
         if ( fabs(dr) > halfLengths_[i] ) {
            if (dr > 0.0) {
               dr -= lengths_[i];
               translate[i] = -1;
            } else {
               dr += lengths_[i];
               translate[i] = +1;
            }
            assert(fabs(dr) <= halfLengths_[i]);
         } else {
            translate[i] = 0;
         }
         norm += dr*dr;
      }
      return norm;
   }

   /* 
   * Return squared distance and separation with minimum image convention.
   */
   inline 
   double OrthorhombicBoundary::distanceSq(const Vector &r1, const Vector &r2) const
   {
      double dr;
      double norm = 0.0;
      for (int i = 0; i < Dimension; ++i) {
         dr = r1[i] - r2[i];
         if (fabs(dr) > halfLengths_[i]) {
            if (dr > 0.0) {
               dr -= lengths_[i];
            } else {
               dr += lengths_[i];
            }
            assert(fabs(dr) <= halfLengths_[i]);
         }
         norm += dr*dr;
      }
      return norm;
   }

   /* 
   * Calculate squared distance between two positions, and separation vector,
   * using the minimum image convention.
   */
   inline 
   double OrthorhombicBoundary::distanceSq(const Vector &r1, const Vector &r2, 
                                           Vector &dr) const
   {
      for (int i = 0; i < Dimension; ++i) {
         dr[i] = r1[i] - r2[i];
         if (fabs(dr[i]) > halfLengths_[i]) {
            if (dr[i] > 0.0) {
               dr[i] -= lengths_[i];
            } else {
               dr[i] += lengths_[i];
            }
            assert(fabs(dr[i]) <= halfLengths_[i]);
         }
      }
      return ( dr[0]*dr[0] + dr[1]*dr[1] + dr[2]*dr[2] );
   }

   /*
   * Return actual lattice system.
   */
   inline LatticeSystem OrthorhombicBoundary::latticeSystem()
   { return lattice_; }

   /*
   * Transform Cartesian Vector Rc to generalized Vector Rg.
   */
   inline void 
   OrthorhombicBoundary::transformCartToGen(const Vector& Rc, Vector& Rg) 
   const
   {
      for (int i = 0; i < Dimension; ++i) {
         Rg[i] = Rc[i] * invLengths_[i];
      }
   }
      
   /*
   * Transform Cartesian Vector Rc to generalized Vector Rg.
   */
   inline void 
   OrthorhombicBoundary::transformGenToCart(const Vector& Rg, Vector& Rc) 
   const
   {
      for (int i = 0; i < Dimension; ++i) {
         Rc[i] = Rg[i] * lengths_[i];
      }
   }

   /*
   * Is a generalized separation vector a minimimum image of itself?
   */
   inline
   bool OrthorhombicBoundary::isMinImageGen(const Vector &dr)
   {
      for (int i = 0; i < Dimension; ++i) {
         if (fabs(dr[i]) > 0.5) {
            return false;
         }
      }
      return true;
   }

   /*
   * Is Cartesian separation vector dr a minimimum image of itself?
   */
   inline
   bool OrthorhombicBoundary::isMinImageCart(const Vector &dr)
   {
      for (int i = 0; i < Dimension; ++i) {
         if (fabs(dr[i]) > halfLengths_[i]) {
            return false;
         }
      }
      return true;
   }

   // Member function template

   /*
   * Serialize an OrthorhombicBoundary to/from an archive.
   */
   template <class Archive> void 
   OrthorhombicBoundary::serialize(Archive& ar, const unsigned int version)
   {
      OrthoRegion::serialize(ar, version);
      ar & bravaisBasisVectors_;
      ar & reciprocalBasisVectors_;
      ar & invLengths_;
      ar & minLength_;
      serializeEnum(ar, lattice_, version);
      if (Archive::is_loading()) {
         isValid();
      }
   }

}
 
#ifdef UTIL_MPI
#include <util/mpi/MpiSendRecv.h>
#include <util/mpi/MpiTraits.h>

namespace Util
{

   /**
   * Send an OrthorhombicBoundary via MPI.
   */
   template <>
   void send<Util::OrthorhombicBoundary>(MPI::Comm& comm, 
             Util::OrthorhombicBoundary& data, int dest, int tag);

   /**
   * Receive an OrthorhombicBoundary via MPI.
   */
   template <>
   void recv<Util::OrthorhombicBoundary>(MPI::Comm& comm, 
             Util::OrthorhombicBoundary& data, int source, int tag);

   /**
   * Broadcast an OrthorhombicBoundary via MPI.
   */
   template <>
   void bcast<Util::OrthorhombicBoundary>(MPI::Intracomm& comm, 
              Util::OrthorhombicBoundary& data, int root);

   /**
   * Explicit specialization MpiTraits<OrthorhombicBoundary>.
   */
   template <>
   class MpiTraits<Util::OrthorhombicBoundary>
   {
   public:
      static MPI::Datatype type;         ///< MPI Datatype
      static bool hasType;               ///< Is the MPI type initialized?
   };

}
#endif // ifdef  UTIL_MPI

#endif // ifndef UTIL_ORTHORHOMBIC_BOUNDARY_H

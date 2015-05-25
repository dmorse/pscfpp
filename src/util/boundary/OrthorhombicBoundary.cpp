/*
* Simpatico - Simulation Package for Polymeric and Molecular Liquids
*
* Copyright 2010 - 2014, The Regents of the University of Minnesota
* Distributed under the terms of the GNU General Public License.
*/

#include "OrthorhombicBoundary.h"
#include <util/space/Vector.h>
#include <util/space/Dimension.h>
#include <util/random/Random.h>
#include <util/math/Constants.h>
#include <util/math/feq.h>
#include <util/format/Dbl.h>
#include <util/global.h>

namespace Util
{

   /*
   * Default constructor.
   */
   OrthorhombicBoundary::OrthorhombicBoundary()
    : OrthoRegion(),
      lattice_(Orthorhombic)
   {
      for (int i = 0; i < Dimension; ++i) {
         invLengths_[i] = 1.0/lengths_[i];
         bravaisBasisVectors_.append(Vector::Zero);
         bravaisBasisVectors_[i][i] = lengths_[i];
         reciprocalBasisVectors_.append(Vector::Zero);
         reciprocalBasisVectors_[i][i] = 2.0*Constants::Pi/lengths_[i];
      }
      minLength_ = lengths_[0];
      for (int i = 1; i < Dimension; ++i) {
         if (lengths_[i] < minLength_) {
            minLength_ = lengths_[i];
         }
      }
   }

   /*
   * Set box lengths and then call reset.
   */
   void OrthorhombicBoundary::setOrthorhombic(const Vector &lengths)
   {
      maxima_  = lengths;
      lattice_ = Orthorhombic;
      reset();
   }

   /*
   * Set box lengths and then call reset.
   */
   void OrthorhombicBoundary::setTetragonal(double a, double bc) 
   {  
      maxima_[0] = a;
      maxima_[1] = bc;
      maxima_[2] = bc;
      lattice_ = Tetragonal;
      reset();
   }

   /*
   * Set box lengths and call reset.
   */
   void OrthorhombicBoundary::setCubic(double a)
   {
      maxima_[0] = a;
      maxima_[1] = a;
      maxima_[2] = a;
      lattice_ = Cubic;
      reset();
   }

   /*
   * Reset all quantities that depend on unit cell lengths.
   */
   void OrthorhombicBoundary::reset()
   {
      resetRegion();
      for (int i = 0; i < Dimension; ++i) {
         invLengths_[i] = 1.0/lengths_[i];
         bravaisBasisVectors_[i][i] = lengths_[i];
         reciprocalBasisVectors_[i][i] = 2.0*Constants::Pi/lengths_[i];
      }

      minLength_ = lengths_[0];
      for (int i = 1; i < Dimension; ++i) {
         if (lengths_[i] < minLength_) {
            minLength_ = lengths_[i];
         }
      }

   }

   /*
   * Generate a random position within the box
   *
   * \param random random number generator object
   * \param r      Vector of random coordinates
   */
   void OrthorhombicBoundary::randomPosition(Random &random, Vector &r) const
   {
     for (int i=0; i < Dimension; ++i) {
        r[i] = random.uniform(minima_[i], maxima_[i]);
     }
   }

   /*
   * Check consistency of data.
   */
   bool OrthorhombicBoundary::isValid()
   {
      double dot;
      double twoPi = 2.0*Constants::Pi;
      int i, j;

      OrthoRegion::isValid();
      for (i = 0; i < Dimension; ++i) {
         if (!feq(minima_[i], 0.0)) {
            UTIL_THROW("minima_[i] != 0");
         }
         if (!feq(lengths_[i]*invLengths_[i], 1.0)) {
            UTIL_THROW("invLengths_[i]*lengths_[i] != 1.0");
         }
         if (!feq(lengths_[i], bravaisBasisVectors_[i][i])) {
            UTIL_THROW("bravaisBasisVectors_[i][i] != lengths_[i]");
         }
         for (j = 0; j < Dimension; ++j) {
            dot = bravaisBasisVectors_[i].dot(reciprocalBasisVectors_[j]);
            if (i == j) {
               if (!feq(dot, twoPi)) {
                  UTIL_THROW("a[i].b[i] != twoPi");
               }
            } else {
               if (!feq(dot, 0.0)) {
                  UTIL_THROW("a[i].b[j] != 0 for i != j");
               }
               if (!feq(bravaisBasisVectors_[i][j], 0.0)) {
                  UTIL_THROW("Nonzero off-diagonal element of bravaisBasisVectors_");
               }
               if (!feq(reciprocalBasisVectors_[i][j], 0.0)) {
                  UTIL_THROW("Nonzero off-diagonal element of reciprocalBasisVectors_");
               }
            }
         }
      }
      return true;
   }

   /*
   * Input a OrthorhombicBoundary from an istream, without line breaks.
   */
   std::istream& operator >> (std::istream& in, OrthorhombicBoundary &boundary)
   {
      LatticeSystem lattice;
      in >> lattice;
      if (lattice == Orthorhombic) {
         in >> boundary.maxima_;
      } else
      if (lattice == Tetragonal) {
         double a, bc;
         in >> a >> bc;
         boundary.maxima_[0] = a;
         boundary.maxima_[1] = bc;
         boundary.maxima_[2] = bc;
      } else 
      if (lattice == Cubic) {
         double a;
         in >> a;
         boundary.maxima_[0] = a;
         boundary.maxima_[1] = a;
         boundary.maxima_[2] = a;
      } else {
         UTIL_THROW("Lattice must be orthorhombic, tetragonal or cubic");
      }
      boundary.lattice_ = lattice;
      boundary.reset();
      return in;
   }

   /*
   * Output an OrthorhombicBoundary to an ostream, without line breaks.
   */
   std::ostream&
   operator << (std::ostream& out, const OrthorhombicBoundary &boundary)
   {
      out << boundary.lattice_ << "   ";
      if (boundary.lattice_ == Orthorhombic) {
         out << boundary.lengths_;
      } else
      if (boundary.lattice_ == Tetragonal) {
         out << Dbl(boundary.lengths_[0]);
         out << Dbl(boundary.lengths_[2]);
      } else
      if (boundary.lattice_ == Cubic) {
         out << Dbl(boundary.lengths_[0]);
      }
      return out;
   }

   #ifdef UTIL_MPI
   template <>
   void send<Util::OrthorhombicBoundary>(MPI::Comm& comm,
             Util::OrthorhombicBoundary& data, int dest, int tag)
   {
      Vector lengths = data.lengths();
      send<Vector>(comm, lengths, dest, tag);
   }

   template <>
   void recv<Util::OrthorhombicBoundary>(MPI::Comm& comm,
             Util::OrthorhombicBoundary& data, int source, int tag)
   {
      Vector lengths;
      recv<Vector>(comm, lengths, source, tag);
      data.setOrthorhombic(lengths);
   }

   template <>
   void bcast<Util::OrthorhombicBoundary>(MPI::Intracomm& comm,
              Util::OrthorhombicBoundary& data, int root)
   {
      Vector lengths;
      int    rank = comm.Get_rank();
      if (rank == root) {
         lengths = data.lengths();
      }
      bcast<Vector>(comm, lengths, root);
      if (rank != root) {
         data.setOrthorhombic(lengths);
      }
   }

   /*
   * Initialize MPI Datatype.
   */
   MPI::Datatype MpiTraits<Util::OrthorhombicBoundary>::type = MPI::BYTE;
   bool MpiTraits<Util::OrthorhombicBoundary>::hasType = false;
   #endif

}

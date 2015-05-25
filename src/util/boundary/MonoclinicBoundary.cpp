/*
* Simpatico - Simulation Package for Polymeric and Molecular Liquids
*
* Copyright 2010, The Regents of the University of Minnesota
* Distributed under the terms of the GNU General Public License.
*/

#include "MonoclinicBoundary.h"
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
   MonoclinicBoundary::MonoclinicBoundary() 
    : lattice_(Monoclinic)
   {
      double twoPi = 2.0*Constants::Pi;
      for (int i = 0; i < Dimension; ++i) {
         minima_[i] = 0.0;
         maxima_[i] = 1.0;
         l_[i] = 1.0;
         invL_[i] = 1.0;
         halfL_[i] = 0.5;
         lengths_[i] = 1.0;
         bravaisBasisVectors_.append(Vector::Zero);
         bravaisBasisVectors_[i][i] = l_[i];
         reciprocalBasisVectors_.append(Vector::Zero);
         reciprocalBasisVectors_[i][i] = twoPi / l_[i];
      }
      d_ = 0.0;
      volume_ = 1.0;
         
      c3_ = 0.0;
      e_ = 1.0;
      minLength_ = 1.0;
   }

   /* 
   * Set unit cell parameters and then call reset.
   */
   void MonoclinicBoundary::setMonoclinic(const Vector &lengths, const double d) 
   {
      l_ = lengths;
      d_ = d;
      reset(); 
   }

   /* 
   * Setup orthorhombic unit cell.
   */
   void MonoclinicBoundary::setOrthorhombic(const Vector &lengths) 
   {
      l_ = lengths;
      d_ = 0.0;
      reset(); 
   }

   /* 
   * Reset all quantities that depend on unit cell lengths.
   */
   void MonoclinicBoundary::reset()
   {
      double twoPi = 2.0*Constants::Pi;
      for (int i = 0; i < Dimension; ++i) {
         minima_[i] = 0.0;
         maxima_[i] = l_[i];
         halfL_[i] = 0.5*l_[i];
         invL_[i] = 1.0/l_[i];
         bravaisBasisVectors_[i][i] = l_[i];
         reciprocalBasisVectors_[i][i] = twoPi/l_[i];
      }
      bravaisBasisVectors_[1][2] = d_;
      reciprocalBasisVectors_[2][1] = -twoPi*d_/(l_[1]*l_[2]);

      volume_ = l_[0] * l_[1] * l_[2];
      e_ = sqrt(l_[1]*l_[1] + d_*d_);          
      c3_ = -d_/l_[1];

      lengths_[0] = l_[0];
      lengths_[1] = l_[1];
      lengths_[2] = l_[2] / sqrt(1.0 + c3_*c3_);

      minLength_ = lengths_[0];
      for (int i = 1; i < Dimension; ++i) {
         if (lengths_[i] < minLength_) {
            minLength_ = lengths_[i];
         }
      }
   }

   /* 
   * Generate a random Cartesian position within the primitive unit cell.
   *
   * \param random random number generator object
   * \param r      Vector of random Cartesian coordinates
   */
   void MonoclinicBoundary::randomPosition(Random &random, Vector &r) const 
   {
     Vector Rg;
     for (int i=0; i < Dimension; ++i) {
        Rg[i] = random.uniform(0.0, 1.0);
     }
     transformGenToCart(Rg, r);
   }

   /* 
   * Check consistency of data.
   */
   bool MonoclinicBoundary::isValid() 
   {  
      double dot;
      double twoPi = 2.0*Constants::Pi;
      int i, j;
      for (i = 0; i < Dimension; ++i) {
         if (maxima_[i] <= minima_[i]) {
            UTIL_THROW("maxima_[i] <= minima_[i]");
         }
         if (!feq(l_[i], maxima_[i] - minima_[i])) {
            UTIL_THROW("l_[i] != maxima_[i] - minima_[i]");
         }
         if (!feq(halfL_[i], 0.5*l_[i])) {
            UTIL_THROW("halfL_[i] != 0.5*l_[i]");
         }
         if (!feq(minima_[i], 0.0)) {
            UTIL_THROW("minima_[i] != 0");
         }
         if (!feq(l_[i]*invL_[i], 1.0)) {
            UTIL_THROW("invL_[i]*l_[i] != 1.0");
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
            }
         }
      }
      if (!feq(volume_, l_[0]*l_[1]*l_[2])) {
         UTIL_THROW("volume_ != product of l_");
      }
      return true;
   }

   /* 
   * Input a MonoclinicBoundary from an istream, without line breaks.
   */
   std::istream& operator >> (std::istream& in, MonoclinicBoundary &boundary)
   {
      LatticeSystem lattice;
      in >> lattice;
      if (lattice != Monoclinic) {
         UTIL_THROW("Lattice must be Monoclinic");
      }
      in >> boundary.l_;
      in >> boundary.d_;
      boundary.reset();
      return in;
   }

   /* 
   * Output an MonoclinicBoundary to an ostream, without line breaks.
   */
   std::ostream& 
   operator << (std::ostream& out, const MonoclinicBoundary &boundary) 
   {
      out << boundary.lattice_ << "   ";
      out << boundary.l_ << "   ";
      out << boundary.d_;
      return out;
   }
 
   #ifdef UTIL_MPI
   template <>
   void send<MonoclinicBoundary>(MPI::Comm& comm, 
             MonoclinicBoundary& data, int dest, int tag)
   {
      send<Vector>(comm, data.l_, dest, tag);
      send<double>(comm, data.d_, dest, tag + 386);
   }

   template <>
   void recv<MonoclinicBoundary>(MPI::Comm& comm, 
             MonoclinicBoundary& data, int source, int tag)
   {
      Vector l;
      double d;
      recv<Vector>(comm, l, source, tag);
      recv<double>(comm, d, source, tag + 386);
      data.setMonoclinic(l, d);
   }

   template <>
   void bcast<MonoclinicBoundary>(MPI::Intracomm& comm, 
              MonoclinicBoundary& data, int root)
   {
      Vector l; 
      double d;
      int rank = comm.Get_rank();
      if (rank == root) {
         l = data.l_;
         d = data.d_;
      }
      bcast<Vector>(comm, l, root);
      bcast<double>(comm, d, root);
      if (rank != root) {
         data.setMonoclinic(l, d);
      }
   }

   /*
   * Initialize MPI Datatype.
   */
   MPI::Datatype MpiTraits<MonoclinicBoundary>::type = MPI::BYTE;
   bool MpiTraits<MonoclinicBoundary>::hasType = false;
   #endif

}

/*
* PSCF - Polymer Self-Consistent Field Theory
*
* Copyright 2016, The Regents of the University of Minnesota
* Distributed under the terms of the GNU General Public License.
*/

#include "Molecule.h"  

namespace Pscf{ 

   /*
   * Constructor.
   */
   Homogeneous::Molecule::Molecule()
    : clumps_(),
      nClump_(0),
      size_(0.0)
   {  setClassName("Molecule"); }

   /*
   * Destructor.
   */
   Homogeneous::Molecule::~Molecule()
   {}

   void Homogeneous::Molecule::readParameters(std::istream& in)
   {
      UTIL_ASSERT(clumps_.capacity() == 0);

      read<int>(in, "nClump", nClump_);

      // Allocate all arrays
      clumps_.allocate(nClump_);

      readDArray<Homogeneous::Clump>(in, "clumps", clumps_, nClump_);

      computeSize();
      // Read ensemble and phi or mu
      ensemble_ = Species::Closed;
      readOptional<Species::Ensemble>(in, "ensemble", ensemble_);
      if (ensemble_ == Species::Closed) {
         read(in, "phi", phi_);
      } else {
         read(in, "mu", mu_);
      }

   }

   /**
   * Allocate memory for specified number of clumps.
   */
   void Homogeneous::Molecule::setNClump(int nClump)
   {
      UTIL_ASSERT(clumps_.capacity() == 0);

      nClump_ = nClump;
      clumps_.allocate(nClump_);
   }

   void Homogeneous::Molecule::setPhi(double phi)
   {  phi_ = phi; }

   /*
   * Compute molecular size, by adding all clump sizes.
   */
   void Homogeneous::Molecule::computeSize()
   {
      UTIL_ASSERT(clumps_.capacity() > 0);
      UTIL_ASSERT(clumps_.capacity() == nClump_);

      size_ = 0.0;
      for (int clumpId = 0; clumpId < nClump_; ++clumpId) {
         size_ += clumps_[clumpId].size();
      }
   }

   #if 0
   /*
   * Compute solution to MDE and concentrations.
   */ 
   void Homogeneous::Molecule::compute(const DArray<WField>& wFields)
   {
      UTIL_ASSERT(clumps_.capacity() > 0);
      UTIL_ASSERT(clumps_.capacity() == nClump_);

      // Compute molecular partition function
      double q = clump(0).propagator(0).computeQ();
      if (ensemble() == Species::Closed) {
         mu_ = log(phi_/q);
      } 
      else if (ensemble() == Species::Open) {
         phi_ = exp(mu_)*q;
      }

      // Compute clump concentration fields
      double prefactor = phi_ / (q*size());
      for (int i = 0; i < nClump(); ++i) {
         clump(i).computeConcentration(prefactor);
      }

   }
   #endif
 
}

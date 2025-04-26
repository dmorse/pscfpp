/*
* PSCF - Polymer Self-Consistent Field Theory
*
* Copyright 2015 - 2025, The Regents of the University of Minnesota
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
      size_(0.0),
      hasSize_(false)
   {  setClassName("Molecule"); }

   /*
   * Destructor.
   */
   Homogeneous::Molecule::~Molecule()
   {}

   /*
   * Read chemical composition from file. 
   */
   void Homogeneous::Molecule::readParameters(std::istream& in)
   {
      UTIL_ASSERT(clumps_.capacity() == 0);

      read<int>(in, "nClump", nClump_);

      // Allocate all arrays
      clumps_.allocate(nClump_);

      readDArray<Homogeneous::Clump>(in, "clumps", clumps_, nClump_);
      computeSize();
   }

   /*
   * Allocate memory for specified number of clumps.
   */
   void Homogeneous::Molecule::setNClump(int nClump)
   {
      UTIL_ASSERT(clumps_.capacity() == 0);

      nClump_ = nClump;
      clumps_.allocate(nClump_);
   }

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
      hasSize_ = true;
   }
 
}

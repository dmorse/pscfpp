/*
* PSCF - Polymer Self-Consistent Field
*
* Copyright 2015 - 2026, The Regents of the University of Minnesota
* Distributed under the terms of the GNU General Public License.
*/

#include "FhMolecule.h"  

namespace Pscf{ 

   /*
   * Constructor.
   */
   FhMolecule::FhMolecule()
    : clumps_(),
      nClump_(0),
      size_(0.0),
      hasSize_(false)
   {  setClassName("FhMolecule"); }

   /*
   * Destructor.
   */
   FhMolecule::~FhMolecule()
   {}

   /*
   * Read chemical composition from file. 
   */
   void FhMolecule::readParameters(std::istream& in)
   {
      UTIL_ASSERT(clumps_.capacity() == 0);

      read<int>(in, "nClump", nClump_);

      // Allocate all arrays
      clumps_.allocate(nClump_);

      readDArray<FhClump>(in, "clumps", clumps_, nClump_);
      computeSize();
   }

   /*
   * Allocate memory for specified number of clumps.
   */
   void FhMolecule::setNClump(int nClump)
   {
      UTIL_ASSERT(clumps_.capacity() == 0);

      nClump_ = nClump;
      clumps_.allocate(nClump_);
   }

   /*
   * Compute molecular size, by adding all clump sizes.
   */
   void FhMolecule::computeSize()
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

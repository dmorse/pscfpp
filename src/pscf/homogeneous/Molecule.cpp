/*
* PSCF - Polymer Self-Consistent Field Theory
*
* Copyright 2016, The Regents of the University of Minnesota
* Distributed under the terms of the GNU General Public License.
*/

#include "Molecule.h"  

namespace Pscf{ 
namespace Homogeneous{ 

   /*
   * Constructor.
   */
   Molecule::Molecule()
    : groups_(),
      nGroup_(0),
      size_(0.0)
   {  setClassName("Homogeneous::Molecule"); }

   /*
   * Destructor.
   */
   Molecule::~Molecule()
   {}

   void Molecule::readParameters(std::istream& in)
   {
      read<int>(in, "nGroup", nGroup_);

      // Allocate all arrays
      groups_.allocate(nGroup_);

      readDArray<Group>(in, "groups", groups_, nGroup_);

      // Compute molecular size / volume
      size_ = 0.0;
      for (int groupId = 0; groupId < nGroup_; ++groupId) {
         size_ += groups_[groupId].size();
      }

      // Read ensemble and phi or mu
      ensemble_ = Species::Closed;
      readOptional<Species::Ensemble>(in, "ensemble", ensemble_);
      if (ensemble_ == Species::Closed) {
         read(in, "phi", phi_);
      } else {
         read(in, "mu", mu_);
      }

   }

   #if 0
   /*
   * Compute solution to MDE and concentrations.
   */ 
   void Molecule::compute(const DArray<WField>& wFields)
   {

      // Setup solvers for all groups
      int monomerId;
      for (int j = 0; j < nGroup(); ++j) {
         monomerId = group(j).monomerId();
         group(j).setupSolver(wFields[monomerId]);
      }

      // Compute molecular partition function
      double q = group(0).propagator(0).computeQ();
      if (ensemble() == Species::Closed) {
         mu_ = log(phi_/q);
      } 
      else if (ensemble() == Species::Open) {
         phi_ = exp(mu_)*q;
      }

      // Compute group concentration fields
      double prefactor = phi_ / (q*size());
      for (int i = 0; i < nGroup(); ++i) {
         group(i).computeConcentration(prefactor);
      }

   }
   #endif
 
}
}

#ifndef PSCF_SOLVENT_DESCRIPTOR_TPP
#define PSCF_SOLVENT_DESCRIPTOR_TPP

/*
* PSCF - Polymer Self-Consistent Field Theory
*
* Copyright 2016 - 2022, The Regents of the University of Minnesota
* Distributed under the terms of the GNU General Public License.
*/

#include "SolventDescriptor.h"

namespace Pscf {

   SolventDescriptor::SolventDescriptor()
   {  setClassName("SolventDescriptor"); }

   SolventDescriptor::~SolventDescriptor()
   {}

   void SolventDescriptor::readParameters(std::istream& in)
   {
      read<int>(in, "monomerId", monomerId_);
      read<double>(in, "size", size_);

      // Read phi or mu (but not both)
      bool hasPhi = readOptional(in, "phi", phi_).isActive();
      if (hasPhi) {
         ensemble_ = Species::Closed;
      } else {
         ensemble_ = Species::Open;
         read(in, "mu", mu_);
      }

      #if 0
      ensemble_ = Species::Closed;
      readOptional<Species::Ensemble>(in, "ensemble", ensemble_);
      if (ensemble_ == Species::Closed) {
         read(in, "phi", phi_);
      } else {
         read(in, "mu", mu_);
      }
      #endif
   }

   void SolventDescriptor::setPhi(double phi)
   {
      UTIL_CHECK(ensemble() == Species::Closed);  
      UTIL_CHECK(phi >= 0.0);  
      UTIL_CHECK(phi <= 1.0);  
      phi_ = phi;
   }

   void SolventDescriptor::setMu(double mu)
   {
      UTIL_CHECK(ensemble() == Species::Open);  
      mu_ = mu; 
   }

   /*
   * Set the id for this solvent.
   */ 
   void SolventDescriptor::setMonomerId(int monomerId)
   {  monomerId_ = monomerId; }
  
   /*
   * Set the id for this solvent.
   */ 
   void SolventDescriptor::setSize(double size)
   {  size_ = size; }

}
#endif

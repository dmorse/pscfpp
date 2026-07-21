#ifndef CPC_SOLVENT_TPP
#define CPC_SOLVENT_TPP

/*
* PSCF - Polymer Self-Consistent Field
*
* Copyright 2015 - 2025, The Regents of the University of Minnesota
* Distributed under the terms of the GNU General Public License.
*/

#include "Solvent.h"
#include <pscf/mesh/Mesh.h>
#include <pscf/cpu/complex.h>

namespace Pscf {
namespace Cpc { 

   /*
   * Constructor
   */
   template <int D>
   Solvent<D>::Solvent()
    : SolventSpecies< std::complex<double> >(),
      meshPtr_(nullptr)
   {  ParamComposite::setClassName("Solvent"); }

   /*
   * Destructor
   */
   template <int D>
   Solvent<D>::~Solvent()
   {}

   /*
   * Create an association with a Mesh.
   */
   template <int D>
   void Solvent<D>::associate(Mesh<D> const & mesh)
   {  meshPtr_ = &mesh; }

   /*
   * Allocate the concentration field (cField).
   */
   template <int D>
   void Solvent<D>::allocate()
   {
      UTIL_CHECK(meshPtr_);
      cField_.allocate(meshPtr_->dimensions());
   }

   /*
   * Compute concentration, q, and phi or mu.
   */ 
   template <int D>
   void Solvent<D>::compute(CField<D, CppTp<D> > const & wField)
   {
      int nx = meshPtr_->size(); // Number of grid points
      UTIL_CHECK(nx > 0);
      UTIL_CHECK(cField_.capacity() == nx);
      UTIL_CHECK(wField.capacity() == nx);

      double r;

      // Initialize cField_ to zero
      r = 0.0;
      for (int i = 0; i < nx; ++i) {
         assign(cField_[i], r);
      }
      fftw_complex Q;
      assign(Q, r);

      // Evaluate unnormalized integral and Q
      r = - size();
      fftw_complex z;
      for (int i = 0; i < nx; ++i) {
         mul(z, wField[i], r);
         assignExp(cField_[i], z); 
         addEq(Q, cField_[i]);
         // cField_[i] = exp(-size*wField[i]);
         // Q += cField_[i];
      }
      r = (double)nx;
      divEq(Q, r);
      // Q /= (double)nx

      // Set Q in Species base class and compute mu or phi 
      std::complex<double> Qstd;
      assign(Qstd, Q);
      Species< std::complex<double> >::setQ(Qstd);

      // Normalize concentration 
      fftw_complex prefactor;
      assign(z, phi());
      div(prefactor, z, Q);
      // prefactor = phi()/Q
      for (int i = 0; i < nx; ++i) {
          mulEq(cField_[i], prefactor);
          //cField_[i] *= prefactor;
      }
 
   }

} // namespace Cpc
} // namespace Rpc
#endif

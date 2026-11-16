#ifndef RP_SOLVENT_TPP
#define RP_SOLVENT_TPP

/*
* PSCF - Polymer Self-Consistent Field
*
* Copyright 2015 - 2026, David Morse (morse012@um.edu) and co-authors
* Distributed under the terms of the GNU General Public License.
*/

#include "Solvent.h"
#include <pscf/field/RField.h>
#include <pscf/backend/Reduce.h>
#include <pscf/backend/VecOp.h>
#include <pscf/mesh/Mesh.h>

#include <cmath>

namespace Pscf {
namespace Rp {

   /*
   * Constructor.
   */
   template <int D, class T>
   Solvent<D,T>::Solvent()
    : meshPtr_(nullptr)
   {  ParamComposite::setClassName("Solvent"); }

   /*
   * Destructor.
   */
   template <int D, class T>
   Solvent<D,T>::~Solvent()
   {}

   /*
   * Create an association with a mesh.
   */
   template <int D, class T>
   void Solvent<D,T>::associate(Mesh<D> const & mesh)
   {
      UTIL_CHECK(!meshPtr_);
      UTIL_CHECK(mesh.size() > 1);
      meshPtr_ = &mesh;
   }

   /*
   * Allocate memory for the concentration field (cField).
   */
   template <int D, class T>
   void Solvent<D,T>::allocate()
   {
      UTIL_CHECK(meshPtr_);
      cField_.allocate(meshPtr_->dimensions());
   }

   /*
   * Compute concentration, q, and phi or mu.
   */
   template <int D, class T>
   void Solvent<D,T>::compute(RField<D,T> const & wField, 
                              double phiTot)
   {
      // Local constants
      const int nx = meshPtr_->size();
      const double size = SolventSpeciesT::size();

      // Evaluate the Boltzmann weight of the field minus its minimum
      // wMin, which lies between 0 and 1 and is 1 at the minimum, so
      // that it can neither overflow nor vanish everywhere for a large
      // solvent size, and its spatial average Q
      const double wMin = Reduce::min(wField);
      VecOp::subVS(cField_, wField, wMin);
      VecOp::expVc(cField_, cField_, -1.0*size);
      double Q = Reduce::sum(cField_);
      Q = Q / double(nx);     // spatial average
      Q /= phiTot;            // correct for partial occupation

      // Note: phiTot = 1.0 except in the case of a mask that confines
      // material to a fraction of the unit cell.

      // Set ln q = ln Q - size*wMin and compute mu or phi
      SpeciesT::setLnQ(std::log(Q) - size*wMin);

      // Normalize concentration (the shift cancels in phi/Q)
      double prefactor = SpeciesT::phi()/Q;
      VecOp::mulEqS(cField_, prefactor);

   }

}
}
#endif

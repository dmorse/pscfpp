/*
* PSCF - Molecule Self-Consistent Field Theory
*
* Copyright 2016, The Regents of the University of Minnesota
* Distributed under the terms of the GNU General Public License.
*/

#include "Mixture.h"
#include <pscf/inter/Interaction.h>

namespace Pscf {
namespace Homogeneous {

   using namespace Util;

   /*
   * Constructor.
   */
   Mixture::Mixture()
    : ParamComposite(),
      molecules_(),
      mu_(),
      phi_(),
      c_(),
      w_(),
      nMolecule_(0),
      nMonomer_(0)
   {  setClassName("Mixture"); }

   /*
   * Destructor.
   */
   Mixture::~Mixture()
   {}

   /*
   * Read all parameters and initialize.
   */
   void Mixture::readParameters(std::istream& in)
   {
      // Precondition
      UTIL_ASSERT(molecules_.capacity() == 0);

      read<int>(in, "nMonomer", nMonomer_);
      c_.allocate(nMonomer_);
      w_.allocate(nMonomer_);

      read<int>(in, "nMolecule", nMolecule_);
      molecules_.allocate(nMolecule_);
      mu_.allocate(nMolecule_);
      phi_.allocate(nMolecule_);
      for (int i = 0; i < nMolecule_; ++i) {
         readParamComposite(in, molecules_[i]);
      }

      validate();
   }

   void Mixture::setNMolecule(int nMolecule)
   {
      UTIL_ASSERT(molecules_.capacity() == 0);
      nMolecule_ = nMolecule;
      molecules_.allocate(nMolecule_);
      mu_.allocate(nMolecule_);
      phi_.allocate(nMolecule_);
   }

   void Mixture::setNMonomer(int nMonomer)
   {
      UTIL_ASSERT(nMonomer_ == 0);
      nMonomer_ = nMonomer;
      c_.allocate(nMonomer_);
      w_.allocate(nMonomer_);
   }

   /*
   * Set molecular and monomer volume fractions.
   */
   void Mixture::setComposition(DArray<double> const & phi)
   {
      validate();
      UTIL_ASSERT(phi.capacity() == nMolecule_);

      // Set molecular volume fractions
      int i, j, k;
      for (i = 0; i < nMolecule_; ++i) {
         phi_[i] = phi[i];
      }

      // Initialize monomer fractions to zero
      for (k = 0; k < nMonomer_; ++k) {
         c_[k] = 0.0;
      }

      // Compute monomer volume fractions
      double concentration;
      for (i = 0; i < nMolecule_; ++i) {
         Molecule& mol = molecules_[i];
         concentration = phi_[i]/mol.size();
         for (j = 0; j < molecules_[i].nClump(); ++j) {
            k = mol.clump(j).monomerId();
            c_[k] += concentration*mol.clump(j).size();
         }
      }
      
   }

   /*
   * Set composition and compute thermodynamic properties.
   */
   void Mixture::computeMu(Interaction const & interaction,
                           DArray<double> const & phi, double xi)
   {
      UTIL_ASSERT(interaction.nMonomer() == nMonomer_);

      setComposition(phi);
      interaction.computeW(c_, w_);

      int imol, iclump, imon;
      double mu, size;
      for (imol = 0; imol < nMolecule_; ++imol) {
         Molecule& mol = molecules_[imol];
         mu = log( phi_[imol] );
         mu += xi*mol.size();
         for (iclump = 0; iclump < mol.nClump(); ++iclump) {
            imon = mol.clump(iclump).monomerId();
            size = mol.clump(iclump).size();
            mu += size*w_[imon];
         }
         mu_[imol] = mu;
      }
   }

   /*
   * Compute composition from chemical potentials.
   */
   void Mixture::computePhi(Interaction const & interaction, 
                           DArray<double> const & mu, 
                           DArray<double> const & phi, 
                           double& xi)
   {
      UTIL_ASSERT(interaction.nMonomer() == nMonomer_);

      // Allocate residual and jacobian on first use.
      if (residual_.capacity() == 0) {
         residual_.allocate(nMonomer_);
         jacobian_.allocate(nMonomer_, nMonomer_);
      }

      // Compute initial state, with assumed xi = 0
      xi = 0.0;
      computeMu(interaction, phi, xi);

      double error;
      double epsilon = 1.0E-6;
      computeResidual(mu, error);
      if (error < epsilon) return;

      //for (int it = 0, converged = false; it < 50 && !converged; ++it) {
         // Compute jacobian
         // Update
         // Set composition
         // computeMu
         // Test error (return if converged?)
      //}
   }

   void Mixture::computeResidual(DArray<double> const & mu, double& error)
   {
      error = 0.0;
      for (int i = 0; i < nMonomer_; ++i) {
         residual_[i] = mu_[i] - mu[i];
         if (residual_[i] > error) {
            error = residual_[i];
         }
      }
   }

   /*
   * Check validity after completing initialization.
   */
   void Mixture::validate() const
   {
      UTIL_ASSERT(nMolecule_ > 0);
      UTIL_ASSERT(nMonomer_ > 0);
      for (int i = 0; i < nMolecule_; ++i) {
         Molecule const & mol = molecules_[i];
         UTIL_ASSERT(mol.nClump() > 0);
         for (int j = 0; j < mol.nClump(); ++j) {
            UTIL_ASSERT(mol.clump(j).monomerId() < nMonomer_);
         }
      }
   }

} // namespace Homogeneous
} // namespace Pscf

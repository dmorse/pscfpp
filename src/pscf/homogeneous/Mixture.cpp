/*
* PSCF - Molecule Self-Consistent Field Theory
*
* Copyright 2016 - 2022, The Regents of the University of Minnesota
* Distributed under the terms of the GNU General Public License.
*/

#include "Mixture.h"
#include <pscf/inter/Interaction.h>
#include <pscf/math/LuSolver.h>
#include <cmath>

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
      residual_(),
      dX_(),
      dWdC_(),
      dWdPhi_(),
      jacobian_(),
      phiOld_(),
      fHelmholtz_(0.0),
      pressure_(0.0),
      solverPtr_(0),
      nMolecule_(0),
      nMonomer_(0),
      hasComposition_(false)
   {  setClassName("Mixture"); }

   /*
   * Destructor.
   */
   Mixture::~Mixture()
   {
      if (solverPtr_) {
         delete solverPtr_;
      }
   }

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
      double sum = 0;
      for (int i = 0; i < nMolecule_ - 1; ++i) {
         UTIL_CHECK(phi[i] >= 0.0);
         UTIL_CHECK(phi[i] <= 1.0);
         phi_[i] = phi[i];
         sum += phi[i];
      }
      UTIL_CHECK(sum <= 1.0);
      phi_[nMolecule_ -1] = 1.0 - sum;

      computeC();
      hasComposition_ = true;
   }

   void Mixture::computeC()
   {
      // Initialize monomer fractions to zero
      int k;
      for (k = 0; k < nMonomer_; ++k) {
         c_[k] = 0.0;
      }

      // Compute monomer volume fractions
      double concentration;
      int i, j;
      for (i = 0; i < nMolecule_; ++i) {
         Molecule& mol = molecules_[i];
         concentration = phi_[i]/mol.size();
         for (j = 0; j < molecules_[i].nClump(); ++j) {
            k = mol.clump(j).monomerId();
            c_[k] += concentration*mol.clump(j).size();
         }
      }

      // Check and normalize monomer volume fractions
      double sum = 0.0;
      for (k = 0; k < nMonomer_; ++k) {
         UTIL_CHECK(c_[k] >= 0.0);
         UTIL_CHECK(c_[k] <= 1.0);
         sum += c_[k];
      }
      UTIL_CHECK(sum > 0.9999);
      UTIL_CHECK(sum < 1.0001);
      for (k = 0; k < nMonomer_; ++k) {
         c_[k] /= sum;
      }

   }

   /*
   * Compute thermodynamic properties.
   */
   void Mixture::computeMu(Interaction const & interaction, 
                           double xi)
   {
      UTIL_CHECK(interaction.nMonomer() == nMonomer_);
      UTIL_CHECK(hasComposition_);

      // Compute monomer excess chemical potentials
      interaction.computeW(c_, w_);

      int m; // molecule index
      int c; // clump index
      int t; // monomer type index
      double mu, size;
      for (m = 0; m < nMolecule_; ++m) {
         Molecule& mol = molecules_[m];
         mu = log( phi_[m] );
         mu += xi*mol.size();
         for (c = 0; c < mol.nClump(); ++c) {
            t = mol.clump(c).monomerId();
            size = mol.clump(c).size();
            mu += size*w_[t];
         }
         mu_[m] = mu;
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
         dX_.allocate(nMolecule_);
         dWdC_.allocate(nMonomer_, nMonomer_);
         dWdPhi_.allocate(nMonomer_, nMolecule_);
         jacobian_.allocate(nMolecule_, nMolecule_);
         phiOld_.allocate(nMolecule_);
         solverPtr_ = new LuSolver();
         solverPtr_->allocate(nMolecule_);
      }

      // Compute initial state
      setComposition(phi);
      computeMu(interaction, xi);
      adjustXi(mu, xi);

      // Compute initial residual
      double error;
      double epsilon = 1.0E-10;
      computeResidual(mu, error);

      #if 0
      std::cout << "mu[0] =" << mu[0] << std::endl;
      std::cout << "mu[1] =" << mu[1] << std::endl;
      std::cout << "mu_[0] =" << mu_[0] << std::endl;
      std::cout << "mu_[1] =" << mu_[1] << std::endl;
      std::cout << "residual[0] =" << residual_[0] << std::endl;
      std::cout << "residual[1] =" << residual_[1] << std::endl;
      std::cout << "error     ="  << error << std::endl;
      #endif

      if (error < epsilon) return;

      double s1;  // clump size
      double v1;  // molecule size
      double f1;  // clump fraction
      int m1, m2; // molecule type indices
      int c1;     // clump index
      int t1, t2; // monomer type indices

      for (int it = 0; it < 50; ++it) {

         // Compute matrix of derivative dWdC (C = monomer fraction)
         interaction.computeDwDc(c_, dWdC_);

         // Compute matrix derivative dWdPhi (Phi = molecule fraction)
         for (t1 = 0; t1 < nMonomer_; ++t1) {
            for (m1 = 0; m1 < nMolecule_; ++m1) {
               dWdPhi_(t1, m1) = 0.0;
            }
         }
         for (m1 = 0; m1 < nMolecule_; ++m1) {
            v1 = molecule(m1).size();
            for (c1 = 0; c1 < molecule(m1).nClump(); ++c1) {
               t1 = molecule(m1).clump(c1).monomerId();
               s1 = molecule(m1).clump(c1).size();
               f1 = s1/v1;
               for (t2 = 0; t2 < nMonomer_; ++t2) {
                  dWdPhi_(t2, m1) += dWdC_(t2, t1)*f1;
               }
            }
         }

         // Compute matrix d(mu)/d(Phi), stored in jacobian_
         for (m1 = 0; m1 < nMolecule_; ++m1) {
            for (m2 = 0; m2 < nMolecule_; ++m2) {
               jacobian_(m1, m2) = 0.0;
            }
            jacobian_(m1, m1) = 1.0/phi_[m1];
         }
         for (m1 = 0; m1 < nMolecule_; ++m1) {
            v1 = molecule(m1).size();
            for (c1 = 0; c1 < molecule(m1).nClump(); ++c1) {
               t1 = molecule(m1).clump(c1).monomerId();
               s1 = molecule(m1).clump(c1).size();
               for (m2 = 0; m2 < nMolecule_; ++m2) {
                  jacobian_(m1, m2) += s1*dWdPhi_(t1, m2);
               }
            }
         }

         // Impose incompressibility
         int mLast = nMolecule_ - 1;
         for (m1 = 0; m1 < nMolecule_; ++m1) {
            for (m2 = 0; m2 < nMolecule_ - 1; ++m2) {
               jacobian_(m1, m2) -= jacobian_(m1, mLast);
            }
            // Derivative of mu_[m1] with respect to xi
            jacobian_(m1, mLast) = molecule(m1).size();
         }

         // Newton Raphson update of phi and xi fields
         solverPtr_->computeLU(jacobian_);
         solverPtr_->solve(residual_, dX_);

         // Store old value of phi
         for (m1 = 0; m1 < nMolecule_; ++m1) {
            phiOld_[m1] = phi_[m1];
         }

         // Apply update
         bool inRange = false;
         for (int j = 0; j < 5; ++j) {

            // Update volume fractions
            double sum = 0.0;
            for (m1 = 0; m1 < nMolecule_ - 1; ++m1) {
               phi_[m1] = phiOld_[m1] - dX_[m1];
               sum += phi_[m1];
            }
            phi_[mLast] = 1.0 - sum;

            // Check if all volume fractions are in [0,1]
            inRange = true;
            for (m1 = 0; m1 < nMolecule_; ++m1) {
               if (phi_[m1] < 0.0) inRange = false;
               if (phi_[m1] > 1.0) inRange = false;
            }

            // Exit loop or reduce increment
            if (inRange) {
               break;
            } else {
               for (m1 = 0; m1 < nMolecule_; ++m1) {
                  dX_[m1] *= 0.5;
               }
            }

         }
         if (inRange) {
            xi = xi - dX_[mLast];
         } else {
            UTIL_THROW("Volume fractions remain out of range");
         }

         // Compute residual
         computeC();
         computeMu(interaction, xi);
         adjustXi(mu, xi);
         computeResidual(mu, error);

         #if 0
         std::cout << "mu[0] =" << mu[0] << std::endl;
         std::cout << "mu[1] =" << mu[1] << std::endl;
         std::cout << "mu_[0] =" << mu_[0] << std::endl;
         std::cout << "mu_[1] =" << mu_[1] << std::endl;
         std::cout << "residual[0] =" << residual_[0] << std::endl;
         std::cout << "residual[1] =" << residual_[1] << std::endl;
         std::cout << "error     ="  << error << std::endl;
         #endif

         if (error < epsilon) return;
      }

      UTIL_THROW("Failed to converge");
   }

   void Mixture::adjustXi(DArray<double> const & mu, double& xi)
   {
      double dxi = 0.0;
      double sum = 0.0;
      for (int i=0; i < nMolecule_; ++i) {
         dxi += mu[i] - mu_[i];
         sum += molecules_[i].size();
      }
      dxi = dxi/sum;
      for (int i=0; i < nMolecule_; ++i) {
         mu_[i] += dxi*molecules_[i].size();
      }
      xi += dxi;
   }

   void Mixture::computeResidual(DArray<double> const & mu, double& error)
   {
      error = 0.0;
      for (int i = 0; i < nMonomer_; ++i) {
         residual_[i] = mu_[i] - mu[i];
         if (std::abs(residual_[i]) > error) {
            error = std::abs(residual_[i]);
         }
      }
   }

   /*
   * Compute Helmoltz free energy and pressure
   */
   void Mixture::computeFreeEnergy(Interaction const & interaction)
   {
      fHelmholtz_ = 0.0;
 
      // Compute ideal gas contributions to fHelhmoltz_
      double size;
      for (int i = 0; i < nMolecule_; ++i) {
         size = molecules_[i].size();
         fHelmholtz_ += phi_[i]*( log(phi_[i]) - 1.0 )/size;
      }

      // Add average interaction free energy density per monomer
      fHelmholtz_ += interaction.fHelmholtz(c_);

      // Compute pressure
      pressure_ = -fHelmholtz_;
      for (int i = 0; i < nMolecule_; ++i) {
         size = molecules_[i].size();
         pressure_ += phi_[i]*mu_[i]/size;
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

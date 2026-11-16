/*
* PSCF - Molecule Self-Consistent Field Theory
*
* Copyright 2015 - 2026, David Morse (morse012@um.edu) and co-authors
* Distributed under the terms of the GNU General Public License.
*/

#include "FhMixture.tpp"
#include <pscf/floryHuggins/FhInteraction.h>
#include <pscf/math/LuSolver.h>
#include <cmath>

namespace Pscf {

   using namespace Util;

   namespace {

      /*
      * Is a species with chemical potential mu absent (mu = -infinity)?
      */
      bool isAbsent(double mu)
      {  return std::isinf(mu) && mu < 0.0; }

   }

   /*
   * Constructor.
   */
   FhMixture::FhMixture()
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
      lnPhi_(),
      lnPhiOld_(),
      fHelmholtz_(0.0),
      pressure_(0.0),
      solverPtr_(0),
      nMolecule_(0),
      nMonomer_(0),
      hasComposition_(false)
   {  setClassName("FhMixture"); }

   /*
   * Destructor.
   */
   FhMixture::~FhMixture()
   {
      if (solverPtr_) {
         delete solverPtr_;
      }
   }

   /*
   * Read all parameters and initialize.
   */
   void FhMixture::readParameters(std::istream& in)
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
      lnPhi_.allocate(nMolecule_);
      for (int i = 0; i < nMolecule_; ++i) {
         readParamComposite(in, molecules_[i]);
      }

      validate();
   }

   void FhMixture::setNMolecule(int nMolecule)
   {
      UTIL_CHECK(molecules_.capacity() == 0);
      UTIL_CHECK(nMolecule > 0);
      nMolecule_ = nMolecule;
      molecules_.allocate(nMolecule_);
      mu_.allocate(nMolecule_);
      phi_.allocate(nMolecule_);
      lnPhi_.allocate(nMolecule_);
   }

   void FhMixture::setNMonomer(int nMonomer)
   {
      UTIL_CHECK(nMonomer_ == 0);
      UTIL_CHECK(nMonomer > 0);
      nMonomer_ = nMonomer;
      c_.allocate(nMonomer_);
      w_.allocate(nMonomer_);
   }
    
   /*
   * Set molecular and monomer volume fractions.
   */
   void FhMixture::setComposition(DArray<double> const & phi)
   {
      validate();
      UTIL_ASSERT(phi.capacity() == nMolecule_);

      // Set molecular volume fractions. The phi of the species with
      // the largest phi is 1 minus the sum of the others: that of any
      // other species would lose its precision if it is dilute.
      int iMax = 0;
      for (int i = 1; i < nMolecule_; ++i) {
         if (phi[i] > phi[iMax]) iMax = i;
      }
      double sum = 0;
      for (int i = 0; i < nMolecule_; ++i) {
         if (i == iMax) continue;
         UTIL_CHECK(phi[i] >= 0.0);
         UTIL_CHECK(phi[i] <= 1.0);
         phi_[i] = phi[i];
         sum += phi[i];
      }
      UTIL_CHECK(sum <= 1.0);
      phi_[iMax] = 1.0 - sum;
      for (int i = 0; i < nMolecule_; ++i) {
         lnPhi_[i] = log(phi_[i]);
      }

      computeC();
      hasComposition_ = true;
   }

   void FhMixture::computeC()
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
         FhMolecule& mol = molecules_[i];
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
   void FhMixture::computeMu(FhInteraction const & interaction, 
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
         FhMolecule& mol = molecules_[m];
         mu = lnPhi_[m];
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
   void FhMixture::computePhi(FhInteraction const & interaction, 
                              DArray<double> const & mu, 
                           DArray<double> const & phi, 
                           double& xi)
   {
      UTIL_ASSERT(interaction.nMonomer() == nMonomer_);

      // Allocate residual and jacobian on first use.
      if (residual_.capacity() == 0) {
         residual_.allocate(nMolecule_);
         dX_.allocate(nMolecule_);
         dWdC_.allocate(nMonomer_, nMonomer_);
         dWdPhi_.allocate(nMonomer_, nMolecule_);
         jacobian_.allocate(nMolecule_, nMolecule_);
         lnPhiOld_.allocate(nMolecule_);
         solverPtr_ = new LuSolver();
         solverPtr_->allocate(nMolecule_);
      }

      // Initial composition: the guess (except for the dependent
      // species, see below). A species with mu = -infinity (a closed
      // species with phi = 0) is absent: its phi stays 0.
      validate();
      int m1, m2; // molecule type indices
      for (m1 = 0; m1 < nMolecule_; ++m1) {
         if (isAbsent(mu[m1])) {
            phi_[m1] = 0.0;
         } else {
            UTIL_CHECK(phi[m1] > 0.0);
            UTIL_CHECK(phi[m1] <= 1.0);
            phi_[m1] = phi[m1];
         }
         lnPhi_[m1] = log(phi_[m1]);
      }
      hasComposition_ = true;

      double error;
      double epsilon = 1.0E-10;
      double s1;  // clump size
      double v1;  // molecule size
      double f1;  // clump fraction
      int c1;     // clump index
      int t1, t2; // monomer type indices
      int mDep = -1; // index of the dependent species
      double sum;

      for (int it = 0; it < 50; ++it) {

         // The dependent species is the species with the largest phi:
         // its phi is 1 minus the sum of the others, and xi replaces it
         // as a Newton variable. Any other choice loses the relative
         // precision of the phi of the dependent species when it is
         // dilute, so that the residual cannot reach the tolerance. It
         // is chosen again in each iteration, since the largest species
         // can change. The residual is computed for the initial state
         // and after a change of the dependent species.
         int mMax = 0;
         for (m1 = 1; m1 < nMolecule_; ++m1) {
            if (phi_[m1] > phi_[mMax]) mMax = m1;
         }
         if (isAbsent(mu[mMax])) {
            UTIL_THROW("No species is present (all mu = -infinity)");
         }
         if (mMax != mDep) {
            mDep = mMax;
            sum = 0.0;
            for (m1 = 0; m1 < nMolecule_; ++m1) {
               if (m1 != mDep) sum += phi_[m1];
            }
            UTIL_CHECK(sum < 1.0);
            phi_[mDep] = 1.0 - sum;
            lnPhi_[mDep] = std::log1p(-sum);
            computeC();
            computeMu(interaction, xi);
            adjustXi(mu, xi);
            computeResidual(mu, error);
            if (error < epsilon) return;
         }

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

         // Compute the interaction part of matrix d(mu)/d(Phi), stored
         // in jacobian_ (the ideal part is added below)
         for (m1 = 0; m1 < nMolecule_; ++m1) {
            for (m2 = 0; m2 < nMolecule_; ++m2) {
               jacobian_(m1, m2) = 0.0;
            }
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
         for (m1 = 0; m1 < nMolecule_; ++m1) {
            for (m2 = 0; m2 < nMolecule_; ++m2) {
               if (m2 != mDep) {
                  jacobian_(m1, m2) -= jacobian_(m1, mDep);
               }
            }
            // Derivative of mu_[m1] with respect to xi
            jacobian_(m1, mDep) = molecule(m1).size();
         }

         // The Newton variables are xi and ln(phi) of every species
         // other than the dependent one, mDep, whose phi is 1 minus
         // the sum of the others. Steps in ln(phi) reach a dilute
         // species, whose mu depends on ln(phi), in a few steps, and
         // ln(phi) remains finite where phi underflows to zero. The
         // columns of these variables are d(mu)/d(ln phi) =
         // phi*d(mu)/d(phi), plus the ideal part:
         // d(ln phi[m2])/d(ln phi[m2]) = 1, and
         // d(ln phi[mDep])/d(ln phi[m2]) = -phi[m2]/phi[mDep].
         for (m2 = 0; m2 < nMolecule_; ++m2) {
            if (m2 == mDep) continue;
            for (m1 = 0; m1 < nMolecule_; ++m1) {
               jacobian_(m1, m2) *= phi_[m2];
            }
            jacobian_(m2, m2) += 1.0;
            jacobian_(mDep, m2) -= phi_[m2]/phi_[mDep];
         }

         // An absent species has a unit row and a zero residual, so
         // that its ln(phi) = -infinity does not change
         for (m1 = 0; m1 < nMolecule_; ++m1) {
            if (isAbsent(mu[m1])) {
               for (m2 = 0; m2 < nMolecule_; ++m2) {
                  jacobian_(m1, m2) = 0.0;
               }
               jacobian_(m1, m1) = 1.0;
            }
         }

         // Newton Raphson update of phi and xi fields. The Jacobian is
         // singular if the composition is at a spinodal of the mixture,
         // e.g., at the critical point of a symmetric blend.
         solverPtr_->computeLU(jacobian_);
         if (solverPtr_->isSingular()) {
            UTIL_THROW("Singular Jacobian (composition at a spinodal?)");
         }
         solverPtr_->solve(residual_, dX_);

         // Store old values of ln(phi)
         for (m1 = 0; m1 < nMolecule_; ++m1) {
            lnPhiOld_[m1] = lnPhi_[m1];
         }

         // Apply update
         bool inRange = false;
         for (int j = 0; j < 30; ++j) {

            // Update volume fractions
            sum = 0.0;
            for (m1 = 0; m1 < nMolecule_; ++m1) {
               if (m1 == mDep) continue;
               lnPhi_[m1] = lnPhiOld_[m1] - dX_[m1];
               phi_[m1] = std::exp(lnPhi_[m1]);
               sum += phi_[m1];
            }
            phi_[mDep] = 1.0 - sum;

            // Check if all volume fractions are in (0,1]
            inRange = (phi_[mDep] > 0.0);
            for (m1 = 0; m1 < nMolecule_; ++m1) {
               if (m1 == mDep) continue;
               if (!(lnPhi_[m1] <= 0.0)) inRange = false;
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
            lnPhi_[mDep] = std::log1p(-sum);
            xi = xi - dX_[mDep];
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

   void FhMixture::adjustXi(DArray<double> const & mu, double& xi)
   {
      double dxi = 0.0;
      double sum = 0.0;
      for (int i=0; i < nMolecule_; ++i) {
         if (isAbsent(mu[i])) continue;
         dxi += mu[i] - mu_[i];
         sum += molecules_[i].size();
      }
      dxi = dxi/sum;
      for (int i=0; i < nMolecule_; ++i) {
         mu_[i] += dxi*molecules_[i].size();
      }
      xi += dxi;
   }

   void FhMixture::computeResidual(DArray<double> const & mu, double& error)
   {
      error = 0.0;
      for (int i = 0; i < nMolecule_; ++i) {
         if (isAbsent(mu[i])) {
            residual_[i] = 0.0;
            continue;
         }
         residual_[i] = mu_[i] - mu[i];
         if (std::abs(residual_[i]) > error) {
            error = std::abs(residual_[i]);
         }
      }
   }

   /*
   * Compute Helmoltz free energy and pressure
   */
   void FhMixture::computeFreeEnergy(FhInteraction const & interaction)
   {
      fHelmholtz_ = 0.0;
 
      // Compute ideal gas contributions to fHelhmoltz_
      double size;
      for (int i = 0; i < nMolecule_; ++i) {
         size = molecules_[i].size();
         // phi*(log(phi) - 1) tends to 0 as phi tends to 0
         if (phi_[i] > 0.0) {
            fHelmholtz_ += phi_[i]*( log(phi_[i]) - 1.0 )/size;
         }
      }

      // Add average interaction free energy density per monomer
      fHelmholtz_ += interaction.fHelmholtz(c_);

      // Compute pressure
      pressure_ = -fHelmholtz_;
      for (int i = 0; i < nMolecule_; ++i) {
         size = molecules_[i].size();
         // phi*mu tends to 0 as phi tends to 0 (mu ~ log(phi))
         if (phi_[i] > 0.0) {
            pressure_ += phi_[i]*mu_[i]/size;
         }
      }

   }

   /*
   * Check validity after completing initialization.
   */
   void FhMixture::validate() const
   {
      UTIL_ASSERT(nMolecule_ > 0);
      UTIL_ASSERT(nMonomer_ > 0);
      for (int i = 0; i < nMolecule_; ++i) {
         FhMolecule const & mol = molecules_[i];
         UTIL_ASSERT(mol.nClump() > 0);
         for (int j = 0; j < mol.nClump(); ++j) {
            UTIL_ASSERT(mol.clump(j).monomerId() < nMonomer_);
         }
      }
   }

   // Explicit instantiation definition for member function template.
   template 
   void FhMixture::initialize<double>(Composition<double> const & );

} // namespace Pscf

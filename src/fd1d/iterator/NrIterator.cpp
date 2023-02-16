/*
* PSCF - Polymer Self-Consistent Field Theory
*
* Copyright 2016 - 2022, The Regents of the University of Minnesota
* Distributed under the terms of the GNU General Public License.
*/

#include "NrIterator.h"
#include <fd1d/System.h>
#include <pscf/inter/Interaction.h>
#include <util/misc/Log.h>

#include <math.h>

namespace Pscf {
namespace Fd1d
{

   using namespace Util;

   NrIterator::NrIterator()
    : Iterator(),
      epsilon_(0.0),
      maxItr_(0),
      isAllocated_(false),
      newJacobian_(false),
      needsJacobian_(true),
      isCanonical_(true)
   {  setClassName("NrIterator"); }

   NrIterator::NrIterator(System& system)
    : Iterator(system),
      epsilon_(0.0),
      maxItr_(0),
      isAllocated_(false),
      newJacobian_(false),
      needsJacobian_(true),
      isCanonical_(true)
   {  setClassName("NrIterator"); }

   NrIterator::~NrIterator()
   {}

   void NrIterator::readParameters(std::istream& in)
   {
      maxItr_ = 400;
      read(in, "epsilon", epsilon_);
      readOptional(in, "maxItr", maxItr_);
      setup();
   }

   void NrIterator::setup()
   {
      int nm = system().mixture().nMonomer();   // number of monomer types
      int nx = domain().nx(); // number of grid points
      UTIL_CHECK(nm > 0);
      UTIL_CHECK(nx > 0);
      int nr = nm*nx;                  // number of residual components
      if (isAllocated_) {
         UTIL_CHECK(cArray_.capacity() == nm);
         UTIL_CHECK(residual_.capacity() == nr);
      } else {
         cArray_.allocate(nm);
         wArray_.allocate(nm);
         residual_.allocate(nr);
         jacobian_.allocate(nr, nr);
         residualNew_.allocate(nr);
         dOmega_.allocate(nr);
         wFieldsNew_.allocate(nm);
         cFieldsNew_.allocate(nm);
         for (int i = 0; i < nm; ++i) {
            wFieldsNew_[i].allocate(nx);
            cFieldsNew_[i].allocate(nx);
         }
         solver_.allocate(nr);
         isAllocated_ = true;
      }
   }

   void NrIterator::computeResidual(Array<WField> const & wFields,
                                    Array<CField> const & cFields,
                                    Array<double>& residual)
   {
      int nm = system().mixture().nMonomer();  // number of monomer types
      int nx = domain().nx();         // number of grid points
      int i;                          // grid point index
      int j;                          // monomer indices
      int ir;                         // residual index

      // Loop over grid points
      for (i = 0; i < nx; ++i) {

         // Copy volume fractions at grid point i to cArray_
         for (j = 0; j < nm; ++j) {
            cArray_[j] = cFields[j][i];
         }

         // Compute w fields, without Langrange multiplier, from c fields
         system().interaction().computeW(cArray_, wArray_);

         // Initial residual = wPredicted(from above) - actual w
         for (j = 0; j < nm; ++j) {
            ir = j*nx + i;
            residual[ir] = wArray_[j] - wFields[j][i];
         }

         // Residuals j = 1, ..., nm-1 are differences from component j=0
         for (j = 1; j < nm; ++j) {
            ir = j*nx + i;
            residual[ir] = residual[ir] - residual[i];
         }

         // Residual for component j=0 then imposes incompressiblity
         residual[i] = -1.0;
         for (j = 0; j < nm; ++j) {
            residual[i] += cArray_[j];
         }
      }

      /*
      * Note: In canonical ensemble, the spatial integral of the
      * incompressiblity residual is guaranteed to be zero, as a result of how
      * volume fractions are computed in SCFT. One of the nx incompressibility
      * constraints is thus redundant. To avoid this redundancy, replace the
      * incompressibility residual at the last grid point by a residual that
      * requires the w field for the last monomer type at the last grid point
      * to equal zero.
      */

      if (isCanonical_) {
         residual[nx-1] = wFields[nm-1][nx-1];
      }

   }

   /*
   * Compute Jacobian matrix numerically, by evaluating finite differences.
   */
   void NrIterator::computeJacobian()
   {
      int nm = system().mixture().nMonomer();   // number of monomer types
      int nx = domain().nx();          // number of grid points
      int i;                           // monomer index
      int j;                           // grid point index

      // Copy system().wFields to wFieldsNew.
      for (i = 0; i < nm; ++i) {
         for (j = 0; j < nx; ++j) {
            UTIL_CHECK(nx == wFieldsNew_[i].capacity());
            UTIL_CHECK(nx == system().wField(i).capacity());
            wFieldsNew_[i][j] = system().wField(i)[j];
         }
      }

      // Compute jacobian, column by column
      double delta = 0.001;
      int nr = nm*nx;              // number of residual elements
      int jr;                      // jacobian row index
      int jc = 0;                  // jacobian column index
      for (i = 0; i < nm; ++i) {
         for (j = 0; j < nx; ++j) {
            wFieldsNew_[i][j] += delta;
            system().mixture().compute(wFieldsNew_, cFieldsNew_);
            computeResidual(wFieldsNew_, cFieldsNew_, residualNew_);
            for (jr = 0; jr < nr; ++jr) {
               jacobian_(jr, jc) =
                    (residualNew_[jr] - residual_[jr])/delta;
            }
            wFieldsNew_[i][j] = system().wField(i)[j];
            ++jc;
         }
      }

      // Decompose Jacobian matrix
      solver_.computeLU(jacobian_);

   }

   void NrIterator::incrementWFields(Array<WField> const & wOld,
                                     Array<double> const & dW,
                                     Array<WField> & wNew)
   {
      int nm = system().mixture().nMonomer(); // number of monomers types
      int nx = domain().nx();        // number of grid points
      int i;                         // monomer index
      int j;                         // grid point index
      int k = 0;                     // residual element index

      // Add dW
      for (i = 0; i < nm; ++i) {
         for (j = 0; j < nx; ++j) {
            wNew[i][j] = wOld[i][j] - dW[k];
            ++k;
         }
      }

      // If canonical, shift such that last element is exactly zero
      if (isCanonical_) {
         double shift = wNew[nm-1][nx-1];
         for (i = 0; i < nm; ++i) {
            for (j = 0; j < nx; ++j) {
               wNew[i][j] -= shift;
            }
         }
      }

   }

   double NrIterator::residualNorm(Array<double> const & residual) const
   {
      int nm = system().mixture().nMonomer();  // number of monomer types
      int nx = domain().nx();         // number of grid points
      int nr = nm*nx;                 // number of residual components
      double value, norm;
      norm = 0.0;
      for (int ir = 0; ir <  nr; ++ir) {
         value = fabs(residual[ir]);
         if (value > norm) {
            norm = value;
         }
      }
      return norm;
   }

   int NrIterator::solve(bool isContinuation)
   {
      int nm = system().mixture().nMonomer();  // number of monomer types
      int np = system().mixture().nPolymer();  // number of polymer species
      int nx = domain().nx();         // number of grid points
      int nr = nm*nx;                 // number of residual elements

      // Determine if isCanonical (iff all species ensembles are closed)
      isCanonical_ = true;
      Species::Ensemble ensemble;
      for (int i = 0; i < np; ++i) {
         ensemble = system().mixture().polymer(i).ensemble();
         if (ensemble == Species::Unknown) {
            UTIL_THROW("Unknown species ensemble");
         }
         if (ensemble == Species::Open) {
            isCanonical_ = false;
         }
      }

      // If isCanonical, shift so that last element is zero.
      // Note: This is one of the residuals in this case.
      if (isCanonical_) {
         double shift = wFields()[nm-1][nx-1];
         int i, j;
         for (i = 0; i < nm; ++i) {
            for (j = 0; j < nx; ++j) {
               wFields()[i][j] -= shift;
            }
         }
      }

      // Compute initial residual vector and norm
      system().mixture().compute(system().wFields(), system().cFields());
      computeResidual(system().wFields(), system().cFields(), residual_);
      double norm = residualNorm(residual_);

      // Set Jacobian status
      newJacobian_ = false;
      if (!isContinuation) {
         needsJacobian_ = true;
      }

      // Iterative loop
      double normNew;
      int i, j, k;
      for (i = 0; i < maxItr_; ++i) {
         Log::file() << "iteration " << i
                     << " , error = " << norm
                     << std::endl;

         if (norm < epsilon_) {
            Log::file() << "Converged" << std::endl;
            system().computeFreeEnergy();
            // Success
            return 0;
         }

         if (needsJacobian_) {
            Log::file() << "Computing jacobian" << std::endl;;
            computeJacobian();
            newJacobian_ = true;
            needsJacobian_ = false;
         }

         // Compute Newton-Raphson increment dOmega_
         solver_.solve(residual_, dOmega_);

         // Try full Newton-Raphson update
         incrementWFields(system().wFields(), dOmega_, wFieldsNew_);
         system().mixture().compute(wFieldsNew_, cFieldsNew_);
         computeResidual(wFieldsNew_, cFieldsNew_, residualNew_);
         normNew = residualNorm(residualNew_);

         // Decrease increment if necessary
         j = 0;
         while (normNew > norm && j < 3) {
            Log::file() << "      decreasing increment,  error = "
                        << normNew << std::endl;
            needsJacobian_ = true;
            for (k = 0; k < nr; ++k) {
               dOmega_[k] *= 0.66666666;
            }
            incrementWFields(system().wFields(), dOmega_, wFieldsNew_);
            system().mixture().compute(wFieldsNew_, cFieldsNew_);
            computeResidual(wFieldsNew_, cFieldsNew_, residualNew_);
            normNew = residualNorm(residualNew_);
            ++j;
         }

         // If necessary, try reversing direction
         if (normNew > norm) {
            Log::file() << "      reversing increment,  norm = "
                        << normNew << std::endl;
            needsJacobian_ = true;
            for (k = 0; k < nr; ++k) {
               dOmega_[k] *= -1.000;
            }
            incrementWFields(system().wFields(), dOmega_, wFieldsNew_);
            system().mixture().compute(wFieldsNew_, cFieldsNew_);
            computeResidual(wFieldsNew_, cFieldsNew_, residualNew_);
            normNew = residualNorm(residualNew_);
         }

         // Accept or reject update
         if (normNew < norm) {

            // Update system fields and residual vector
            for (j = 0; j < nm; ++j) {
               for (k = 0; k < nx; ++k) {
                  system().wField(j)[k] = wFieldsNew_[j][k];
                  system().cField(j)[k] = cFieldsNew_[j][k];
               }
            }
            for (j = 0; j < nr; ++j) {
               residual_[j] = residualNew_[j];
            }
            newJacobian_ = false;
            if (!needsJacobian_) {
               if (normNew/norm > 0.5) {
                  needsJacobian_ = true;
               }
            }
            norm = normNew;
         } else {
            Log::file() << "Iteration failed, norm = "
                      << normNew << std::endl;
            if (newJacobian_) {
               return 1;
               Log::file() << "Unrecoverable failure " << std::endl;
            } else {
               Log::file() << "Try rebuilding Jacobian" << std::endl;
               needsJacobian_ = true;
            }
         }

      }

      // Failure
      return 1;
   }

} // namespace Fd1d
} // namespace Pscf

/*
* PSCF - Polymer Self-Consistent Field Theory
*
* Copyright 2016, The Regents of the University of Minnesota
* Distributed under the terms of the GNU General Public License.
*/

#include "NrIterator.h"
#include "System.h"
#include "Mixture.h"
#include <pscf/Interaction.h>

#include <math.h>

namespace Pscf {
namespace Fd1d
{

   using namespace Util;

   NrIterator::NrIterator()
    : Iterator(),
      epsilon_(0.0),
      isAllocated_(false)
   {  setClassName("NrIterator"); }

   NrIterator::~NrIterator()
   {}

   void NrIterator::readParameters(std::istream& in)
   {
      read(in, "epsilon", epsilon_);
      if (domain().nx() > 0) {
         allocate();
      }
   }

   void NrIterator::allocate()
   {
      int nm = mixture().nMonomer();   // number of monomer types
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
      int nm = mixture().nMonomer();  // number of monomer types
      int nx = domain().nx();         // number of grid points
      int i;                          // grid point index
      int j;                          // monomer indices
      int ir;                         // residual index
      for (i = 0; i < nx; ++i) {
         for (j = 0; j < nm; ++j) {
            cArray_[j] = cFields[j][i];
         }
         system().interaction().computeW(cArray_, 0, wArray_);
         for (j = 0; j < nm; ++j) {
            ir = j*nx + i;
            residual[ir] = wArray_[j] - wFields[j][i];
         }
         for (j = 1; j < nm; ++j) {
            ir = j*nx + i;
            residual[ir] = residual[ir] - residual[i];
         }
         residual[i] = -1.0;
         for (j = 0; j < nm; ++j) {
            residual[i] += cArray_[j];
         }
      }

      // Use last incompressibility residual to set last field point
      residual[nx-1] = wFields[nm-1][nx-1];
   }

   void NrIterator::computeJacobian()
   {
      // std::cout << "Begin computeJacobian ... ";
      int nm = mixture().nMonomer();   // number of monomer types
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
            mixture().compute(wFieldsNew_, cFieldsNew_);
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

      // std::cout << "Finish computeJacobian" << std::endl;
   }

   #if 0
   void NrIterator::update()
   {
      // std::cout << "Begin update .. ";
      computeJacobian();

      // Compute increment dOmega_
      solver_.solve(residual_, dOmega_);

      // Increment wFields;
      incrementWFields(system().wFields(), dOmega_);

      #if 0
      int nm = mixture().nMonomer();   // number of monomers types
      int nx = domain().nx();          // number of grid points
      int i;                           // monomer index
      int j;                           // grid point index
      int k = 0;                       // residual element index
      for (i = 0; i < nm; ++i) {
         for (j = 0; j < nx; ++j) {
            system().wField(i)[j] -= dOmega_[k];
            ++k;
         }
      }
      #endif

      mixture().compute(system().wFields(), system().cFields());
      computeResidual(system().wFields(), system().cFields(), residual_);

      //std::cout << "Finish update" << std::endl;
   }
   #endif

   void NrIterator::incrementWFields(Array<WField> const & wOld, 
                                     Array<double> const & dW, 
                                     Array<WField> & wNew)
   {
      int nm = mixture().nMonomer(); // number of monomers types
      int nx = domain().nx();        // number of grid points
      int i;                         // monomer index
      int j;                         // grid point index
      int k = 0;                     // residual element index
      for (i = 0; i < nm; ++i) {
         for (j = 0; j < nx; ++j) {
            wNew[i][j] = wOld[i][j] - dW[k];
            ++k;
         }
      }
   }
   
   double NrIterator::residualNorm(Array<double> const & residual) const
   {
      int nm = mixture().nMonomer();  // number of monomer types
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

   int NrIterator::solve()
   {
      // Allocate memory if needed or, if allocated, check array sizes.
      allocate();

      // Compute initial residual vector and norm
      mixture().compute(system().wFields(), system().cFields());
      system().computeFreeEnergy();
      computeResidual(system().wFields(), system().cFields(), residual_);
      double norm = residualNorm(residual_);
      double normNew;

      int nm = mixture().nMonomer();  // number of monomer types
      int nx = domain().nx();         // number of grid points
      int nr = nm*nx;                 // number of residual elements
      int i, j, k;
      bool hasJacobian = false;
      bool needsJacobian = true;

      // Iterative loop
      for (i = 0; i < 100; ++i) {
         std::cout << "Begin iteration " << i;

         #if 0
         std::cout << "\n";
         for (j = 0; j < nx; ++j) {
             for (k = 0; k < nm; ++k) {
                std::cout << system().wField(k)[j] << "  ";
             }
             for (k = 0; k < nm; ++k) {
                std::cout << system().cField(k)[j] << "  ";
             }
             for (k = 0; k < nm; ++k) {
                std::cout << residual_[k*nx + j] << "  ";
             }
             std::cout << "\n";
         }
         #endif

         std::cout << " , residual norm = " << norm << std::endl;

         if (norm < epsilon_) {
            std::cout << "Converged" << std::endl;
            return 0;
         } 

         if (needsJacobian) {
            std::cout << "      computing jacobian" << std::endl;;
            computeJacobian();
            hasJacobian = true;
            needsJacobian = false;
         }

         // Compute Newton-Raphson increment dOmega_
         solver_.solve(residual_, dOmega_);

         // Try full Newton-Raphson update
         incrementWFields(system().wFields(), dOmega_, wFieldsNew_);
         mixture().compute(wFieldsNew_, cFieldsNew_);
         computeResidual(wFieldsNew_, cFieldsNew_, residualNew_);
         normNew = residualNorm(residualNew_);

         // Decrease increment if necessary
         j = 0;
         while (normNew > norm && j < 5) {
            std::cout << "      decreasing increment,  norm = " 
                      << normNew << std::endl;
            needsJacobian = true;
            for (k = 0; k < nr; ++k) {
               dOmega_[k] *= 0.66666666;
            }
            incrementWFields(system().wFields(), dOmega_, wFieldsNew_);
            mixture().compute(wFieldsNew_, cFieldsNew_);
            computeResidual(wFieldsNew_, cFieldsNew_, residualNew_);
            normNew = residualNorm(residualNew_);
            ++j;
         }
 
         // Accept or reject update
         if (normNew < norm) {
            for (j = 0; j < nm; ++j) {
               for (k = 0; k < nx; ++k) {
                  system().wField(j)[k] = wFieldsNew_[j][k];
                  system().cField(j)[k] = cFieldsNew_[j][k];
               }
            }
            system().computeFreeEnergy();
            hasJacobian = false;
            if (normNew/norm < 0.25) {
               needsJacobian = false;
            } else {
               needsJacobian = true;
            }
            for (j = 0; j < nr; ++j) {
               residual_[j] = residualNew_[j];
            }
            norm = normNew;
         } else {
            std::cout << "Iteration failed, norm = " 
                      << normNew << std::endl;
            if (hasJacobian) {
               return 1;
               std::cout << "Unrecoverable failure " << std::endl;
            } else {
               std::cout << "Try rebuilding Jacobian" << std::endl;
               needsJacobian = true;
            }
         }

      }

      // Failure 
      return 1;
   }

} // namespace Fd1d
} // namespace Pscf

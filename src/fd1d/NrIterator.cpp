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
      if (system().domain().nx() > 0) {
         allocate();
      }
   }

   void NrIterator::allocate()
   {
      int nm = mixture().nMonomer();   // number of monomer types
      int nx = system().domain().nx(); // number of grid points
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
      //std::cout << "Begin computeResidual .. ";
      int nm = mixture().nMonomer();   // number of monomer types
      int nx = system().domain().nx(); // number of grid points
      int i;                           // grid point index
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

      //std::cout << "Finish computingRresidual" << std::endl;
   }

   void NrIterator::computeJacobian()
   {
      // std::cout << "Begin computeJacobian ... ";
      int nm = mixture().nMonomer();   // number of monomer types
      int nx = system().domain().nx(); // number of grid points
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

      // std::cout << "Finishing computeJacobian" << std::endl;
   }

   void NrIterator::update()
   {
      // std::cout << "Begin update .. ";
      computeJacobian();

      // Compute increment dOmega_
      solver_.solve(residual_, dOmega_);

      // Increment wFields;
      int nm = mixture().nMonomer();   // number of monomers types
      int nx = system().domain().nx(); // number of grid points
      int i;                           // monomer index
      int j;                           // grid point index
      int k = 0;                       // residual element index
      for (i = 0; i < nm; ++i) {
         for (j = 0; j < nx; ++j) {
            system().wField(i)[j] -= dOmega_[k];
            ++k;
         }
      }

      #if 0
      std::cout << "WFields" << std::endl;
      for (j = 0; j < nx; ++j) {
         std::cout << j << "  ";
         for (i = 0; i < nm; ++i) {
            std::cout << system().wField(i)[j] << "  ";
         }
         std::cout << std::endl;
      }
      #endif

      mixture().compute(system().wFields(), system().cFields());
      computeResidual(system().wFields(), system().cFields(), residual_);

      #if 0
      std::cout << "CFields" << std::endl;
      for (j = 0; j < nx; ++j) {
         std::cout << j << "  ";
         for (i = 0; i < nm; ++i) {
            std::cout << system().cField(i)[j] << "  ";
         }
         std::cout << std::endl;
      }
      #endif

      #if 0
      std::cout << "Residuals" << std::endl;
      int nr = nm*nx;
      for (int ir = 0; ir <  nr; ++ir) {
         std::cout << ir << "  " << residual_[ir] << std::endl;
      }
      #endif

      //std::cout << "Finish update" << std::endl;
   }

   bool NrIterator::isConverged()
   {
      int nm = mixture().nMonomer();    // number of monomer types
      int nx = system().domain().nx();  // number of grid points
      int nr = nm*nx;                   // number of residual components
      for (int ir = 0; ir <  nr; ++ir) {
         if (fabs(residual_[ir]) > epsilon_) {
            return false;
         }
      }
      std::cout << "Converged" << std::endl;
      return true;
   }

   int NrIterator::solve()
   {
      // Allocate memory if needed or, if allocated, check array sizes.
      allocate();

      // Compute current residual
      mixture().compute(system().wFields(), system().cFields());
      computeResidual(system().wFields(), system().cFields(), residual_);

      // Iterative loop
      
      for (int i = 0; i < 100; ++i) {
         std::cout << "Begin iteration " << i << std::endl;
         if (isConverged()) {
            return 0;
         } else {
            update();
         }
      }
      return 1;
   }

} // namespace Fd1d
} // namespace Pscf

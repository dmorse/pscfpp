/*
* PFTS - Polymer Field Theory Simulator
*
* Copyright 2013, David Morse (morse012@.umn.edu)
* Distributed under the terms of the GNU General Public License.
*/

#include "NrIterator.h"
#include "System.h"
#include "Mixture.h"

namespace Pscf {
namespace Fd1d
{

   using namespace Util;

   NrIterator::NrIterator()
    : Iterator(),
      epsilon_(0.0)
   {  setClassName("NrIterator"); }

   NrIterator::~NrIterator()
   {}

   void NrIterator::readParameters(std::istream& in)
   {
      read(in, "epsilon", epsilon_);

      // Allocate work space
      int nm = mixture().nMonomer();   // number of  monomers
      int nx = system().grid().nx();   // number of grid points
      int nr = nm*nx;                  // number of residual components
      residual_.allocate(nr);
      jacobian_.allocate(nr, nr);
      residualNew_.allocate(nr);
      dOmega_.allocate(nr);
      wFieldsNew_.allocate(nm);
      for (int i = 0; i < nm; ++i) {
         wFieldsNew_[i].allocate(nx);
      }
   }

   void NrIterator::computeResidual(Array<WField> const & wFields, 
                                    Array<double>& residual)
   {}

   void NrIterator::computeJacobian()
   {
      // Compute residual for current fields
      computeResidual(system().wFields(), residual_);

      int nm = mixture().nMonomer();  // number of  monomers
      int nx = system().grid().nx();  // number of grid points
      int i;                          // monomer index
      int j;                          // grid point index

      // Copy system().wFields to wFieldsNew.
      for (i = 0; i < nm; ++i) {
         for (j = 0; i < nx; ++j) {
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
            computeResidual(wFieldsNew_, residualNew_);
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

   void NrIterator::update()
   {
      computeJacobian();

      // Compute increment dOmega_
      solver_.solve(residual_, dOmega_);

      // Increment wFields;
      int nm = mixture().nMonomer();  // number of  monomers
      int nx = system().grid().nx();  // number of grid points
      int i;                          // monomer index
      int j;                          // grid point index
      int k = 0;                      // residual element index
      for (i = 0; i < nm; ++i) {
         for (j = 0; j < nx; ++j) {
            system().wField(i)[j] = dOmega_[k];
            ++k;
         }
      }

   }

   bool NrIterator::isConverged()
   {
      bool criterion = true;
      int nm = mixture().nMonomer();  // number of  monomers
      int nx = system().grid().nx();  // number of grid points
      int nr = nm*nx;  // number of residual components
      for (int ir = 0; ir <  nr; ++ir) {
         if (abs(residual_[ir]) > epsilon_) {
            criterion = false;
            break;
         }
      }
      return criterion;
   }

   int NrIterator::solve()
   {
      computeResidual(system().wFields(), residual_);
      return 0;
   }

} // namespace Fd1d
} // namespace Pscf

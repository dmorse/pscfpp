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
      int i, j;
      int nm = mixture().nMonomer();         // number of  monomers
      int nx = system().grid().nx();   // number of grid points
      int nr = nm*nx;                  // number of residual components
      computeResidual(system().wFields(), residual_);

      // Copy system().wFields to wFieldsNew.
      for (i = 0; i < nm; ++i) {
         for (j = 0; i < nx; ++j) {
            wFieldsNew_[i][j] = system().wField(i)[j];
         }
      }

      double delta = 0.001;
      int jr, jc;
      jc = 0;
      for (i = 0; i < nm; ++i) {
         for (j = 0; i < nx; ++j) {
             wFieldsNew_[i][j] += delta;
             computeResidual(wFieldsNew_, residualNew_);
             for (jr=0; jr < nr; ++jr) {
                jacobian_(jr, jc) = (residualNew_[jr] - residual_[jr])/delta;
             }
             wFieldsNew_[i][j] = system().wField(i)[j];
         }
      }
      #if 0
      // LU decompose Jacobian
      #endif
   }

   void NrIterator::update()
   {
      #if 0
      Solve for update;
      Increment wFields;
      #endif
   }

   int NrIterator::solve()
   {
      return 0;
   }

} // namespace Fd1d
} // namespace Pscf

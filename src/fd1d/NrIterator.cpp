/*
* PFTS - Polymer Field Theory Simulator
*
* Copyright 2013, David Morse (morse012@.umn.edu)
* Distributed under the terms of the GNU General Public License.
*/

#include "NrIterator.h"

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
   }

   void NrIterator::computeResidual(Array<WField> const & wFields, 
                                    Array<double>& residual)
   {}

   void NrIterator::computeJacobian()
   {
      #if 0
      Copy mixture().wFields to wFieldsNew.
      for (each monomer type) {
         for (each grid point) {
             Set relevant element of wFields new.
             computeResidual(wFieldsNew, residualNew);
             Copy difference in residuals to Jacobian row;
             Reset relevant element to original value;
         }
      }
      LU decompose Jacobian
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

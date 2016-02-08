#ifndef FD1D_NR_ITERATOR_H
#define FD1D_NR_ITERATOR_H

/*
* PFTS - Polymer Field Theory Simulator
*
* Copyright 2013, David Morse (morse012@.umn.edu)
* Distributed under the terms of the GNU General Public License.
*/

#include "Iterator.h"
#include <util/containers/DArray.h>
#include <util/containers/DMatrix.h>

namespace Pscf {
namespace Fd1d
{

   using namespace Util;

   /**
   * Newton-Raphson Iterator.
   *
   * \ingroup Pscf_Fd1d_Module
   */
   class NrIterator : public Iterator
   {

   public:

      /**
      * Constructor.
      */
      NrIterator();

      /**
      * Destructor.
      */
      virtual ~NrIterator();

      /**
      * Read all parameters and initialize.
      *
      * \param in input parameter stream
      */
      void readParameters(std::istream& in);

      /**
      * Iterate self-consistent field equations to solution.
      * 
      * \return error code: 0 for success, 1 for failure.
      */
      int solve();

      /**
      * Get error tolerance.
      */
      double epsilon();

   private:

      /// Residual vector. size = (# monomers)x(# grid points).
      DArray<double> residual_;

      /// Jacobian matrix
      DMatrix<double> Jacobian_;

      /// Error tolerance.
      double epsilon_;

   };

   // Inline function

   inline double NrIterator::epsilon()
   {  return epsilon_; }

} // namespace Fd1d
} // namespace Pscf
#endif

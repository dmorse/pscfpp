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
      */
      void readParameters(std::istream& in);

      /**
      * Iterate to solution.
      * 
      * \return error code: 0 for success, 1 for failure.
      */
      int solve();

      /**
      * Get error tolerance.
      */
      double epsilon();

   private:

      DArray<double> residual_;

      DMatrix<double> Jacobian_;

      // Error tolerance.
      double epsilon_;

   };

   inline double NrIterator::epsilon()
   {  return epsilon_; }

} // namespace Fd1d
} // namespace Pscf
#endif

#ifndef PRDC_CPU_REDUCE_H
#define PRDC_CPU_REDUCE_H

/*
* PSCF Package
*
* Copyright 2015 - 2025, The Regents of the University of Minnesota
* Distributed under the terms of the GNU General Public License.
*/

#include <util/containers/Array.h>

using namespace Util;

namespace Pscf {
namespace Prdc {
namespace Cpu {

   /**
   * Functions that perform array reductions on the Cpu.
   *
   * A reduction is any operation that involves reducing all of 
   * the elements of an array or set of arrays to a single scalar.  
   * Examples include taking the sum or finding the maximum of all 
   * array elements, or taking an inner product of two arrays.
   *
   * \ingroup Prdc_Cpu_Module
   * @{
   */
   namespace Reduce {

      /**
      * Get maximum of array elements .
      *
      * \param in  input array
      */
      double max(Array<double> const & in);

      /**
      * Get maximum absolute magnitude of array elements .
      *
      * \param in  input array
      */
      double maxAbs(Array<double> const & in);

      /**
      * Get minimum of array elements .
      *
      * \param in  input array
      */
      double min(Array<double> const & in);

      /**
      * Get minimum absolute magnitude of array elements .
      *
      * \param in  input array
      */
      double minAbs(Array<double> const & in);

      /**
      * Compute sum of array elements .
      *
      * \param in  input array
      */
      double sum(Array<double> const & in);

      /**
      * Compute inner product of two real arrays .
      *
      * \param a  first input array
      * \param b  second input array
      */
      double innerProduct(Array<double> const & a,
                          Array<double> const & b);

   /** @} */

} // Reduce
} // Cpu
} // Prdc
} // Pscf
#endif

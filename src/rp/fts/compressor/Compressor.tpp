#ifndef RP_COMPRESSOR_TPP
#define RP_COMPRESSOR_TPP

/*
* PSCF - Polymer Self-Consistent Field
*
* Copyright 2015 - 2025, The Regents of the University of Minnesota
* Distributed under the terms of the GNU General Public License.
*/

#include "Compressor.h"

namespace Pscf {
namespace Rp {

   /*
   * Default constructor.
   */
   template <int D, class T>
   Compressor<D,T>::Compressor()
    : mdeCounter_(0),
      sysPtr_(nullptr)
   {  setClassName("Compressor"); }

   /*
   * Constructor (creates association with parent system)
   */
   template <int D, class T>
   Compressor<D,T>::Compressor(System<D,T>& system)
    : mdeCounter_(0),
      sysPtr_(&system)
   {  setClassName("Compressor"); }

   /*
   * Create association with the parent system.
   */
   template <int D, class T>
   void Compressor<D,T>::setSystem(System<D,T>& system)
   {  sysPtr_ = &system; }

   /*
   * Get number of times MDE has been solved.
   */
   template <int D, class T>
   int Compressor<D,T>::mdeCounter() const
   {  return mdeCounter_; }

} // namespace Rp
} // namespace Pscf
#endif

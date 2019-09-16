/*
* PSCF - Polymer Self-Consistent Field Theory
*
* Copyright 2016, The Regents of the University of Minnesota
* Distributed under the terms of the GNU General Public License.
*/

#include "Iterator.h"

namespace Pscf {
namespace Pssp_gpu
{

   using namespace Util;

   template<int D>
   Iterator<D>::Iterator()
   {  setClassName("Iterator"); }

   template<int D>
   Iterator<D>::Iterator(System<D>* system)
    : systemPtr_(system)
   {  setClassName("Iterator"); }

   template<int D>
   Iterator<D>::~Iterator()
   {}

} // namespace Pssp_gpu
} // namespace Pscf

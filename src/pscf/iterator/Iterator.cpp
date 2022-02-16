#ifndef PSCF_ITERATOR_CPP
#define PSCF_ITERATOR_CPP

/*
* PSCF - Polymer Self-Consistent Field Theory
*
* Copyright 2016 - 2019, The Regents of the University of Minnesota
* Distributed under the terms of the GNU General Public License.
*/

#include "Iterator.h"

namespace Pscf {
namespace Pspc
{

   using namespace Util;

   Iterator::Iterator(IteratorMediator<T>& iterMed)
    : isFlexible_(false),
      iterMed_(&iterMed)
   {  setClassName("Iterator"); }

   Iterator::~Iterator()
   {}

} // namespace Pspc
} // namespace Pscf
#endif

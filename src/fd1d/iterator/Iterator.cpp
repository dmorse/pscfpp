/*
* PSCF - Polymer Self-Consistent Field Theory
*
* Copyright 2016 - 2019, The Regents of the University of Minnesota
* Distributed under the terms of the GNU General Public License.
*/

#include "Iterator.h"
#include <fd1d/System.h>
#include <fd1d/domain/Domain.h>
#include <fd1d/solvers/Mixture.h>

namespace Pscf {
namespace Fd1d
{

   using namespace Util;

   Iterator::Iterator()
   {  setClassName("Iterator"); }

   Iterator::Iterator(System& system)
    : SystemAccess(system)
   {  setClassName("Iterator"); }

   Iterator::~Iterator()
   {}

} // namespace Fd1d
} // namespace Pscf

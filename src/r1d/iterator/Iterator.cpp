/*
* PSCF - Polymer Self-Consistent Field Theory
*
* Copyright 2015 - 2025, The Regents of the University of Minnesota
* Distributed under the terms of the GNU General Public License.
*/

#include "Iterator.h"
#include <r1d/System.h>
#include <r1d/domain/Domain.h>
#include <r1d/solvers/Mixture.h>

namespace Pscf {
namespace R1d
{

   using namespace Util;

   Iterator::Iterator()
   {  setClassName("Iterator"); }

   Iterator::Iterator(System& system)
    : SystemAccess(system)
   {  setClassName("Iterator"); }

   Iterator::~Iterator()
   {}

} // namespace R1d
} // namespace Pscf

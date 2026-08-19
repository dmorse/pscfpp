/*
* PSCF - Polymer Self-Consistent Field
*
* Copyright 2015 - 2026, The Regents of the University of Minnesota
* Distributed under the terms of the GNU General Public License.
*/

#include "Iterator.h"
#include <r1d/system/System.h>
#include <r1d/field/Domain.h>
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

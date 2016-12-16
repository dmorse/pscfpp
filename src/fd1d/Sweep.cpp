/*
* PSCF - Polymer Self-Consistent Field Theory
*
* Copyright 2016, The Regents of the University of Minnesota
* Distributed under the terms of the GNU General Public License.
*/

#include "Sweep.h"
#include "System.h"
#include "Mixture.h"
#include "Domain.h"
#include "Iterator.h"

namespace Pscf {
namespace Fd1d
{

   using namespace Util;

   Sweep::Sweep()
    : systemPtr_(0),
      mixturePtr_(0),
      domainPtr_(0),
      iteratorPtr_(0)
   {  setClassName("Sweep"); }

   Sweep::~Sweep()
   {}

   void Sweep::setSystem(System& system)
   {
      systemPtr_  = &system;
      mixturePtr_ = &(system.mixture());
      domainPtr_ = &(system.domain());
      iteratorPtr_ = &(system.iterator());
   }

} // namespace Fd1d
} // namespace Pscf

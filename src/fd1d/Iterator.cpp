/*
* PSCF - Polymer Self-Consistent Field Theory
*
* Copyright 2016, The Regents of the University of Minnesota
* Distributed under the terms of the GNU General Public License.
*/

#include "Iterator.h"
#include "System.h"
#include "Mixture.h"

namespace Pscf {
namespace Fd1d
{

   using namespace Util;

   Iterator::Iterator()
    : systemPtr_(0),
      mixturePtr_(0)
   {  setClassName("Iterator"); }

   Iterator::~Iterator()
   {}

   void Iterator::setSystem(System& system)
   {  
      systemPtr_  = &system; 
      mixturePtr_ = &(system.mixture()); 
   }

} // namespace Fd1d
} // namespace Pscf

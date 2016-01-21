/*
* PFTS - Polymer Field Theory Simulator
*
* Copyright 2013, David Morse (morse012@.umn.edu)
* Distributed under the terms of the GNU General Public License.
*/

#include "Iterator.h"

namespace Pscf {
namespace Fd1d
{

   using namespace Util;

   Iterator::Iterator()
    : mixturePtr_(0)
   {  setClassName("Iterator"); }

   void Iterator::setMixture(Mixture & mixture )
   {  mixturePtr_ = &mixture; }

} // namespace Fd1d
} // namespace Pscf

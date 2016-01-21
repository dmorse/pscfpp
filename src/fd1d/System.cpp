/*
* PFTS - Polymer Field Theory Simulator
*
* Copyright 2013, David Morse (morse012@.umn.edu)
* Distributed under the terms of the GNU General Public License.
*/

#include "System.h"
#include "NrIterator.h"

namespace Pscf {
namespace Fd1d
{

   using namespace Util;

   System::System()
    : mixture_(),
      iteratorPtr_(0)
   {  
      setClassName("System"); 
      iteratorPtr_ = new NrIterator(); 
   }

   void System::readParameters(std::istream& in)
   {
   }  

} // namespace Fd1d
} // namespace Pscf

/*
* Simpatico - Simulation Package for Polymeric and Molecular Liquids
*
* Copyright 2010 - 2014, The Regents of the University of Minnesota
* Distributed under the terms of the GNU General Public License.
*/

#include "initStatic.h"
#include <util/global.h>

#include <util/misc/Log.h>
#include <util/misc/Memory.h>
#include <util/format/Format.h>
#include <util/param/ParamComponent.h>
#include <util/math/Constants.h>
#include <util/space/Vector.h>
#include <util/space/IntVector.h>
#include <util/space/Tensor.h>

namespace Util
{

   /*
   * Initialize all static class members in Util namespace.
   */
   void initStatic()
   {
      // Precondition: This function can only be called once.
      static int nCall = 0;
      if (nCall == 0) {
         // Call initStatic() methods of all relevant classes.
         Log::initStatic();
         Format::initStatic();
         Memory::initStatic();
         ParamComponent::initStatic();
         Constants::initStatic();
         Vector::initStatic();
         IntVector::initStatic();
         Tensor::initStatic();
      }
      ++nCall;
   }

}

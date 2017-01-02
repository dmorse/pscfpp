/*
* Simpatico - Simulation Package for Polymeric and Molecular Liquids
*
* Copyright 2010 - 2014, The Regents of the University of Minnesota
* Distributed under the terms of the GNU General Public License.
*/

#include "End.h"
#include <util/global.h>

namespace Util
{

   /* 
   * Constructor.
   */
   End::End()
    : label_("}")
   {}

   /* 
   * Destructor.
   */
   End::~End()
   {}

   /* 
   * Read and check end bracket.
   */
   void End::readParam(std::istream &in)
   {
      if (isIoProcessor()) {
         in >> label_;
         if (ParamComponent::echo()) {
            writeParam(Log::file());
         }
      }
   }

   /* 
   * End::writeParam()
   */
   void End::writeParam(std::ostream &out)
   {  out << indent() << "}" << std::endl; }

   /* 
   * Empty implementation of virtual resetParam() method.
   */
   void End::resetParam()
   {}

} 

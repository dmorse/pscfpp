/*
* Simpatico - Simulation Package for Polymeric and Molecular Liquids
*
* Copyright 2010 - 2014, The Regents of the University of Minnesota
* Distributed under the terms of the GNU General Public License.
*/

#include "Begin.h"
#ifdef UTIL_MPI
#include <util/mpi/MpiSendRecv.h>
#endif

#include <util/global.h>

namespace Util
{

   /* 
   * Begin constructor
   */
   Begin::Begin(const char *label, bool isRequired)
    : label_(isRequired),
      isActive_(false)
   {
      std::string expected = label;
      expected += "{";
      label_.setString(expected);
   }

   /* 
   * Read label, check against expected value
   */
   void Begin::readParam(std::istream &in)
   {
      if (isIoProcessor()) {
         in >> label_;

         // If this parameter is required and the label string
         // doesn't match, an Exception will be thrown by the
         // operator >> for a Label, terminating this function.

         if (Label::isClear()) {
            // If the label string matches
            isActive_ = true;
            if (ParamComponent::echo()) {
               writeParam(Log::file());
            }
         } else {
            // If label does not match and this isOptional
            isActive_ = false;
            if (ParamComponent::echo()) {
               Log::file() << indent() 
                           << label_.string() << " [absent] }"
                           << std::endl; 
            }
         }
      } else {
         #ifdef UTIL_MPI
         if (!hasIoCommunicator()) {
            UTIL_THROW("Error: not isIoProcessor and not hasIoCommunicator");
         }
         #else
         UTIL_THROW("Error: not isIoProcessor and no MPI");
         #endif
      }
      #ifdef UTIL_MPI
      if (hasIoCommunicator()) {
         if (isRequired()) {
            isActive_ = true;
         } else {
            bcast<bool>(ioCommunicator(), isActive_, 0); 
         }
      }
      #endif
   }

   /* 
   * Begin::writeParam() template
   */
   void Begin::writeParam(std::ostream &out)
   {  out << indent() << label_.string() << std::endl; }

   /*
   * Do-nothing implementation of virtual resetIo function.
   */
   void Begin::resetParam()
   {}

} 

/*
* Simpatico - Simulation Package for Polymeric and Molecular Liquids
*
* Copyright 2010 - 2014, The Regents of the University of Minnesota
* Distributed under the terms of the GNU General Public License.
*/

#include "Parameter.h"

#ifdef UTIL_MPI
#include <util/mpi/MpiSendRecv.h>
#endif

#include <iomanip>

namespace Util
{

   /*
   * Constructor. 
   */
   Parameter::Parameter(const char *label, bool isRequired)
    : label_(label, isRequired),
      isActive_(isRequired)
   {}

   /*
   * Destructor.
   */
   Parameter::~Parameter()
   {}

   /*
   * Read a parameter.
   */
   void Parameter::readParam(std::istream &in)
   {
      if (isIoProcessor()) {
         in >> label_;

         // If this parameter is required and the 
         // label doesn't match an exception will 
         // be thrown by the >> operator.

         if (Label::isClear()) {
            // If the label string matches
            readValue(in);
            isActive_ = true;
            if (ParamComponent::echo()) {
               writeParam(Log::file());
            }
         } else {
            // If label does not match and this isOptional
            isActive_ = false;
            if (ParamComponent::echo() && !isRequired()) {
               Log::file() << indent() 
                           << label_ << std::right
                           << std::setw(Parameter::Width)
                           << "[ absent ]" << std::endl;
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
            bcastValue();
            isActive_ = true;
         } else {
            bcast<bool>(ioCommunicator(), isActive_, 0); 
            if (isActive_) {
               bcastValue();
            }
         }
      }
      #endif
   }

   /*
   * Load from an archive.
   */
   void Parameter::load(Serializable::IArchive& ar)
   {
      if (isIoProcessor()) {
         if (isRequired()) { 
            isActive_ = true;
         } else {
            ar >> isActive_;
         }
         if (isActive_) {
            loadValue(ar);
            if (ParamComponent::echo()) {
               writeParam(Log::file());
            }
         } else {
            if (ParamComponent::echo() && !isRequired()) {
               Log::file() << indent() 
                           << label_ << std::right
                           << std::setw(Parameter::Width)
                           << "[ absent ]" << std::endl;
            }
         }
      } else {
         #ifdef UTIL_MPI
         if (!hasIoCommunicator()) {
            UTIL_THROW("Error: not isIoProcessor and !hasIoCommunicator");
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
         if (isActive_) {
            bcastValue();
         }
      }
      #endif
   }

   /*
   * Save to an archive.
   */
   void Parameter::save(Serializable::OArchive& ar)
   {
      if (!isRequired()) {
         ar << isActive_;
      }
      if (isActive_) {
         saveValue(ar);
      }
   }

   /*
   * Return label string.
   */
   std::string Parameter::label() const
   {  return label_.string(); }

   /*
   * Is this a required parameter?
   */
   bool Parameter::isRequired() const
   {  return label_.isRequired(); }

   /*
   * Is this an active parameter?
   */
   bool Parameter::isActive() const
   {  return isActive_; }

} 

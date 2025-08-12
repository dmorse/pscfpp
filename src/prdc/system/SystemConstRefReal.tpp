#ifndef PRDC_SYSTEM_CONST_REF_REAL_TPP
#define PRDC_SYSTEM_CONST_REF_REAL_TPP

/*
* PSCF - Polymer Self-Consistent Field
*
* Copyright 2015 - 2025, The Regents of the University of Minnesota
* Distributed under the terms of the GNU General Public License.
*/

#include "SystemConstRefReal.h"

namespace Pscf {
namespace Prdc {

   /*
   * Default constructor.
   */
   template <class ST>
   SystemConstRefReal<ST>::SystemConstRefReal()
    : systemPtr_(nullptr),
      mixturePtr_(nullptr),
      interactionPtr_(nullptr),
      domainPtr_(nullptr),
      cPtr_(nullptr),
      wPtr_(nullptr),
      hPtr_(nullptr),
      maskPtr_(nullptr),
      fileMasterPtr_(nullptr)
   {}

   /*
   * Constructor (creates associations).
   */
   template <class ST>
   SystemConstRefReal<ST>::SystemConstRefReal(SystemT const & system)
    : systemPtr_(nullptr),
      mixturePtr_(nullptr),
      interactionPtr_(nullptr),
      domainPtr_(nullptr),
      cPtr_(nullptr),
      wPtr_(nullptr),
      hPtr_(nullptr),
      maskPtr_(nullptr),
      fileMasterPtr_(nullptr)
   {  associate(system); }

   /*
   * Destructor.
   */
   template <class ST>
   SystemConstRefReal<ST>::~SystemConstRefReal()
   {}

   template <class ST>
   void SystemConstRefReal<ST>::associate(SystemT const & system)
   {
      systemPtr_ = &system;
      mixturePtr_ = &(system.mixture());
      interactionPtr_ = &(system.interaction());
      domainPtr_ = &(system.domain());
      cPtr_ = &(system.c());
      wPtr_ = &(system.w());
      hPtr_ = &(system.h());
      maskPtr_ = &(system.mask());
      fileMasterPtr_ = &(system.fileMaster());
   }

} // namespace Prdc
} // namespace Pscf
#endif

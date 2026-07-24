#ifndef RP_SYSTEM_CONST_REF_TPP
#define RP_SYSTEM_CONST_REF_TPP

/*
* PSCF - Polymer Self-Consistent Field
*
* Copyright 2015 - 2025, The Regents of the University of Minnesota
* Distributed under the terms of the GNU General Public License.
*/

#include <rp/system/System.h>
#include <rp/system/SystemConstRef.h>

namespace Pscf {
namespace Rp {

   /*
   * Default constructor.
   */
   template <int D, class T>
   SystemConstRef<D,T>::SystemConstRef()
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
   template <int D, class T>
   SystemConstRef<D,T>::SystemConstRef(System<D,T> const & system)
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
   template <int D, class T>
   SystemConstRef<D,T>::~SystemConstRef()
   {}

   template <int D, class T>
   void SystemConstRef<D,T>::associate(System<D,T> const & system)
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

} // namespace Rp
} // namespace Pscf
#endif

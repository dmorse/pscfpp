#ifndef PSCF_ABSTRACT_SYSTEM_H
#define PSCF_ABSTRACT_SYSTEM_H

/*
* PSCF - Polymer Self-Consistent Field Theory
*
* Copyright 2016 - 2019, The Regents of the University of Minnesota
* Distributed under the terms of the GNU General Public License.
*/

#include <util/param/ParamComposite.h>     // base class

namespace Pscf {

   using namespace Util;

   /**
   * Abstract system class for SCFT simulation of one system,
   * independent of implementation.
   *
   * \ingroup Pscf_Module
   */

   class AbstractSystem : public ParamComposite
   {

   public:

   private:
   
   };

} // namespace Pscf
#endif

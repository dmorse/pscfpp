#ifndef PSCF_SOLVENT_TMPL_H
#define PSCF_SOLVENT_TMPL_H

/*
* PSCF - Polymer Self-Consistent Field Theory
*
* Copyright 2013, David Morse (morse012@.umn.edu)
* Distributed under the terms of the GNU General Public License.
*/

#include <pscf/Species.h>
#include <util/param/ParamComposite.h>

namespace Pscf
{ 

   using namespace Util;

   /**
   * Template for a class representing a solvent species.
   *
   * \ingroup Pscf_Base_Module
   */
   template <class CField>
   class SolventTmpl : public Species, public ParamComposite
   {
   public:
   
      /**
      * Get monomer concentration field for this solvent.
      */
      const CField& concentration() const;
   
   protected:
   
      CField concentration_;
   
   };

} 
#endif 

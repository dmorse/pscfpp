#ifndef PFTS_SOLVENT_TMPL_H
#define PFTS_SOLVENT_TMPL_H

/*
* PFTS - Polymer Field Theory Simulator
*
* Copyright 2013, David Morse (morse012@.umn.edu)
* Distributed under the terms of the GNU General Public License.
*/

#include <pfts/Species.h>
#include <util/param/ParamComposite.h>

namespace Pfts
{ 

   using namespace Util;

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

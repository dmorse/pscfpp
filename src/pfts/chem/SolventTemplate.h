#ifndef PFTS_SOLVENT_TEMPLATE_H
#define PFTS_SOLVENT_TEMPLATE_H

/*
* PFTS - Polymer Field Theory Simulator
*
* Copyright 2013, David Morse (morse012@.umn.edu)
* Distributed under the terms of the GNU General Public License.
*/

namespace Pfts{ 

   using namespace Util;


   template <class CField>
   class SolventTemplate : public Species, public SolventDescriptor
   {
   public:
   
      /**
      * Get monomer concentration field for this solvent.
      */
      const C& concentration() const;
   
   protected:
   
      CField concentration_;
   
   };

} 
#endif 

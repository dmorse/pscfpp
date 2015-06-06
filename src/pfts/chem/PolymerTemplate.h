#ifndef PFTS_TY_H
#define PFTS_TY_H

/*
* PFTS - Polymer Field Theory Simulator
*
* Copyright 2013, David Morse (morse012@.umn.edu)
* Distributed under the terms of the GNU General Public License.
*/

namespace Pfts{ 

   using namespace Util;

   template <class Propagator, class CField>
   class PolymerTemplate : public Species, public PolymerDescriptor
   {
   public:
   
      /**
      * Read parameters and initialize.
      */
      virtual void readParameters();
   
      /**
      * Get monomer concentration field for specific block.
      */
      const CField& blockCField(int blockId) const;
   
      /**
      * Get propagator for a specific block and direction.
      */
      const Propagator& propagator(int blockId, int direction) const;
   
   protected:
   
      void init();
   
   private:
   
      DArray<Pair<Propagator>> propagators_;
   
      DArray<CField> blockCFields_;
   
   };

} 
#endif 

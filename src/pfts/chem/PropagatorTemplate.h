#ifndef PFTS_PROPAGATOR_TEMPLATE_H
#define PFTS_PROPAGATOR_TEMPLATE_H

/*
* PFTS - Polymer Field Theory Simulator
*
* Copyright 2013, David Morse (morse012@.umn.edu)
* Distributed under the terms of the GNU General Public License.
*/

namespace Pfts{ 

   using namespace Util;

   template <class W, class Q>
   class PropagatorTemplate
   {
   public:
   
      void setBlock(Block& block, int sign);
      void addSource(Propagator<W, Q>& other);
      void init();
      void compute(W& w);
   
      const Block& block() const;
      const Q& tail() const;
   
   private:
   
      Block* blockPtr_;
      int sign_;
      GArray<Propagator<W, Q>*> sourcesPtrs_;
   
   };

} 
#endif 

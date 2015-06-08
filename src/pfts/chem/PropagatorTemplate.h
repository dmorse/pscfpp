#ifndef PFTS_PROPAGATOR_TEMPLATE_H
#define PFTS_PROPAGATOR_TEMPLATE_H

/*
* PFTS - Polymer Field Theory Simulator
*
* Copyright 2013, David Morse (morse012@.umn.edu)
* Distributed under the terms of the GNU General Public License.
*/

#include "Block.h"
#include <util/containers/GArray.h>

namespace Pfts{ 

   using namespace Util;

   template <class TPropagator>
   class PropagatorTemplate
   {

   public:
   
      void setBlock(Block& block, int directionId);

      void addSource(TPropagator& other);

      virtual void clear();
   
      // virtual void compute(TWField& w) = 0;
   
      // virtual const TQField& tail() const = 0;
   
      const Block& block() const;

      int directionId() const;

      int nSource() const;

      const TPropagator& source(int id) const;
 
   private:
  
      /// Pointer to associated block.
      Block* blockPtr_;

      /// Direction of propagation (0 or 1)
      int directionId_;

      /// Pointers to propagators that feed source vertex.
      GArray<TPropagator*> sourcePtrs_;

      /// Set true by compute function and false by clear.
      bool isComputed_;
   
   };

   template <class TPropagator>
   const Block& PropagatorTemplate<TPropagator>::block() const
   {  return *blockPtr_; }

   template <class TPropagator>
   int PropagatorTemplate<TPropagator>::directionId() const
   {  return directionId_; }

   template <class TPropagator>
   void PropagatorTemplate<TPropagator>::clear()
   {  isComputed = false_; }

   template <class TPropagator>
   int PropagatorTemplate<TPropagator>::nSource() const
   {  return sourcePtrs_.size(); }

   template <class TPropagator>
   const TPropagator& 
   PropagatorTemplate<TPropagator>::source(int id) const
   {  return *(sourcePtrs_[id]); }

   template <class TPropagator>
   void PropagatorTemplate<TPropagator>::setBlock(Block& block, 
                                                  int directionId);
   {  
      blockPtr_ = &block;
      directionId_ = directionId; 
   }

   template <class TPropagator>
   void PropagatorTemplate<TPropagator>::addSource(TPropagator& other)
   {  sourcePtrs_.append(&other); }

} 
#endif 

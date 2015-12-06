#ifndef PFTS_PROPAGATOR_TMPL_H
#define PFTS_PROPAGATOR_TMPL_H

/*
* PFTS - Polymer Field Theory Simulator
*
* Copyright 2013, David Morse (morse012@.umn.edu)
* Distributed under the terms of the GNU General Public License.
*/

#include <util/containers/GArray.h>
#include <pfts/chem/Block.h>

namespace Pfts{ 

   using namespace Util;

   template <class TPropagator>
   class PropagatorTmpl
   {

   public:
  
      /**
      * Associate this propagator with a block and direction.
      *
      * \param block associated Block object.
      * \param directionId direction = 0 or 1.
      */ 
      void setBlock(Block& block, int directionId);

      /**
      * Add another propagator to the list of sources for this one.
      *
      * A source is a propagator that is needed to compute the initial
      * condition for this one, and that thus must be computed before 
      * this one.
      */
      void addSource(const TPropagator& other);

      /**
      * Set the isComplete flag to false.
      */
      virtual void clear();
 
      /**
      * Get associated Block object by reference.
      */
      const Block& block() const;

      /**
      * Get direction index for this propagator.
      */
      int directionId() const;

      /**
      * Number of source / prerequisite propagators.
      */
      int nSource() const;

      /**
      * Get a source propagator.
      */
      const TPropagator& source(int id) const;

      /**
      * Has the modified diffusion equation been solved?
      */
      bool isComplete() const;
 
   private:
  
      /// Pointer to associated block.
      Block* blockPtr_;

      /// Direction of propagation (0 or 1)
      int directionId_;

      /// Pointers to propagators that feed source vertex.
      GArray<TPropagator const *> sourcePtrs_;

      /// Set true by compute function and false by clear.
      bool isComplete_;
   
   };

   template <class TPropagator>
   inline void PropagatorTmpl<TPropagator>::clear()
   {  isComplete_ = false; }

   template <class TPropagator>
   inline const Block& PropagatorTmpl<TPropagator>::block() const
   {  return *blockPtr_; }

   template <class TPropagator>
   inline int PropagatorTmpl<TPropagator>::directionId() const
   {  return directionId_; }

   template <class TPropagator>
   inline int PropagatorTmpl<TPropagator>::nSource() const
   {  return sourcePtrs_.size(); }

   template <class TPropagator>
   inline const TPropagator& 
   PropagatorTmpl<TPropagator>::source(int id) const
   {  return *(sourcePtrs_[id]); }

   template <class TPropagator>
   bool PropagatorTmpl<TPropagator>::isComplete() const
   {  return isComplete_; }

   template <class TPropagator>
   void PropagatorTmpl<TPropagator>::setBlock(Block& block, 
                                              int directionId)
   {
      blockPtr_ = &block;
      directionId_ = directionId; 
   }

   template <class TPropagator>
   void PropagatorTmpl<TPropagator>::addSource(const TPropagator& other)
   {  sourcePtrs_.append(&other); }

} 
#endif 

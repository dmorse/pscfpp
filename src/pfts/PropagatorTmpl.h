#ifndef PFTS_PROPAGATOR_TMPL_H
#define PFTS_PROPAGATOR_TMPL_H

/*
* PFTS - Polymer Field Theory Simulator
*
* Copyright 2013, David Morse (morse012@.umn.edu)
* Distributed under the terms of the GNU General Public License.
*/

#include <util/containers/GArray.h>
#include <pfts/Block.h>

namespace Pfts{ 

   using namespace Util;

   /**
   * Template for propagator classes.
   *
   * The template argument Propagator should be a concrete class that is
   * derived from the template PropagatorTmpl<Propagator>, with a class
   * definition something like:
   * \code
   *    class Propagator : public PropagatorTmpl<Propagator>
   *    {} 
   * \endcode
   * This usage is an example of the so-called "curiously recurring template 
   * pattern" (CRTP). It is used here to allow the template to have a member
   * array that stores pointers to other instances of class Propagator.
   */
   template <class Propagator>
   class PropagatorTmpl
   {

   public:

      /**
      * Constructor.
      */ 
      PropagatorTmpl();
 
      /**
      * Associate this propagator with a block and direction.
      *
      * \param block associated Block object.
      * \param directionId direction = 0 or 1.
      */ 
      void setBlock(const Block& block, int directionId);

      /**
      * Add another propagator to the list of sources for this one.
      *
      * A source is a propagator that is needed to compute the initial
      * condition for this one, and that thus must be computed before 
      * this one.
      */
      void addSource(const Propagator& other);

      /**
      * Check if all sources are completed.
      */
      bool isReady();
 
      /**
      * Set the isComplete flag to true or false.
      */
      void setIsComplete(bool isComplete);
 
      /**
      * Get the associated Block object by reference.
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
      * 
      * \param id index of source propagator, < nSource
      */
      const Propagator& source(int id) const;

      /**
      * Has the modified diffusion equation been solved?
      */
      bool isComplete() const;
 
   private:
  
      /// Pointer to associated block.
      Block const * blockPtr_;

      /// Direction of propagation (0 or 1)
      int directionId_;

      /// Pointers to propagators that feed source vertex.
      GArray<Propagator const *> sourcePtrs_;

      /// Set true after solving modified diffusion equation.
      bool isComplete_;
  
   };

   // Inline member functions

   /*
   * Get the associated Block object.
   */
   template <class Propagator>
   inline const Block& PropagatorTmpl<Propagator>::block() const
   {  return *blockPtr_; }

   /*
   * Get the direction index.
   */
   template <class Propagator>
   inline int PropagatorTmpl<Propagator>::directionId() const
   {  return directionId_; }

   /*
   * Get the number of source propagators.
   */
   template <class Propagator>
   inline int PropagatorTmpl<Propagator>::nSource() const
   {  return sourcePtrs_.size(); }

   /*
   * Get a source propagator.
   */
   template <class Propagator>
   inline const Propagator& 
   PropagatorTmpl<Propagator>::source(int id) const
   {  return *(sourcePtrs_[id]); }

   /*
   * Is the computation of this propagator completed?
   */
   template <class Propagator>
   inline bool PropagatorTmpl<Propagator>::isComplete() const
   {  return isComplete_; }

   // Noninline member functions

   /*
   * Constructor.
   */
   template <class Propagator>
   PropagatorTmpl<Propagator>::PropagatorTmpl()
    : blockPtr_(0),
      directionId_(-1),
      sourcePtrs_(),
      isComplete_(false)
   {}

   /*
   * Associate this propagator with a block and direction
   */
   template <class Propagator>
   void PropagatorTmpl<Propagator>::setBlock(const Block& block, int directionId)
   {
      blockPtr_ = &block;
      directionId_ = directionId;
   }

   /*
   * Add a source propagator to the list.
   */
   template <class Propagator>
   void PropagatorTmpl<Propagator>::addSource(const Propagator& other)
   {  sourcePtrs_.append(&other); }

   /*
   * Check if all source propagators are marked completed.
   */
   template <class Propagator>
   bool PropagatorTmpl<Propagator>::isReady()
   {
      for (int i=0; i < sourcePtrs_.size(); ++i) {
         if (!sourcePtrs_[i]->isComplete()) {
            return false;
         }
      }
      return true;
   }

   /*
   * Mark this propagator as complete (true) or incomplete (false).
   */
   template <class Propagator>
   void PropagatorTmpl<Propagator>::setIsComplete(bool isComplete)
   {  isComplete_ = isComplete; }

}
#endif 

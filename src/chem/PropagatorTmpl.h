#ifndef CHEM_PROPAGATOR_TMPL_H
#define CHEM_PROPAGATOR_TMPL_H

/*
* PFTS - Polymer Field Theory Simulator
*
* Copyright 2013, David Morse (morse012@.umn.edu)
* Distributed under the terms of the GNU General Public License.
*/

#include <util/containers/GArray.h>
#include <chem/Block.h>

namespace Chem{ 

   using namespace Util;

   /**
   * Template for propagator classes.
   *
   * The template argument Propagator should be a concrete class that is
   * derived from the template PropagatorTmpl<Propagator>, as in the following
   * example: 
   * \code
   *    class Propagator : public PropagatorTmpl<Propagator>
   *    {
   *     ...
   *    }; 
   * \endcode
   * This usage is an example of the so-called "curiously recurring template 
   * pattern" (CRTP). It is used here to allow the template to have a members
   * that stores pointers to other instances of derived class Propagator.
   *
   * The Propagator class is used in templates PolymerTmpl and SystemTmpl that 
   * require that it define the following public typedefs and member functions:
   * \code
   *    class Propagator : public PropagatorTmpl<Propagator> 
   *    {
   *    public:
   * 
   *        // Chemical potential field type.
   *        typedef DArray<double> WField;
   *
   *        // Monomer concentration field type.
   *        typedef DArray<double> CField;
   *
   *        // Solve the modified diffusion equation for this direction.
   *        void solve(const WField& wField);
   *
   *        // Compute the integral \int q(r,s)q^{*}(r,s) for this block
   *        void integrate(const CField& integral);
   * 
   *    }; 
   * \endcode
   * The typedefs WField and CField define the types of the objects used
   * to represent a chemical potential field for a particular monomer type
   * and a monomer concentration field. In the above example, both of these
   * typenames are defined to be synonyms for DArrray<double>, i.e., for 
   * dynamically allocated arrays of double precision floating point 
   * numbers. Other implementations may use more specialized types.
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
      * Set monomer statistical segment length.
      *
      * \param kuhn monomer statistical segment length
      */
      void setKuhn(double kuhn);

      /**
      * Set the partner of this propagator.
      *
      * A partner of a propagator is the propagator for the same block
      * that propagates in the opposite direction.
      *
      * \param partner reference to partner propagator
      */
      void setPartner(const Propagator& partner);

      /**
      * Add a propagator to the list of sources for this one.
      *
      * A source is a propagator that is needed to compute the initial
      * condition for this one, and that thus must be computed before 
      * this one.
      * 
      * \param source reference to source propagator
      */
      void addSource(const Propagator& source);

      /**
      * Check if all sources are completed.
      */
      bool isReady();
 
      /**
      * Set the isSolved flag to true or false.
      */
      void setIsSolved(bool isSolved);
 
      /**
      * Get the associated Block object by reference.
      */
      const Block& block() const;

      /**
      * Get direction index for this propagator.
      */
      int directionId() const;

      /**
      * Get monomer statistical segment length.
      */
      double kuhn() const;

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
      * Does this have a partner propagator?
      */
      bool hasPartner() const;

      /**
      * Get partner propagator.
      */
      const Propagator& partner() const;

      /**
      * Has the modified diffusion equation been solved?
      */
      bool isSolved() const;
 
   private:
  
      /// Monomer statistical segment length.
      double kuhn_;

      /// Pointer to associated block.
      Block const * blockPtr_;

      /// Direction of propagation (0 or 1).
      int directionId_;

      /// Pointer to partner - same block,opposite direction.
      Propagator const * partnerPtr_;

      /// Pointers to propagators that feed source vertex.
      GArray<Propagator const *> sourcePtrs_;

      /// Set true after solving modified diffusion equation.
      bool isSolved_;
  
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
   * Get the monomer statistical segment length. 
   */
   template <class Propagator>
   inline double PropagatorTmpl<Propagator>::kuhn() const
   {  return kuhn_; }

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

   /**
   * Does this have a partner propagator?
   */
   template <class Propagator>
   inline
   bool PropagatorTmpl<Propagator>::hasPartner() const
   {  return partnerPtr_; }

   /*
   * Is the computation of this propagator completed?
   */
   template <class Propagator>
   inline bool PropagatorTmpl<Propagator>::isSolved() const
   {  return isSolved_; }

   // Noninline member functions

   /*
   * Constructor.
   */
   template <class Propagator>
   PropagatorTmpl<Propagator>::PropagatorTmpl()
    : blockPtr_(0),
      directionId_(-1),
      sourcePtrs_(),
      isSolved_(false)
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
   * Set the monomer statistical segment length.
   */
   template <class Propagator>
   void PropagatorTmpl<Propagator>::setKuhn(double kuhn)
   {  kuhn_ = kuhn; }

   /*
   * Set the partner propagator.
   */
   template <class Propagator>
   void PropagatorTmpl<Propagator>::setPartner(const Propagator& partner)
   {  partnerPtr_ = &partner; }

   /*
   * Add a source propagator to the list.
   */
   template <class Propagator>
   void PropagatorTmpl<Propagator>::addSource(const Propagator& source)
   {  sourcePtrs_.append(&source); }

   /*
   * Check if all source propagators are marked completed.
   */
   template <class Propagator>
   bool PropagatorTmpl<Propagator>::isReady()
   {
      for (int i=0; i < sourcePtrs_.size(); ++i) {
         if (!sourcePtrs_[i]->isSolved()) {
            return false;
         }
      }
      return true;
   }

   /*
   * Get partner propagator.
   */
   template <class Propagator>
   const Propagator& PropagatorTmpl<Propagator>::partner() 
   const
   {
      UTIL_CHECK(partnerPtr_);
      return *partnerPtr_;
   }

   /*
   * Mark this propagator as solved (true) or not (false).
   */
   template <class Propagator>
   void PropagatorTmpl<Propagator>::setIsSolved(bool isSolved)
   {  isSolved_ = isSolved; }

}
#endif 

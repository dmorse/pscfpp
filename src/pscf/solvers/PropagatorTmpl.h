#ifndef PSCF_PROPAGATOR_TMPL_H
#define PSCF_PROPAGATOR_TMPL_H

/*
* PSCF - Polymer Self-Consistent Field Theory
*
* Copyright 2016 - 2022, The Regents of the University of Minnesota
* Distributed under the terms of the GNU General Public License.
*/

#include <pscf/chem/PolymerModel.h>
#include <util/containers/GArray.h>

namespace Pscf
{ 

   using namespace Util;

   /**
   * Template for propagator classes.
   *
   * The template argument TP should be a concrete propagator class that 
   * is derived from the template PropagatorTmpl<TP>. By convention, each
   * implementation of SCFT is defined in a different sub-namespace of
   * namespace Pscf. For each such implementation, there is a concrete
   * propagator class, named Propagator by convention, that is a subclass
   * of the template instance PropagatorTmpl<Propagator>, using the syntax 
   * shown below:
   * \code
   *
   *    class Propagator : public PropagatorTmpl<Propagator>
   *    {
   *     ...
   *    }; 
   *
   * \endcode
   * This usage is an example of the so-called "curiously recurring 
   * template pattern" (CRTP). It is used here to allow the template 
   * PropagatorTmpl<Propagator> to have a member variables that store 
   * pointers to other instances of derived class Propagator (or TP).
   *
   * The concrete Propagator class is used in templates BlockTmpl, 
   * PolymerTmpl and SystemTmpl. The usage in those templates require 
   * that this class define typedefs named WField and CField and member
   * functions named solve() and computeQ(), neither of which takes any 
   * arguments.
   * 
   * The typedefs WField and CField must be aliases for a type container 
   * used represent chemical potential and concentration fields, 
   * respectively.
   * 
   * The function void solve() must solve the modified diffusion equation 
   * for this propapagator, using information that is accessible to the
   * concrete propagator class. Upon return, the solution must be stored
   * in internal data structures that store the propagator field at 
   * different values of a contour variable, and the function isSolved() 
   * must return true.
   *
   * The function computeQ() must compute and return the molecular 
   * partition function Q using information available to this propagator. 
   * An example of the required interface is shown below:
   * 
   * \code
   *
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
   *        void solve();
   *
   *        // Compute and return the molecular partition function Q.
   *        double computeQ();
   *
   *    };
   *
   * \endcode
   *
   * In the above example, the field container typenames WField and 
   * CField are defined to be synonyms for DArrray<double>, i.e., for 
   * dynamically allocated arrays of double precision floating point 
   * numbers. Other implementations may use more specialized types.
   *
   * \ingroup Pscf_Solver_Module
   */
   template <class TP>
   class PropagatorTmpl
   {

   public:

      /**
      * Constructor.
      */ 
      PropagatorTmpl();

      /// \name Mutators
      ///@{
 
      /**
      * Associate this propagator with a direction index.
      *
      * \param directionId direction = 0 or 1.
      */ 
      void setDirectionId(int directionId);

      /**
      * Set the partner of this propagator.
      *
      * The partner of a propagator is the propagator for the same block
      * that propagates in the opposite direction.
      *
      * \param partner reference to partner propagator
      */
      void setPartner(const TP& partner);

      /**
      * Add a propagator to the list of sources for this one.
      *
      * A source is a propagator that terminates at the root vertex 
      * of this one and is needed to compute the initial condition 
      * for this one, and that thus must be computed before this.
      * 
      * \param source reference to source propagator
      */
      void addSource(const TP& source);

      /**
      * Set the isSolved flag to true or false.
      */
      void setIsSolved(bool isSolved);

      /**
      * Set whether the propagator owns the head and tail (bead model).
      *
      * The concept of "ownership" of an attached vertex is meaningful
      * in a bead-spring model, in which we require that each vertex
      * bead be treated as part of exactly one of the attached blocks.
      * Because it is meaningless in the context of a thread model, it
      * is an error to call this function when PolymerModel::isThread().
      *
      * Precondition: PolymerModel::isBead() must be true.
      *
      * \param ownsHead  Does this propagator own the head vertex bead?
      * \param ownsTail  Does this propagator own the tail vertex bead?
      */
      void setVertexOwnership(bool ownsHead, bool ownsTail);
 
      ///@}
      /// \name Accessors
      ///@{

      /**
      * Get a source propagator.
      * 
      * \param id index of source propagator, < nSource
      */
      const TP& source(int id) const;

      /**
      * Get partner propagator.
      */
      const TP& partner() const;

      /**
      * Get direction index for this propagator.
      */
      int directionId() const;

      /**
      * Number of source / prerequisite propagators.
      */
      int nSource() const;

      /**
      * Does this have a partner propagator?
      */
      bool hasPartner() const;

      /**
      * Does this propagator own the attached head vertex bead?
      *
      * Precondition: PolymerModel::isBead()
      */
      bool ownsHead() const;

      /**
      * Does this propagator own the attached tail vertex bead?
      *
      * Precondition: PolymerModel::isBead()
      */
      bool ownsTail() const;

      /**
      * Has the modified diffusion equation been solved?
      */
      bool isSolved() const;
 
      /**
      * Are all source propagators solved?
      */
      bool isReady() const;
 
      ///@}

   private:
  
      /// Direction of propagation (0 or 1).
      int directionId_;

      /// Pointer to partner - same block,opposite direction.
      TP const * partnerPtr_;

      /// Pointers to propagators that feed source vertex.
      GArray<TP const *> sourcePtrs_;

      /// Does this propagator own the head vertex bead (bead model)?
      /// Only used or meaningful for the bead model
      bool ownsHead_;

      /// Does this propagator own the tail vertex bead (bead model)?
      /// Only used or meaningful for the bead model
      bool ownsTail_;

      /// Set true after solving modified diffusion equation.
      bool isSolved_;
  
   };

   // Inline member functions

   /*
   * Get the direction index.
   */
   template <class TP>
   inline int PropagatorTmpl<TP>::directionId() const
   {  return directionId_; }

   /*
   * Get the number of source propagators.
   */
   template <class TP>
   inline int PropagatorTmpl<TP>::nSource() const
   {  return sourcePtrs_.size(); }

   /*
   * Get a source propagator.
   */
   template <class TP>
   inline const TP& 
   PropagatorTmpl<TP>::source(int id) const
   {  return *(sourcePtrs_[id]); }

   /**
   * Does this have a partner propagator?
   */
   template <class TP>
   inline
   bool PropagatorTmpl<TP>::hasPartner() const
   {  return partnerPtr_; }

   /*
   * Is the computation of this propagator completed?
   */
   template <class TP>
   inline bool PropagatorTmpl<TP>::isSolved() const
   {  return isSolved_; }

   // Noninline member functions

   /*
   * Constructor.
   */
   template <class TP>
   PropagatorTmpl<TP>::PropagatorTmpl()
    : directionId_(-1),
      partnerPtr_(0),
      sourcePtrs_(),
      isSolved_(false)
   {}

   /*
   * Set the directionId.
   */
   template <class TP>
   void PropagatorTmpl<TP>::setDirectionId(int directionId)
   {  directionId_ = directionId; }

   /*
   * Set the partner propagator.
   */
   template <class TP>
   void PropagatorTmpl<TP>::setPartner(const TP& partner)
   {  partnerPtr_ = &partner; }

   /*
   * Add a source propagator to the list.
   */
   template <class TP>
   void PropagatorTmpl<TP>::addSource(const TP& source)
   {  sourcePtrs_.append(&source); }

   /*
   * Set vertex ownership.
   */
   template <class TP>
   void PropagatorTmpl<TP>::setVertexOwnership(bool ownsHead, bool ownsTail)
   {
      UTIL_CHECK(PolymerModel::isBead());
      ownsHead_ = ownsHead;  
      ownsTail_ = ownsTail;  
   }

   // Accessors

   /*
   * Get partner propagator.
   */
   template <class TP>
   const TP& PropagatorTmpl<TP>::partner() 
   const
   {
      UTIL_CHECK(partnerPtr_);
      return *partnerPtr_;
   }

   /*
   * Mark this propagator as solved (true) or not (false).
   */
   template <class TP>
   void PropagatorTmpl<TP>::setIsSolved(bool isSolved)
   {  isSolved_ = isSolved; }

   /*
   * Check if all source propagators are marked completed.
   */
   template <class TP>
   bool PropagatorTmpl<TP>::isReady() const
   {
      for (int i=0; i < sourcePtrs_.size(); ++i) {
         if (!sourcePtrs_[i]->isSolved()) {
            return false;
         }
      }
      return true;
   }

   /*
   * Does this propagator own the head vertex bead?
   */
   template <class TP>
   inline bool PropagatorTmpl<TP>::ownsHead() const
   {
      UTIL_CHECK(PolymerModel::isBead());
      return ownsHead_; 
   }

   /*
   * Does this propagator own the tail vertex bead?
   */
   template <class TP>
   inline bool PropagatorTmpl<TP>::ownsTail() const
   {
      UTIL_CHECK(PolymerModel::isBead());
      return ownsTail_; 
   }

}
#endif 

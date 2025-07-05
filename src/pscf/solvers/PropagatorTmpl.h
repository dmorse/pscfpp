#ifndef PSCF_PROPAGATOR_TMPL_H
#define PSCF_PROPAGATOR_TMPL_H

/*
* PSCF - Polymer Self-Consistent Field Theory
*
* Copyright 2015 - 2025, The Regents of the University of Minnesota
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
   * The template argument QT should be a concrete propagator class that 
   * is derived from the template PropagatorTmpl<QT>. By convention, each
   * implementation of field theory is defined in a different sub-namespace
   * of namespace Pscf. For each such implementation, there is a concrete
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
   * template pattern" (CRQT). It is used here to allow the template 
   * PropagatorTmpl<Propagator> to have member variables that store 
   * pointers to other instances of the derived class Propagator (or QT).
   *
   * The concrete Propagator class is used in templates BlockTmpl, 
   * PolymerTmpl and MixtureTmpl. The usage in those templates require that
   * the Propagator class define an alias named FieldT for the field type.
   * Propagator must also provide member functions named solve() and 
   * computeQ(), neither of which takes any arguments.
   * 
   * The type FieldT must be an aliases for the type of containers used 
   * to store a single chemical potential or concentration field, or a
   * slice of a propagator.
   * 
   * The function void solve() must solve the modified diffusion equation 
   * for this propapagator, using information that is accessible to the
   * concrete propagator class. Upon return, the solution must be stored
   * in an internal data structures that is read-accessible through the 
   * public class interface, and the function isSolved() must return true.
   *
   * The function computeQ() must compute and return the molecular 
   * partition function Q using information available to this propagator. 
   *
   * A simple example of the required interface is shown below:
   * \code
   *
   *    class Propagator : public PropagatorTmpl<Propagator> 
   *    {
   *    public:
   * 
   *        // Field type.
   *        using FieldT = DArray<double>;
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
   * In the above example, the field container typenames FieldT is
   * an alias for DArrray<double>, i.e., for a dynamically allocated 
   * arrays of double precision floating point numbers. Other 
   * implementations may use more specialized types.
   *
   * \ingroup Pscf_Solver_Module
   */
   template <class QT>
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
      void setPartner(const QT& partner);

      /**
      * Add a propagator to the list of sources for this one.
      *
      * A source is a propagator that terminates at the root vertex 
      * of this one and is needed to compute the initial condition 
      * for this one, and that thus must be computed before this.
      * 
      * \param source reference to source propagator
      */
      void addSource(const QT& source);

      /**
      * Set the isSolved flag to true or false.
      */
      void setIsSolved(bool isSolved);

      /**
      * Set flags indicating whether vertices are chain ends.
      *
      * \param isHeadEnd  Does this propagator own the head vertex bead?
      * \param isTailEnd  Does this propagator own the tail vertex bead?
      */
      void setEndFlags(bool isHeadEnd, bool isTailEnd);
 
      ///@}
      /// \name Accessors
      ///@{

      /**
      * Get a source propagator.
      * 
      * \param id index of source propagator, < nSource
      */
      const QT& source(int id) const;

      /**
      * Get partner propagator.
      */
      const QT& partner() const;

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
      * Is the head vertex a chain end?
      *
      * Precondition: PolymerModel::isBead()
      */
      bool isHeadEnd() const;

      /**
      * Is the tail vertex a chain end?
      *
      * Precondition: PolymerModel::isBead()
      */
      bool isTailEnd() const;

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
      QT const * partnerPtr_;

      /// Pointers to propagators that feed source vertex.
      GArray<QT const *> sourcePtrs_;

      /// True iff the head vertex is chain end.
      bool isHeadEnd_;

      /// True iff the tail vertex is chain end.
      bool isTailEnd_;

      /// Set true after solving modified diffusion equation.
      bool isSolved_;
  
   };

   // Inline member functions

   /*
   * Get the direction index.
   */
   template <class QT>
   inline int PropagatorTmpl<QT>::directionId() const
   {  return directionId_; }

   /*
   * Get the number of source propagators.
   */
   template <class QT>
   inline int PropagatorTmpl<QT>::nSource() const
   {  return sourcePtrs_.size(); }

   /*
   * Get a source propagator.
   */
   template <class QT>
   inline const QT& 
   PropagatorTmpl<QT>::source(int id) const
   {  return *(sourcePtrs_[id]); }

   /**
   * Does this have a partner propagator?
   */
   template <class QT>
   inline bool PropagatorTmpl<QT>::hasPartner() const
   {  return partnerPtr_; }

   /*
   * Does this propagator own the head vertex bead?
   */
   template <class QT>
   inline bool PropagatorTmpl<QT>::isHeadEnd() const
   {  return isHeadEnd_; }

   /*
   * Does this propagator own the tail vertex bead?
   */
   template <class QT>
   inline bool PropagatorTmpl<QT>::isTailEnd() const
   {  return isTailEnd_; }

   /*
   * Is the computation of this propagator completed?
   */
   template <class QT>
   inline bool PropagatorTmpl<QT>::isSolved() const
   {  return isSolved_; }

   // Noninline member functions

   /*
   * Constructor.
   */
   template <class QT>
   PropagatorTmpl<QT>::PropagatorTmpl()
    : directionId_(-1),
      partnerPtr_(0),
      sourcePtrs_(),
      isSolved_(false)
   {}

   /*
   * Set the directionId.
   */
   template <class QT>
   void PropagatorTmpl<QT>::setDirectionId(int directionId)
   {  directionId_ = directionId; }

   /*
   * Set the partner propagator.
   */
   template <class QT>
   void PropagatorTmpl<QT>::setPartner(const QT& partner)
   {  partnerPtr_ = &partner; }

   /*
   * Add a source propagator to the list.
   */
   template <class QT>
   void PropagatorTmpl<QT>::addSource(const QT& source)
   {  sourcePtrs_.append(&source); }

   /*
   * Set flags indicate whether vertices are are chain ends.
   */
   template <class QT>
   void PropagatorTmpl<QT>::setEndFlags(bool isHeadEnd, bool isTailEnd)
   {
      isHeadEnd_ = isHeadEnd;  
      isTailEnd_ = isTailEnd;  
   }

   /*
   * Check if all source propagators are marked completed.
   */
   template <class QT>
   bool PropagatorTmpl<QT>::isReady() const
   {
      for (int i=0; i < sourcePtrs_.size(); ++i) {
         if (!sourcePtrs_[i]->isSolved()) {
            return false;
         }
      }
      return true;
   }

   /*
   * Mark this propagator as solved (true) or not (false).
   */
   template <class QT>
   void PropagatorTmpl<QT>::setIsSolved(bool isSolved)
   {  isSolved_ = isSolved; }

   /*
   * Get partner propagator.
   */
   template <class QT>
   const QT& PropagatorTmpl<QT>::partner() 
   const
   {
      UTIL_CHECK(partnerPtr_);
      return *partnerPtr_;
   }

}
#endif 

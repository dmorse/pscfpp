#ifndef PSCF_PROPAGATOR_TMPL_H
#define PSCF_PROPAGATOR_TMPL_H

/*
* PSCF - Polymer Self-Consistent Field Theory
*
* Copyright 2016 - 2022, The Regents of the University of Minnesota
* Distributed under the terms of the GNU General Public License.
*/

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
   * The TP propagator class is used in templates BlockTmpl, PolymerTmpl 
   * and SystemTmpl. The usage in those templates require that it define 
   * the following public typedefs and member functions:
   * \code
   *
   *    class TP : public PropagatorTmpl<TP> 
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
   * The typedefs WField and CField define the types of the objects used 
   * to represent a chemical potential field for a particular monomer type
   * and a monomer concentration field. In the above example, both of these
   * typenames are defined to be synonyms for DArrray<double>, i.e., for 
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
      //@{
 
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
 
      //@}
      /// \name Accessors
      //@{

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
      * Has the modified diffusion equation been solved?
      */
      bool isSolved() const;
 
      /**
      * Are all source propagators are solved?
      */
      bool isReady() const;
 
      //@}

   private:
  
      /// Direction of propagation (0 or 1).
      int directionId_;

      /// Pointer to partner - same block,opposite direction.
      TP const * partnerPtr_;

      /// Pointers to propagators that feed source vertex.
      GArray<TP const *> sourcePtrs_;

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

}
#endif 

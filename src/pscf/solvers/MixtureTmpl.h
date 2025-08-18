#ifndef PSCF_MIXTURE_TMPL_H
#define PSCF_MIXTURE_TMPL_H

/*
* PSCF - Polymer Self-Consistent Field
*
* Copyright 2015 - 2025, The Regents of the University of Minnesota
* Distributed under the terms of the GNU General Public License.
*/

#include <pscf/chem/MixtureBase.h>
#include <util/param/ParamComposite.h>
#include <util/containers/DArray.h>

namespace Pscf
{

   using namespace Util;

   /**
   * Solvers for a mixture of polymer and solvent species.
   *
   * \ingroup Pscf_Solver_Module
   */
   template <class PT, class ST>
   class MixtureTmpl : public MixtureBase, public ParamComposite
   {
   public:

      // Public typedefs

      /**
      * Solvent species solver type.
      */
      using SolventT = ST;

      /**
      * Polymer species solver type.
      */
      using PolymerT = PT;

      /**
      * Block polymer block type.
      */
      using BlockT = typename PT::BlockT;

      /**
      * Polymer block propagator type.
      */
      using PropagatorT = typename BlockT::PropagatorT;

      // Public member functions

      /**
      * Constructor.
      */
      MixtureTmpl();

      /**
      * Destructor.
      */
      ~MixtureTmpl();

      /**
      * Read parameters from file and initialize.
      *
      * \param in  input parameter file
      */
      virtual void readParameters(std::istream& in);

      /// \name Accessors 
      ///@{

      /*
      * Public member functions inherited from MixtureBase:
      *
      *   int nMonomer() const;
      *   int nPolymer() const;
      *   int nSolvent() const;
      *   int nBlock() const;
      *   double vMonomer() const;
      *   Monomer const & monomer(int id) const;
      *
      *   void setVmonomer(double);
      */

      /**
      * Get a polymer solver object by non-const reference.
      *
      * \param id  integer polymer species index (0 <= id < nPolymer)
      */
      PolymerT& polymer(int id);

      /**
      * Get a polymer solver by const reference.
      *
      * \param id  integer polymer species index (0 <= id < nPolymer)
      */
      PolymerT const & polymer(int id) const;

      /**
      * Get a PolymerSpecies descriptor by const reference.
      *
      * Defines function declared pure virtual by MixtureBase.
      *
      * \param id  integer polymer species index (0 <= id < nPolymer)
      */
      PolymerSpecies const & polymerSpecies(int id) const final;

      /**
      * Get a solvent solver object.
      *
      * \param id  integer solvent species index (0 <= id < nSolvent)
      */
      SolventT& solvent(int id);

      /**
      * Get a solvent solver object by constant reference.
      *
      * \param id  integer solvent species index (0 <= id < nSolvent)
      */
      SolventT const & solvent(int id) const;

      /**
      * Set a SolventSpecies descriptor object by const reference.
      *
      * Defines function declared pure virtual by MixtureBase.
      *
      * \param id integer solvent species index (0 <= id < nSolvent)
      */
      SolventSpecies const & solventSpecies(int id) const final;

      ///@}

   private:

      /**
      * Array of polymer species solver objects.
      *
      * Array capacity = nPolymer.
      */
      DArray<PolymerT> polymers_;

      /**
      * Array of solvent species objects.
      *
      * Array capacity = nSolvent.
      */
      DArray<SolventT> solvents_;

      // Restrict access to inherited protected data
      using MixtureBase::monomers_;
      using MixtureBase::nMonomer_;
      using MixtureBase::nPolymer_;
      using MixtureBase::nSolvent_;
      using MixtureBase::nBlock_;
      using MixtureBase::vMonomer_;

   };

   // Inline member functions

   template <class PT, class ST>
   inline PT& MixtureTmpl<PT,ST>::polymer(int id)
   {  
      UTIL_CHECK(id < nPolymer_);
      return polymers_[id];
   }

   template <class PT, class ST>
   inline PT const & MixtureTmpl<PT,ST>::polymer(int id) const
   {  
      UTIL_CHECK(id < nPolymer_);
      return polymers_[id];
   }

   template <class PT, class ST>
   inline ST& MixtureTmpl<PT,ST>::solvent(int id)
   {  
      UTIL_CHECK(id < nSolvent_);
      return solvents_[id]; 
   }

   template <class PT, class ST>
   inline ST const & MixtureTmpl<PT,ST>::solvent(int id) const
   {  
      UTIL_CHECK(id < nSolvent_);
      return solvents_[id]; 
   }

}
#include "MixtureTmpl.tpp"
#endif

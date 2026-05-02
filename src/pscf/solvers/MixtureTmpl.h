#ifndef PSCF_MIXTURE_TMPL_H
#define PSCF_MIXTURE_TMPL_H

/*
* PSCF - Polymer Self-Consistent Field
*
* Copyright 2015 - 2025, The Regents of the University of Minnesota
* Distributed under the terms of the GNU General Public License.
*/

#include <pscf/chem/MixtureBase.h>       // base class
#include <util/param/ParamComposite.h>   // base class
#include <util/containers/DArray.h>      // member template

namespace Pscf
{

   using namespace Util;

   /**
   * Solvers for a mixture of polymer and solvent species.
   *
   * \ingroup Pscf_Solver_Module
   */
   template <class PT, class ST, typename WT = double>
   class MixtureTmpl : public MixtureBase<WT>, public ParamComposite
   {
   public:

      // Public type aliases

      /**
      * Solvent species solver type.
      */
      using SolventT = ST;

      /**
      * Polymer species solver type.
      */
      using PolymerT = PT;

      // Public member functions

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
      * Read parameters from file and initialize.
      *
      * \param in  input parameter file
      */
      void readParameters(std::istream& in) override;

      /// \name Accessors 
      ///@{

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
      PolymerSpecies<WT> const & polymerSpecies(int id) const final;

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
      SolventSpecies<WT> const & solventSpecies(int id) const final;

      ///@}

   protected:

      /**
      * Constructor.
      */
      MixtureTmpl();

      /**
      * Destructor.
      */
      ~MixtureTmpl() override = default;

      /**
      * Alias for base class.
      */
      using  MixtureBaseT = MixtureBase<WT>;

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
      using MixtureBase<WT>::monomers_;
      using MixtureBase<WT>::nMonomer_;
      using MixtureBase<WT>::nPolymer_;
      using MixtureBase<WT>::nSolvent_;
      using MixtureBase<WT>::nBlock_;
      using MixtureBase<WT>::vMonomer_;

   };

   // Inline member functions

   template <class PT, class ST, typename WT>
   inline PT& MixtureTmpl<PT,ST,WT>::polymer(int id)
   {  
      UTIL_CHECK(id < nPolymer_);
      return polymers_[id];
   }

   template <class PT, class ST, typename WT>
   inline PT const & MixtureTmpl<PT,ST,WT>::polymer(int id) const
   {  
      UTIL_CHECK(id < nPolymer_);
      return polymers_[id];
   }

   template <class PT, class ST, typename WT>
   inline ST& MixtureTmpl<PT,ST,WT>::solvent(int id)
   {  
      UTIL_CHECK(id < nSolvent_);
      return solvents_[id]; 
   }

   template <class PT, class ST, typename WT>
   inline ST const & MixtureTmpl<PT,ST,WT>::solvent(int id) const
   {  
      UTIL_CHECK(id < nSolvent_);
      return solvents_[id]; 
   }

}
#endif

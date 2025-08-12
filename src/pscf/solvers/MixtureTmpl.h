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

   // Non-inline member functions

   /*
   * Constructor.
   */
   template <class PT, class ST>
   MixtureTmpl<PT,ST>::MixtureTmpl()
    : MixtureBase(),
      ParamComposite(),
      polymers_(),
      solvents_()
   {}

   /*
   * Destructor.
   */
   template <class PT, class ST>
   MixtureTmpl<PT,ST>::~MixtureTmpl()
   {}

   /*
   * Get a PolymerSpecies descriptor by non-const reference.
   */
   template <class PT, class ST>
   PolymerSpecies const & MixtureTmpl<PT,ST>::polymerSpecies(int id) const
   {  
      UTIL_CHECK(id < nPolymer_);
      return polymers_[id];
   }

   /*
   * Get a SolventSpecies descriptor by const reference.
   */
   template <class PT, class ST>
   SolventSpecies const & MixtureTmpl<PT,ST>::solventSpecies(int id) const
   {  
      UTIL_CHECK(id < nSolvent_);
      return solvents_[id]; 
   }

   /*
   * Read all parameters and initialize.
   */
   template <class PT, class ST>
   void MixtureTmpl<PT,ST>::readParameters(std::istream& in)
   {
      // Read nMonomer and monomers array
      read<int>(in, "nMonomer", nMonomer_);
      monomers_.allocate(nMonomer_);
      for (int i = 0; i < nMonomer_; ++i) {
         monomers_[i].setId(i);
      }
      readDArray< Monomer >(in, "monomers", monomers_, nMonomer_);

      /*
      * The input format for a single monomer is defined in the istream
      * extraction operation (operator >>) for a Pscf::Monomer, in file
      * pscf/chem/Monomer.cpp. The text representation contains only the
      * value for the monomer statistical segment Monomer::kuhn.
      */

      // Read nPolymer
      read<int>(in, "nPolymer", nPolymer_);
      UTIL_CHECK(nPolymer_ > 0);

      // Optionally read nSolvent, with nSolvent=0 by default
      nSolvent_ = 0;
      readOptional<int>(in, "nSolvent", nSolvent_);

      // Read polymers and compute nBlock
      nBlock_ = 0;
      if (nPolymer_ > 0) {

         polymers_.allocate(nPolymer_);
         for (int i = 0; i < nPolymer_; ++i) {
            readParamComposite(in, polymer(i));
            nBlock_ = nBlock_ + polymer(i).nBlock();
         }
   
         // Set statistical segment lengths for all blocks
         double kuhn;
         int monomerId;
         for (int i = 0; i < nPolymer_; ++i) {
            for (int j = 0; j < polymer(i).nBlock(); ++j) {
               monomerId = polymer(i).block(j).monomerId();
               kuhn = monomer(monomerId).kuhn();
               polymer(i).block(j).setKuhn(kuhn);
            }
         }

      }

      // Read solvents
      if (nSolvent_ > 0) {

         solvents_.allocate(nSolvent_);
         for (int i = 0; i < nSolvent_; ++i) {
            readParamComposite(in, solvent(i));
         }

      }

      // Optionally read monomer reference value
      vMonomer_ = 1.0; // Default value
      readOptional(in, "vMonomer", vMonomer_);

   }

}
#endif

#ifndef PSCF_MIXTURE_TMPL_H
#define PSCF_MIXTURE_TMPL_H

/*
* PSCF - Polymer Self-Consistent Field Theory
*
* Copyright 2016 - 2022, The Regents of the University of Minnesota
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
   template <class TP, class TS>
   class MixtureTmpl : public MixtureBase, public ParamComposite
   {
   public:

      // Public typedefs

      /**
      * Polymer species solver typename.
      */
      typedef TP Polymer;

      /**
      * Solvent species solver typename.
      */
      typedef TS Solvent;

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
      * Get a Polymer solver object.
      *
      * \param id  integer polymer species index (0 <= id < nPolymer)
      */
      Polymer& polymer(int id);

      /**
      * Get a Polymer solver by const reference.
      *
      * \param id  integer polymer species index (0 <= id < nPolymer)
      */
      Polymer const & polymer(int id) const;

      /**
      * Get a PolymerSpecies descriptor by const reference.
      *
      * Implements pure virtual function
      *
      * \param id  integer polymer species index (0 <= id < nPolymer)
      */
      PolymerSpecies const & polymerSpecies(int id) const final;

      /**
      * Set a Solvent solver object.
      *
      * \param id  integer solvent species index (0 <= id < nSolvent)
      */
      Solvent& solvent(int id);

      /**
      * Set a Solvent solver object by constant reference.
      *
      * \param id  integer solvent species index (0 <= id < nSolvent)
      */
      Solvent const & solvent(int id) const;

      /**
      * Set a SolventSpecies descriptor object by const reference.
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
      DArray<Polymer> polymers_;

      /**
      * Array of solvent species objects.
      *
      * Array capacity = nSolvent.
      */
      DArray<Solvent> solvents_;

      // Change inherited data member access from protected to private
      using MixtureBase::monomers_;
      using MixtureBase::nMonomer_;
      using MixtureBase::nPolymer_;
      using MixtureBase::nSolvent_;
      using MixtureBase::nBlock_;
      using MixtureBase::vMonomer_;

   };

   // Inline member functions

   template <class TP, class TS>
   inline TP& MixtureTmpl<TP,TS>::polymer(int id)
   {  
      UTIL_CHECK(id < nPolymer_);
      return polymers_[id];
   }

   template <class TP, class TS>
   inline TP const & MixtureTmpl<TP,TS>::polymer(int id) const
   {  
      UTIL_CHECK(id < nPolymer_);
      return polymers_[id];
   }

   template <class TP, class TS>
   inline TS& MixtureTmpl<TP,TS>::solvent(int id)
   {  
      UTIL_CHECK(id < nSolvent_);
      return solvents_[id]; 
   }

   template <class TP, class TS>
   inline TS const & MixtureTmpl<TP,TS>::solvent(int id) const
   {  
      UTIL_CHECK(id < nSolvent_);
      return solvents_[id]; 
   }

   // Non-inline member functions

   /*
   * Constructor.
   */
   template <class TP, class TS>
   MixtureTmpl<TP,TS>::MixtureTmpl()
    : MixtureBase(),
      ParamComposite(),
      polymers_(),
      solvents_()
   {}

   /*
   * Destructor.
   */
   template <class TP, class TS>
   MixtureTmpl<TP,TS>::~MixtureTmpl()
   {}

   template <class TP, class TS>
   PolymerSpecies const & MixtureTmpl<TP,TS>::polymerSpecies(int id) const
   {  
      UTIL_CHECK(id < nPolymer_);
      return polymers_[id];
   }

   template <class TP, class TS>
   SolventSpecies const & MixtureTmpl<TP,TS>::solventSpecies(int id) const
   {  
      UTIL_CHECK(id < nSolvent_);
      return solvents_[id]; 
   }

   /*
   * Read all parameters and initialize.
   */
   template <class TP, class TS>
   void MixtureTmpl<TP,TS>::readParameters(std::istream& in)
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

#ifndef PFTS_SYSTEM_TMPL_H
#define PFTS_SYSTEM_TMPL_H

/*
* PFTS - Polymer Field Theory Simulator
*
* Copyright 2013, David Morse (morse012@.umn.edu)
* Distributed under the terms of the GNU General Public License.
*/

#include <pfts/chem/Monomer.h>
#include <util/param/ParamComposite.h>
#include <util/containers/DArray.h>

namespace Pfts{

   using namespace Util;

   template <class TPolymer, class TSolvent,
             class TWField, class TCField>
   class SystemTmpl : public ParamComposite
   {
   public:

      /**
      * Read parameters from file and initialize.
      */
      virtual void readParameters(std::istream& in);

      /**
      * Compute ideal gas properties for all species.
      */
      virtual void compute()
      {}

      /**
      * Get a Monomer type descriptor.
      *
      * \param id integer monomer type index
      */
      Monomer& monomer(int id);

      /**
      * Get a polymer object.
      *
      * \param id integer polymer species index
      */
      TPolymer& polymer(int id);

      /**
      * Set a solvent solver object.
      *
      * \param id integer solvent species index
      */
      TSolvent& solvent(int id);

      /**
      * Return W (chemical potential) field for a specific monomer type.
      *
      * \param monomerId integer monomer type index
      */
      TWField& wField(int monomerId);

      /**
      * Return concentration field for specific monomer type.
      *
      * \param monomerId integer monomer type index
      */
      TCField& cField(int id);

      /**
      * Get number of monomers.
      */
      int nMonomer() const;

      /**
      * Get number of polymer species.
      */
      int nPolymer() const;

      /**
      * Get number of solvent (point particle) species.
      */
      int nSolvent() const;

   private:

      /**
      * Array of monomer type descriptors.
      */
      DArray<Monomer> monomers_;

      /**
      * Array of chemical potential fields for monomer types.
      *
      * Indexed by monomer typeId, size = nMonomer.
      */
      DArray<TWField> wFields_;

      /**
      * Array of concentration fields for monomer types.
      *
      * Indexed by monomer typeId, size = nMonomer.
      */
      DArray<TCField> cFields_;

      /**
      * Array of polymer species solvers objects.
      *
      * Size = nPolymer.
      */
      DArray<TPolymer> polymers_;

      /**
      * Array of solvent species solvers.
      *
      * Size = nSolvent.
      */
      DArray<TSolvent> solvents_;

      int nMonomer_; 

      int nPolymer_;

      int nSolvent_;

   };

   template <class TP, class TS, class TW, class TC>
   inline int SystemTmpl<TP,TS,TW,TC>::nMonomer() const
   {  return nMonomer_; }

   template <class TP, class TS, class TW, class TC>
   inline int SystemTmpl<TP,TS,TW,TC>::nPolymer() const
   {  return nPolymer_; }

   template <class TP, class TS, class TW, class TC>
   inline int SystemTmpl<TP,TS,TW,TC>::nSolvent() const
   {  return nSolvent_; }

   template <class TP, class TS, class TW, class TC>
   inline Monomer& SystemTmpl<TP,TS,TW,TC>::monomer(int id)
   {  return monomers_[id]; }

   template <class TP, class TS, class TW, class TC>
   inline TP& SystemTmpl<TP,TS,TW,TC>::polymer(int id)
   {  return polymers_[id]; }

   template <class TP, class TS, class TW, class TC>
   inline TS& SystemTmpl<TP,TS,TW,TC>::solvent(int id)
   {  return solvents_[id]; }

   template <class TP, class TS, class TW, class TC>
   inline TW& SystemTmpl<TP,TS,TW,TC>::wField(int id)
   {  return wFields_[id]; }

   template <class TP, class TS, class TW, class TC>
   inline TC& SystemTmpl<TP,TS,TW,TC>::cField(int id)
   {  return cFields_[id]; }

   template <class TP, class TS, class TW, class TC>
   void SystemTmpl<TP,TS,TW,TC>::readParameters(std::istream& in)
   {
      read<int>(in, "nMonomer", nMonomer_);
      monomers_.allocate(nMonomer_);
      readDArray< Monomer >(in, "monomers", monomers_, nMonomer_);
      read<int>(in, "nPolymer", nPolymer_);
      polymers_.allocate(nPolymer_);
      readBegin(in, "Polymers");
      for (int i = 0; i < nPolymer_; ++i) {
         polymers_[i].readParameters(in);
      }
      readEnd(in);
   }

}
#endif

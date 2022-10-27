#ifndef PSPC_MC_MOVE_MANAGER_H
#define PSPC_MC_MOVE_MANAGER_H

/*
* Simpatico - Simulation Package for Polymeric and Molecular Liquids
*
* Copyright 2010 - 2017, The Regents of the University of Minnesota
* Distributed under the terms of the GNU General Public License.
*/

#include "McMove.h"                      // base class template parameter
#include <util/param/Manager.h>          // base class template
#include <util/containers/DArray.h>      // member template

namespace Util { class Random; }

namespace Pscf {
namespace Pspc {

   using namespace Util;

   template <int D> class System;

   /**
   * Manager for a set of McMove objects.
   *
   * \sa \ref mcMd_mcMove_McMoveManager_page "parameter file format"
   *
   * \ingroup Pspc_Manager_Module
   * \ingroup Pspc_McMove_Module
   */
   template <int D>
   class McMoveManager : public Manager< McMove<D> >
   {

   public:

      /**
      * Constructor.
      *
      * \param system parent System
      */
      McMoveManager(System<D>& system);

      /**
      * Destructor.
      */
      ~McMoveManager();

      /**
      * Read instructions for creating McMove objects.
      *
      * \param in input parameter stream
      */
      virtual void readParameters(std::istream &in);

      /**
      * Initialize at beginning of system run.
      *
      * This method calls the initialize method for every McMove.
      */
      void setup();

      /**
      * Choose an McMove at random, using specified probabilities.
      *
      * \return chosen McMove
      */
      McMove<D>& chooseMove();

      /**
      * Return probability of move i.
      *
      * \param i index for McMove
      * \return probability of McMove number i
      */
      double probability(int i) const;

      /**
      * Output statistics for all moves.
      */
      void output();

      using Manager< McMove<D> >::size;

   protected:

      using Manager< McMove<D> >::setClassName;
       
   private:

      /// Array of McMove probabilities.
      DArray<double>  probabilities_;

      /// Pointer to parent System.
      System<D>* systemPtr_;

      /// Pointer to random number generator.
      Random* randomPtr_;

      /// Return pointer to a new McMoveFactory.
      virtual Factory< McMove<D> >* newDefaultFactory() const;

   };

   // Inline functions

   /*
   * Return probability of move number i
   */
   template <int D>
   inline double McMoveManager<D>::probability(int i) const
   {
      assert(i >= 0);  
      assert(i < size());  
      return probabilities_[i];
   }

}
}
#endif

#ifndef PSPC_MC_MOVE_MANAGER_H
#define PSPC_MC_MOVE_MANAGER_H

/*
* PSCF - Polymer Self-Consistent Field Theory
*
* Copyright 2016 - 2022, The Regents of the University of Minnesota
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
   template <int D> class McSimulator;

   /**
   * Manager for a set of McMove objects.
   *
   * \sa \ref mcMd_mcMove_McSimulator_page "parameter file format"
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
      * \param simulation parent MCsimulator
      */
      McMoveManager(McSimulator<D>& simulator, System<D>& system);

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
      * Output statistics for all moves.
      */
      void output();

      /**
      * Return probability of move i.
      *
      * \param i index for McMove
      * \return probability of McMove number i
      */
      double probability(int i) const;
      
      using Manager< McMove<D> >::size;
  
   protected:
      
      using Manager< McMove<D> >::setClassName;
   
   private:

      // Private data members

      /**
      * Array of McMove probabilities.
      */
      DArray<double>  probabilities_;
      
      /**
       * Pointer to parent Simulator
       */
      McSimulator<D>* mcSimulatorPtr_;

      /**
      * Pointer to parent System.
      */
      System<D>* systemPtr_;

      /**
      * Pointer to random number generator.
      */
      Random* randomPtr_;

      // Private member functions

      /**
      * Return pointer to a new McMoveFactory.
      */
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

   #ifndef PSPC_MC_MOVE_MANAGER_TPP
   // Suppress implicit instantiation
   extern template class McMoveManager<1>;
   extern template class McMoveManager<2>;
   extern template class McMoveManager<3>;
   #endif

}
}
#endif

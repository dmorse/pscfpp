#ifndef RP_MC_MOVE_MANAGER_TPP
#define RP_MC_MOVE_MANAGER_TPP

/*
* PSCF - Polymer Self-Consistent Field
*
* Copyright 2015 - 2025, The Regents of the University of Minnesota
* Distributed under the terms of the GNU General Public License.
*/

#include <util/random/Random.h>
#include <util/global.h>

namespace Pscf {
namespace Rp {

   using namespace Util;

   /*
   * Constructor.
   */
   template <int D, class T>
   McMoveManager<D,T>::McMoveManager(McSimulator<D,T>& simulator,
                                     System<D,T>& system)
    : Base(),
      simulatorPtr_(&simulator),
      systemPtr_(&system),
      randomPtr_(&simulator.random())
   {  ParamComposite::setClassName("McMoveManager"); }

   /*
   * Return a pointer to a new McMoveFactory object.
   */
   template <int D, class T>
   Factory<McMove<D,T> >* McMoveManager<D,T>::newDefaultFactory()
   const
   {  return new typename T::McMoveFactory(*simulatorPtr_); }

   /*
   * Read instructions for creating objects from file.
   */
   template <int D, class T>
   void McMoveManager<D,T>::readParameters(std::istream &in)
   {
      // Read parameters for all McMove objects
      Base::readParameters(in);

      // Allocate and store probabilities
      probabilities_.allocate(Base::size());
      double  totalProbability = 0.0;
      int  iMove;
      for (iMove = 0; iMove < Base::size(); ++iMove) {
         probabilities_[iMove] = (*this)[iMove].probability();
         totalProbability += probabilities_[iMove];
      }

      // Validate the sum of MC move probabilities
      if (std::abs(totalProbability - 1.0) > 1.0E-4) {
         std::string msg = "Error: ";
         msg += "Sum of MC move proababilities differs too much from 1.0";
         std::cout << msg;
         UTIL_THROW(msg.c_str());
      }

      // Store normalized probabilities for use in MC move selection
      for (iMove = 0; iMove < Base::size(); ++iMove) {
         probabilities_[iMove] = probabilities_[iMove]/totalProbability;
         (*this)[iMove].setProbability(probabilities_[iMove]);
      }
   }

   /*
   * Initialize all moves just prior to a run.
   */
   template <int D, class T>
   void McMoveManager<D,T>::setup()
   {
      for (int iMove = 0; iMove < Base::size(); ++iMove) {
         (*this)[iMove].setup();
      }
   }

   /*
   * Choose a McMove at random.
   */
   template <int D, class T>
   McMove<D,T>& McMoveManager<D,T>::chooseMove()
   {
      int iMove;
      iMove = randomPtr_->drawFrom(&probabilities_[0], Base::size());
      return (*this)[iMove];
   }

   /*
   * Output statistics for every move.
   */
   template <int D, class T>
   void McMoveManager<D,T>::output() const
   {
      for (int i = 0; i < Base::size(); ++i) {
         (*this)[i].output();
      }
   }

   /*
   * Log output timing results
   */
   template <int D, class T>
   void McMoveManager<D,T>::outputTimers(std::ostream& out) const
   {
      for (int i = 0; i < Base::size(); ++i) {
         (*this)[i].outputTimers(out);
      }
   }

   /*
   * Clear timers
   */
   template <int D, class T>
   void McMoveManager<D,T>::clearTimers()
   {
      for (int i = 0; i < Base::size(); ++i) {
         (*this)[i].clearTimers();
      }
   }

   /*
   * Decide whether any move needs to store cc fields.
   */
   template <int D, class T>
   bool McMoveManager<D,T>::needsCc()
   {
      for (int i = 0; i < Base::size(); ++i) {
         if ((*this)[i].needsCc()) {
            return true;
         }
      }
      return false;
   }

   /*
   * Decide whether any move needs to store dc fields.
   */
   template <int D, class T>
   bool McMoveManager<D,T>::needsDc()
   {
      for (int i = 0; i < Base::size(); ++i) {
         if ((*this)[i].needsDc()) {
            return true;
         }
      }
      return false;
   }

}
}
#endif

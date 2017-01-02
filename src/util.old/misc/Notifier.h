#ifndef UTIL_NOTIFIER_H
#define UTIL_NOTIFIER_H

/*
* Simpatico - Simulation Package for Polymeric and Molecular Liquids
*
* Copyright 2010 - 2014, The Regents of the University of Minnesota
* Distributed under the terms of the GNU General Public License.
*/

#include <list>

namespace Util
{

   using std::list;

   // Base class for component observers
   template <typename Event> class Observer;
   
   /**
   * Abstract template for a notifier (or subject) in the Observer design 
   * pattern.
   *
   * In the observer design pattern, a Notifier manages a list of registered
   * Observer objects, and provides a method to notify all observers when 
   * some event occurs. A list of observer objects is maintained as a list
   * of Observer pointers. The method Notifier::notifyObservers(Event&) 
   * method calls the the update(Event&) method of every Observer in the 
   * list. 
   * 
   * The typename parameter Event is the type of the object that must be 
   * passed to the update() method of each observer. This type can name
   * either a primitive C data type or a class, but must encode whatever 
   * information is required for any Observer to respond appropriately 
   * when notified. 
   *
   * \ingroup Misc_Module
   */
   template <typename Event>
   class Notifier
   {
   
   public:

      // Compiler default constructor.

      // Compiler destructor.
  
      /**
      * Register an observer.
      *
      * \param observer observer object
      */
      void registerObserver(Observer<Event>& observer);
   
      /**
      * Remove an analyzer observer from the container list.
      *
      * \param observer observer object
      */
      void removeObserver(Observer<Event>& observer);
   
      /**
      * Notify the list of observers about an Event.
      */
      void notifyObservers(const Event& event);
   
   private:
   
      /// A linked list containing component observers.
      std::list<Observer<Event>*> observerPtrs_;

   };
  
   /* 
   * Register an observer (add to list).
   */
   template <typename Event>
   void Notifier<Event>::registerObserver(Observer<Event>& observer)
   {
      observerPtrs_.push_back(&observer);
   }
   
   /* 
   * Remove an observer from the list.
   */
   template <typename Event>
   void Notifier<Event>::removeObserver(Observer<Event>& observer)
   {
      observerPtrs_.remove(&observer);
   }

   /* 
   * Notify observers by calling their update methods.
   */
   template <typename Event>
   void Notifier<Event>::notifyObservers(const Event& event)
   {
      typename std::list< Observer<Event>* >::iterator pos;
      pos = observerPtrs_.begin();
      while (pos != observerPtrs_.end())
      {
         (*pos)->update(event);
         ++pos;
      }
   }


}
#endif 

#ifndef UTIL_SIGNAL_H
#define UTIL_SIGNAL_H

/*
* Simpatico - Simulation Package for Polymeric and Molecular Liquids
*
* Copyright 2010 - 2014, The Regents of the University of Minnesota
* Distributed under the terms of the GNU General Public License.
*/

#include "IFunctor.h"
#include "MethodFunctor.h"

#include <list>

namespace Util
{

   using std::list;

   // -----------------------------------------------------------------
   // Signal with one argument.

   /**
   * Notifier (or subject) in the Observer design pattern.
   *
   * A Signal manages a list of registered functor objects, and provides
   * a void Signal<T>::notify(const T&) method that calls them all with 
   * the same argument.
   *
   * The explicit specialization Signal<void>, or Signal<>, has a notify
   * method void Signal<>::notify() that takes no parameters, which calls
   * a method of each observer that takes no parameters.
   * 
   * \ingroup Util_Signal_Module
   */
   template <typename T=void>
   class Signal
   {
   
   public:

      // Compiler default constructor.

      /**
      * Default constructor.
      */
      Signal(){};
  
      /**
      * Destructor.
      */
      ~Signal();
  
      /**
      * Register an observer.
      *
      * \param observer  observer object (invokes method)
      * \param methodPtr pointer to relevant method
      */
      template <class Observer>
      void addObserver(Observer& observer, void (Observer::*methodPtr)(const T&));

      /**
      * Clear all observerse from list.
      */
      void clear();   

      /**
      * Get number of registered observers.
      */
      int nObserver() const;
   
      /**
      * Notify all observers.
      *
      * This method notifies all registered observers by calling the appropriate
      * method of each observer, passing each the parameter t as argument. The
      * explicit specialization Signal<>, with T=void, is used for notification 
      * methods that take
      *
      * \param t Argument passed to notification methods of all observers.
      */
      void notify(const T& t);

   private:
   
      /// A linked list of functors associated with member functions.
      std::list<IFunctor<T>*> functorPtrs_;

      /// Copy constructor - private to prevent copying.
      Signal(const Signal<T>& other);

      /// Assignment operator - private to prevent assignment.
      Signal<T>& operator = (const Signal<T>& other);

   };

   /*
   * Destructor.
   */
   template <typename T>
   Signal<T>::~Signal()
   {  clear(); }

   /* 
   * Register an observer (add to list).
   */
   template <typename T>
   template <class Observer> void 
   Signal<T>::addObserver(Observer& observer, void (Observer::*methodPtr)(const T&))
   {  functorPtrs_.push_back(new MethodFunctor<Observer, T>(observer, methodPtr)); }

   /* 
   * Notify observers (call associated methods).
   */
   template <typename T>
   void Signal<T>::notify(const T& t)
   {
      typename std::list< IFunctor<T>* >::iterator pos;
      pos = functorPtrs_.begin();
      while (pos != functorPtrs_.end())
      {
         (**pos)(t);
         ++pos;
      }
   }

   /* 
   * Clear all observers.
   *
   * Destroy associated functors.
   */
   template <typename T>
   void Signal<T>::clear()
   {
      typename std::list< IFunctor<T>* >::iterator pos;
      pos = functorPtrs_.begin();
      while (pos != functorPtrs_.end())
      {
         delete *pos;
         ++pos;
      }
      functorPtrs_.clear();
   }

   /* 
   * Get number of registered observers.
   */
   template <typename T>
   int Signal<T>::nObserver() const
   { return functorPtrs_.size(); }


   // -----------------------------------------------------------------
   // Signal with no arguments.

   /**
   * Notifier (or subject) in the Observer design pattern (zero parameters).
   *
   * This explicit specialization of Signal<T> provides a notify method
   * that takes no parameters, and that calls methods of each observer
   * object that take no parameters.
   * 
   * \ingroup Util_Module
   */
   template <>
   class Signal<void>
   {
   
   public:

      /**
      * Default constructor.
      */
      Signal(){};

      /**
      * Destructor.
      */
      ~Signal();
  
      /**
      * Register an observer.
      *
      * \param observer  observer object (invokes method)
      * \param methodPtr pointer to relevant method
      */
      template <class Observer>
      void addObserver(Observer& observer, void (Observer::*methodPtr)());
   
      /**
      * Clear all observerse from list.
      */
      void clear();   

      /**
      * Get number of registered observers.
      */
      int nObserver() const;

      /**
      * Notify all observers.
      */
      void notify();
   
   private:
   
      /// A linked list of functors associated with member functions.
      std::list<IFunctor<>*> functorPtrs_;

      /// Copy constructor - private to prevent copying.
      Signal(const Signal<>& other);

      /// Assignment operator - private to prevent assignment.
      Signal<>& operator = (const Signal<>& other);

   };

   /* 
   * Register an observer (add to list).
   */
   template <class Observer> void 
   Signal<>::addObserver(Observer& observer, void (Observer::*methodPtr)())
   {  functorPtrs_.push_back(new MethodFunctor<Observer>(observer, methodPtr)); }

}
#endif 

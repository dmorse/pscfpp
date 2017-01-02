#ifndef UTIL_OBSERVER_H
#define UTIL_OBSERVER_H 

/*
* Simpatico - Simulation Package for Polymeric and Molecular Liquids
*
* Copyright 2010 - 2014, The Regents of the University of Minnesota
* Distributed under the terms of the GNU General Public License.
*/

namespace Util
{

   /**
   * Abstract class template for observer in the observer design pattern.
   *
   * An Observer is notified of an event by calling its update method.
   * The template class parameter Event is the type of object that is 
   * passed to the update() method as a message about an event.
   *
   * \ingroup Misc_Module
   */
   template <typename Event>
   class Observer
   {

   public:

      // Default constructor.

      /**
      * Destructor.
      */
      virtual ~Observer();
   
      /**
      * Respond to news of an event.
      *
      * \param event Object containing information about the event.
      */
      virtual void update(const Event& event) = 0;

   };

   // Destructor
   template <typename Event>
   Observer<Event>::~Observer()
   {}

}

#endif 

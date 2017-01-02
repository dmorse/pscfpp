#ifndef UTIL_METHOD_FUNCTOR_H
#define UTIL_METHOD_FUNCTOR_H

/*
* Simpatico - Simulation Package for Polymeric and Molecular Liquids
*
* Copyright 2010 - 2014, The Regents of the University of Minnesota
* Distributed under the terms of the GNU General Public License.
*/

#include "IFunctor.h"

namespace Util
{

   /**
   * Functor that wraps a one-argument class member function.
   *
   * The constructor to MethodFunctor<T> takes pointers to an invoking instance
   * of class Object and a member function that takes one T argument. The
   * operator () (const T&) invokes that method on that object.
   *
   * \ingroup Util_Signal_Module
   */
   template <class Object, typename T=void>
   class MethodFunctor : public IFunctor<T>
   {
   public:

      typedef void (Object::*Method1Ptr)(const T&);

      /**
      * Constructor.
      *
      * \param object    invoking object
      * \param methodPtr pointer to member function
      */
      MethodFunctor(Object& object, Method1Ptr methodPtr) 
       : objectPtr_(&object),
         methodPtr_(methodPtr)
      {}

      /**
      * Destructor.
      */
      virtual ~MethodFunctor(){}

      /**
      * Operator ().
      *
      * \param t Parameter passed to method of associated T object.
      */
      virtual void operator () (const T& t)
      {  (objectPtr_->*methodPtr_)(t); }

   private:

      Object*     objectPtr_;
      Method1Ptr  methodPtr_;

   };

   /**
   * Functor that wraps a class member function with no arguments.
   */
   template <class Object>
   class MethodFunctor<Object, void> : public IFunctor<void>
   {
   public:

      typedef void (Object::*Method0Ptr)();

      /**
      * Constructor.
      *
      * \param object    invoking object
      * \param methodPtr pointer to member function
      */
      MethodFunctor(Object& object, Method0Ptr methodPtr) 
       : objectPtr_(&object),
         methodPtr_(methodPtr)
      {}

      /**
      * Destructor.
      */
      virtual ~MethodFunctor(){}

      virtual void operator () ()
      {  (objectPtr_->*methodPtr_)(); }

   private:

      Object*    objectPtr_;
      Method0Ptr methodPtr_;

   };

}
#endif 

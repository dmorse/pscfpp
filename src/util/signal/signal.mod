namespace Util{

   /**
   * \defgroup Util_Signal_Module Signals (Observer Pattern)
   * \ingroup Util_NS_Module
   *
   * Classes used to implement the observer design pattern. A Signal
   * maintains a list of registered "observers" and "notifies" each 
   * observer when Signal::notify() is called, by calling a specific
   * method of each observer object. Observers are stored internally 
   * as a list pointers to IFunctor objects, each of which can be
   * called using an overloaded () operator. Each Functor is created
   * as an instance of the MethodFunctor<T>, which stores a pointer
   * to a T object and to pointer to a method of class T, and which
   * uses the () operator to call a specific method of a specific
   * object.
   *
   * The Signal, IFunctor, and MethodFunctor classes are all templates
   * that take an optional parameter T that represents the typename of
   * of a parameter that should be passed to the notify method of the
   * Signal<T>, which then passes it to the void (const T&) operator 
   * of the IFunctor<T>. In each template, setting typename T to the
   * the default value of T=void invokes a explicit specialization in 
   * which the void Signal<>::notify() and void IFunctor<>::operator () 
   * take no parameters.  An instance of Signal<> is thus a signal that 
   * notifies observers by calling methods that take no arguments, while 
   * a Signal<T> is a signal that notifies observers by calling methods
   * with a signature void (const &T).  MethodFunctor takes two 
   * template parameters: MethodFunctor<ObserverClass, typename T=void> 
   * is a subclass of IFunctor<T> for which the (const T&) operator 
   * calls a specific void (const T&) methodof an observer of type 
   * class ObserverObject. 
   */

}

#ifndef METHOD_FUNCTOR_TEST_H
#define METHOD_FUNCTOR_TEST_H

#include <test/UnitTest.h>
#include <test/UnitTestRunner.h>
#include <util/signal/MethodFunctor.h>
#include <util/global.h>

using namespace Util;

class MethodFunctorTest : public UnitTest 
{

private:

   class Observer0 
   {
    
   public: 
   
      Observer0() : isNotified_(false) {}
   
      void update() { isNotified_ = true; }
   
      void clear() { isNotified_ = false; }
   
      bool isNotified(){ return isNotified_; }
   
   private:
   
      bool isNotified_;
   
   };
   
   class Observer1 
   {
   
   public: 
   
      Observer1() : isNotified_(false), value_(0) {}
   
      void update(const int& value) { isNotified_ = true; value_ = value; }
   
      void clear() { isNotified_ = false; }
   
      bool isNotified(){ return isNotified_; }
   
      int  value(){ return value_; }
   
   private:
   
      bool isNotified_;
      int  value_;
   
   };

public:

   void setUp()
   {};

   void tearDown()
   {};

   void testObserver1() 
   {
      printMethod(TEST_FUNC);
      Observer1        observer;
      MethodFunctor<Observer1, int>  functor(observer, &Observer1::update);
      IFunctor<int>* functorPtr = &functor;

      TEST_ASSERT(!observer.isNotified());
      TEST_ASSERT(observer.value() == 0);

      int value = 3;
      (*functorPtr)(value);

      TEST_ASSERT(observer.isNotified());
      TEST_ASSERT(observer.value() == 3);
   }

   void testObserver0() 
   {
      printMethod(TEST_FUNC);

      Observer0 observer;
      MethodFunctor<Observer0>  functor(observer, &Observer0::update);
      IFunctor<>* functorPtr = &functor;

      TEST_ASSERT(!observer.isNotified());
      (*functorPtr)();
      TEST_ASSERT(observer.isNotified());
   }

};

TEST_BEGIN(MethodFunctorTest)
TEST_ADD(MethodFunctorTest, testObserver1)
TEST_ADD(MethodFunctorTest, testObserver0)
TEST_END(MethodFunctorTest)

#endif

#ifndef SIGNAL_TEST_H
#define SIGNAL_TEST_H

#include <test/UnitTest.h>
#include <test/UnitTestRunner.h>
#include <util/signal/MethodFunctor.h>
#include <util/signal/Signal.h>
#include <util/global.h>

using namespace Util;

class SignalTest : public UnitTest 
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
      Signal<int> signal;
      Observer1 observerA;
      Observer1 observerB;

      TEST_ASSERT(signal.nObserver() == 0);
      signal.addObserver(observerA, &Observer1::update);
      TEST_ASSERT(signal.nObserver() == 1);

      TEST_ASSERT(!observerA.isNotified());
      TEST_ASSERT(observerA.value() == 0);

      TEST_ASSERT(!observerB.isNotified());
      TEST_ASSERT(observerB.value() == 0);

      int value = 3;
      signal.notify(value);

      TEST_ASSERT(observerA.isNotified());
      TEST_ASSERT(observerA.value() == 3);

      TEST_ASSERT(!observerB.isNotified());
      TEST_ASSERT(observerB.value() == 0);

      signal.addObserver(observerB, &Observer1::update);
      TEST_ASSERT(signal.nObserver() == 2);

      value = 4;
      signal.notify(value);

      TEST_ASSERT(observerA.isNotified());
      TEST_ASSERT(observerA.value() == 4);

      TEST_ASSERT(observerB.isNotified());
      TEST_ASSERT(observerB.value() == 4);

      signal.clear();
      TEST_ASSERT(signal.nObserver() == 0);

      TEST_ASSERT(observerA.isNotified());
      TEST_ASSERT(observerA.value() == 4);
      TEST_ASSERT(observerB.isNotified());
      TEST_ASSERT(observerB.value() == 4);

   }

   void testObserver0() 
   {
      printMethod(TEST_FUNC);
      Signal<> signal;
      Observer0 observerA;
      Observer0 observerB;

      TEST_ASSERT(signal.nObserver() == 0);
      signal.addObserver(observerA, &Observer0::update);
      TEST_ASSERT(signal.nObserver() == 1);

      TEST_ASSERT(!observerA.isNotified());
      TEST_ASSERT(!observerB.isNotified());

      signal.notify();

      TEST_ASSERT(observerA.isNotified());
      TEST_ASSERT(!observerB.isNotified());

      signal.addObserver(observerB, &Observer0::update);
      TEST_ASSERT(signal.nObserver() == 2);

      signal.notify();

      TEST_ASSERT(observerA.isNotified());
      TEST_ASSERT(observerB.isNotified());

      signal.clear();
      TEST_ASSERT(signal.nObserver() == 0);

      TEST_ASSERT(observerA.isNotified());
      TEST_ASSERT(observerB.isNotified());

   }

};

TEST_BEGIN(SignalTest)
TEST_ADD(SignalTest, testObserver1)
TEST_ADD(SignalTest, testObserver0)
TEST_END(SignalTest)

#endif

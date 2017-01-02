#ifndef OFFSET_TEST_H
#define OFFSET_TEST_H

#include <test/UnitTest.h>
#include <test/UnitTestRunner.h>

#include <util/misc/Offset.h>
#include <util/global.h>

using namespace Util;

class OffsetTest : public UnitTest 
{

public:

   class A 
   {
   public:
      int  x;
      bool f;

      virtual void action()
      {}

      void outputOffsets(std::ostream& out)
      {
         A a;
         out << "x " << memberOffset<A>(a, &A::x) << std::endl;
         out << "f " << memberOffset<A>(a, &A::f) << std::endl;
      } 

   };

   class B : public A 
   {
   public:
      int y;

      virtual void action() {}

      void outputOffsets(std::ostream& out)
      {
         B b;
         out << "x " << memberOffset<A>(b, &A::x) << std::endl;
         out << "f " << memberOffset<A>(b, &A::f) << std::endl;
         out << "y " << memberOffset<B>(b, &B::y) << std::endl;
         out << "z " << memberOffset<B>(b, &B::w) << std::endl;
      } 

   private:

      int w;

   };

   class C : public B
   {
   public:
      int z;

      virtual void action() {}
   };

   void setUp()
   {};

   void tearDown()
   {};

   void testMemberOffset1() 
   {
      printMethod(TEST_FUNC);

      B b;
      std::cout << std::endl;
      std::cout << memberOffset<B>(b, &B::y) << std::endl;
      b.outputOffsets(std::cout);
   }

   void testMemberOffset2() 
   {
      printMethod(TEST_FUNC);

      B b;
      std::cout << std::endl;
      std::cout << memberOffset<B>(b, &A::x) << std::endl;
   }

   void testBaseOffset() 
   {
      printMethod(TEST_FUNC);

      B b;
      std::cout << std::endl;
      std::cout << baseOffset<B, A>(b) << std::endl;

      C c;
      std::cout << std::endl;
      std::cout << baseOffset<C, B>(c) << std::endl;
   }

};

TEST_BEGIN(OffsetTest)
TEST_ADD(OffsetTest, testMemberOffset1)
TEST_ADD(OffsetTest, testMemberOffset2)
TEST_ADD(OffsetTest, testBaseOffset)
TEST_END(OffsetTest)

#endif

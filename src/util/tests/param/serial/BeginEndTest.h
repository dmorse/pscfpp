#ifndef BEGIN_END_TEST_H
#define BEGIN_END_TEST_H

#include <test/UnitTest.h>
#include <test/UnitTestRunner.h>

#include <util/param/Begin.h>
#include <util/param/End.h>

#include <iostream>
#include <fstream>

using namespace Util;

class BeginEndTest : public UnitTest 
{

public:

   void setUp()
   {
      Label::clear();
   }

   void tearDown()
   {}

   void testBeginConstructor1() {
      printMethod(TEST_FUNC);
      Begin* param;
      param = new Begin("ClassName");
      TEST_ASSERT(param->isRequired());
      TEST_ASSERT(!param->isActive());
      delete param;
   }

   void testBeginConstructor2() {
      printMethod(TEST_FUNC);
      Begin* param;
      bool isRequired = false;
      param = new Begin("ClassName", isRequired);
      TEST_ASSERT(!param->isRequired());
      TEST_ASSERT(!param->isActive());
      delete param;
   }

   void testBeginWrite() {
      printMethod(TEST_FUNC);
      Begin* param;
      param = new Begin("ClassName");
      if (verbose() > 0) {
         param->writeParam(std::cout);
      }
      delete param;
   }

   void testBeginRead1() {
      printMethod(TEST_FUNC);
      Begin *param;
      param = new Begin("ClassName");
      TEST_ASSERT(param->isRequired());
      TEST_ASSERT(!param->isActive());
      std::ifstream in;
      openInputFile("in/Begin", in);
      param->readParam(in);
      TEST_ASSERT(param->isActive());
      if (verbose() > 0) {
         param->writeParam(std::cout);
      }
      delete param;
   }

   void testBeginRead2() {
      printMethod(TEST_FUNC);
      Begin *param;
      bool isRequired = false;
      param = new Begin("ClassName", isRequired);
      TEST_ASSERT(!param->isRequired());
      TEST_ASSERT(!param->isActive());
      std::ifstream in;
      openInputFile("in/Begin", in);
      param->readParam(in);
      TEST_ASSERT(!param->isRequired());
      TEST_ASSERT(param->isActive());
      if (verbose() > 0) {
         param->writeParam(std::cout);
      }
      delete param;
   }

   void testBeginRead3() {
      printMethod(TEST_FUNC);
      Begin *param;
      bool isRequired = false;
      param = new Begin("OtherClass", isRequired);
      TEST_ASSERT(!param->isRequired());
      TEST_ASSERT(!param->isActive());
      std::ifstream in;
      openInputFile("in/Begin", in);
      param->readParam(in);
      TEST_ASSERT(!param->isRequired());
      TEST_ASSERT(!param->isActive());
      delete param;
      param = new Begin("ClassName", isRequired);
      param->readParam(in);
      TEST_ASSERT(!param->isRequired());
      TEST_ASSERT(param->isActive());
      if (verbose() > 0) {
         param->writeParam(std::cout);
      }
      delete param;
   }

   void testEndWrite() {
      printMethod(TEST_FUNC);
      End* param;
      param = new End();
      if (verbose() > 0) {
         param->writeParam(std::cout);
      }
      delete param;
   }

   void testEndRead() {
      printMethod(TEST_FUNC);
      End *param;
      param = new End();
      std::ifstream in;
      openInputFile("in/End", in);
      param->readParam(in);
      if (verbose() > 0) {
         param->writeParam(std::cout);
      }
      delete param;
   }

};

TEST_BEGIN(BeginEndTest)
TEST_ADD(BeginEndTest, testBeginConstructor1)
TEST_ADD(BeginEndTest, testBeginConstructor2)
TEST_ADD(BeginEndTest, testBeginWrite)
TEST_ADD(BeginEndTest, testBeginRead1)
TEST_ADD(BeginEndTest, testBeginRead2)
TEST_ADD(BeginEndTest, testBeginRead3)
TEST_ADD(BeginEndTest, testEndWrite)
TEST_ADD(BeginEndTest, testEndRead)
TEST_END(BeginEndTest)

#endif

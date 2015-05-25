#ifndef LABEL_TEST_H
#define LABEL_TEST_H

#include <test/UnitTest.h>
#include <test/UnitTestRunner.h>

#include <util/param/Label.h>
#include <util/global.h>

#include <iostream>
#include <fstream>

using namespace Util;

class LabelTest : public UnitTest 
{

public:

   void setUp()
   { 
      Label::clear();
   }

   void tearDown()
   {}

   void testLabelConstructor1() 
   {
      printMethod(TEST_FUNC);
      Label label("MyLabel");
   }

   void testLabelConstructor2() 
   {
      printMethod(TEST_FUNC);
      Label label("MyLabel", false);
   }

   void testExtractor1() 
   {
      printMethod(TEST_FUNC);
      Label label("MyLabel");
      std::ifstream in;
      openInputFile("in/Label", in);
      in >> label;
      in >> Label("YourLabel");
      in.close();
   }

   void testExtractor2() 
   {
      printMethod(TEST_FUNC);
      Label label0("WrongLabel", false);
      Label label1("AnotherLabel", false);
      Label label2("MyLabel", false);
      Label label3("YourLabel", true);
      std::ifstream in;
      openInputFile("in/Label", in);
      in >> label0;
      in >> label1;
      in >> label2;
      in >> label3;
      in.close();
   }

   #ifndef UTIL_MPI
   void testExtractor3() 
   {
      printMethod(TEST_FUNC);
      Label label0("WrongLabel", false);
      Label label1("AnotherLabel", false);
      Label label2("MyLabel", false);
      Label label3("YourLabel", true);
      std::ifstream in;
      openInputFile("in/Label", in);
      in >> label0;
      in >> label1;
      try {
         in >> label3;
      } catch (Util::Exception& e) {
         std::cout << "Caught expected Exception" << std::endl;
         Label::clear();
      }
      in.close();
   }
   #endif

   void testInserter() 
   {
      printMethod(TEST_FUNC);
      printEndl();
      Label label("MyLabel");
      std::cout << label;
      std::cout << "\n"; 
   }

};


TEST_BEGIN(LabelTest)
TEST_ADD(LabelTest, testLabelConstructor1)
TEST_ADD(LabelTest, testLabelConstructor2)
TEST_ADD(LabelTest, testExtractor1)
TEST_ADD(LabelTest, testExtractor2)
#ifndef UTIL_MPI
TEST_ADD(LabelTest, testExtractor3)
#endif
TEST_ADD(LabelTest, testInserter)
TEST_END(LabelTest)

#endif

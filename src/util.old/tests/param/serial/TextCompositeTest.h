#ifndef TEXT_COMPOSITE_TEST_H
#define TEXT_COMPOSITE_TEST_H

#include <test/UnitTest.h>
#include <test/UnitTestRunner.h>

#include <util/param/TextComposite.h>

#include <iostream>
#include <fstream>

using namespace Util;

class TextCompositeTest : public UnitTest 
{

public:

   void setUp()
   {}

   void tearDown()
   {}

   void testRead() {
      printMethod(TEST_FUNC);
      TextComposite param;
      std::ifstream in;
      openInputFile("in/TextComposite", in);
      param.readParam(in);
      std::cout << std::endl;
      param.writeParam(std::cout);
   }

};

TEST_BEGIN(TextCompositeTest)
TEST_ADD(TextCompositeTest, testRead)
TEST_END(TextCompositeTest)

#endif

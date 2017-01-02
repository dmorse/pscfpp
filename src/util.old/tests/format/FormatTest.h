#ifndef FORMAT_TEST_H
#define FORMAT_TEST_H

#include <test/UnitTest.h>
#include <test/UnitTestRunner.h>

#include <util/format/Format.h>
#include <util/format/Dbl.h>
#include <util/format/Int.h>
#include <util/format/Str.h>

using namespace Util;

class FormatTest : public UnitTest 
{

public:

   void setUp()
   {};

   void tearDown()
   {};

   void testDbl() 
   { 
      printMethod(TEST_FUNC);
      printEndl();
      std::cout << Dbl(2.0) << Dbl(4.0, 15) 
                << Dbl(3.0, 15, 7) 
                << Dbl(3.5, 15, 7, true) 
                << std::endl; 
   }

   void testInt() 
   { 
      printMethod(TEST_FUNC);
      printEndl();
      std::cout << Int(2) << Int(4, 15) << std::endl; 
   }

   void testStr() 
   { 
      printMethod(TEST_FUNC);
      printEndl();
      std::cout << Str("Hello") << Str("World", 15) << std::endl; 
   }

};

TEST_BEGIN(FormatTest)
TEST_ADD(FormatTest, testDbl)
TEST_ADD(FormatTest, testInt)
TEST_ADD(FormatTest, testStr)
TEST_END(FormatTest)

#endif

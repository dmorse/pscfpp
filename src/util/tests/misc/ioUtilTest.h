#ifndef IO_UTIL_TEST_H
#define IO_UTIL_TEST_H

#include <util/misc/ioUtil.h>
#include <util/global.h>

#ifdef UTIL_MPI
#ifndef TEST_MPI
#define TEST_MPI
#endif
#endif

#include <test/UnitTest.h>
#include <test/UnitTestRunner.h>

using namespace Util;

class ioUtilTest : public UnitTest 
{

public:

   void setUp()
   {};

   void tearDown()
   {};

   void testRStrip() 
   {
      printMethod(TEST_FUNC);

      std::string full("  This has white space   ");
      std::string lean("  This has white space");
      int len = rStrip(full);
      TEST_ASSERT(len == 22);
      TEST_ASSERT(full == lean);
   }

   /*
   * Test of getNextLine for string argument.
   */
   void testGetLine() 
   {
      printMethod(TEST_FUNC);

      std::ifstream in;
      openInputFile("in/getLine", in);

      bool flag;
      std::stringstream line;
      flag = getLine(in, line);
      TEST_ASSERT(flag);
      //std::cout << "|" << line.str() << "|" << std::endl;
      std::string expected("   This is the first line   ");
      TEST_ASSERT(line.str() == expected);

      flag = true;
      while (flag) {
         flag = getLine(in, line);
         TEST_ASSERT(line.str() == "");
      }
      TEST_ASSERT(!flag);
   }

   /*
   * Test of getNextLine for string argument.
   */
   void testGetNextLine1() 
   {
      printMethod(TEST_FUNC);

      std::ifstream in;
      openInputFile("in/getNextLine", in);

      bool flag;
      std::string line;
      flag = getNextLine(in, line);
      TEST_ASSERT(flag);
      std::string expected("   This is the first line");
      // std::cout << "|" << line << "|" << std::endl;
      TEST_ASSERT(line == expected);
      flag = getNextLine(in, line);
      TEST_ASSERT(!flag);
   }

   /*
   * Test of getNextLine for stringstream buffer argument.
   */
   void testGetNextLine2() 
   {
      printMethod(TEST_FUNC);

      std::ifstream in;
      openInputFile("in/getNextLine", in);

      bool flag;
      std::stringstream line;
      flag = getNextLine(in, line);
      TEST_ASSERT(flag);
      std::string expected("   This is the first line");
      // std::cout << "|" << line << "|" << std::endl;
      TEST_ASSERT(line.str() == expected);

      std::string word;
      line >> word;
      TEST_ASSERT(word == "This");
      line >> word;
      TEST_ASSERT(word == "is");
      line >> word;
      TEST_ASSERT(word == "the");
      line >> word;
      TEST_ASSERT(word == "first");
      line >> word;
      TEST_ASSERT(word == "line");
      
      flag = getNextLine(in, line);
      TEST_ASSERT(!flag);
   }

   /*
   * Test getNextLine and checkString.
   */
   void testCheckString() 
   {
      printMethod(TEST_FUNC);

      std::ifstream in;
      openInputFile("in/getNextLine", in);

      bool flag;
      std::stringstream line;
      flag = getNextLine(in, line);
      TEST_ASSERT(flag);
      std::string expected("   This is the first line");
      TEST_ASSERT(line.str() == expected);

      std::string word;
      checkString(line, "This");
      checkString(line, "is");
      checkString(line, "the");
      checkString(line, "first");
      checkString(line, "line");
      
      flag = getNextLine(in, line);
      TEST_ASSERT(!flag);
   }

};

TEST_BEGIN(ioUtilTest)
TEST_ADD(ioUtilTest, testRStrip)
TEST_ADD(ioUtilTest, testGetLine)
TEST_ADD(ioUtilTest, testGetNextLine1)
TEST_ADD(ioUtilTest, testGetNextLine2)
TEST_ADD(ioUtilTest, testCheckString)
TEST_END(ioUtilTest)

#endif

#ifndef MPI_LOGGER_TEST_H
#define MPI_LOGGER_TEST_H

#include <util/mpi/MpiLogger.h>

#define TEST_MPI

#include <mpi.h>
#include <test/UnitTest.h>
#include <test/UnitTestRunner.h>

#include <iostream>

using namespace Util;

class MpiLoggerTest : public UnitTest
{

public:

   MpiLoggerTest()
    : UnitTest()
   {}

   void testLogger() 
   {
      printMethod(TEST_FUNC);
      printEndl();

      MpiLogger logger;
      logger.begin();

      std::cout << "I am processor " << communicator().Get_rank() << std::endl;

      logger.end();
   }

};

TEST_BEGIN(MpiLoggerTest)
TEST_ADD(MpiLoggerTest, testLogger)
TEST_END(MpiLoggerTest)

#endif

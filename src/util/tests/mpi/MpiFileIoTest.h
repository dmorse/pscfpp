#ifndef MPI_FILE_IO_TEST_H
#define MPI_FILE_IO_TEST_H

#include <util/mpi/MpiFileIo.h>

#ifndef TEST_MPI
#define TEST_MPI
#endif

#include <mpi.h>
#include <test/UnitTest.h>
#include <test/UnitTestRunner.h>

#include <iostream>

using namespace Util;

class MpiFileIoTest : public UnitTest
{

public:

   MpiFileIoTest()
    : UnitTest(),
      fileIo_()
   {}

   void testSetCommunicator() 
   {
      printMethod(TEST_FUNC);
      TEST_ASSERT(!fileIo().hasIoCommunicator());
      fileIo().setIoCommunicator(communicator());
      TEST_ASSERT(fileIo().hasIoCommunicator());
      TEST_ASSERT(&fileIo().ioCommunicator() == &communicator());
      fileIo().clearCommunicator();
      TEST_ASSERT(!fileIo().hasIoCommunicator());
   }

   void testIsIoProcessor1() 
   {
      printMethod(TEST_FUNC);
      if (mpiRank() == 0) {
         TEST_ASSERT(fileIo().isIoProcessor());
      } else
      if (mpiRank() == 1) {
         TEST_ASSERT(fileIo().isIoProcessor());
      }
   }

   void testIsIoProcessor2() 
   {
      printMethod(TEST_FUNC);
      fileIo().setIoCommunicator(communicator());
      if (mpiRank() == 0) {
         TEST_ASSERT(fileIo().isIoProcessor());
      } else
      if (mpiRank() == 1) {
         TEST_ASSERT(!fileIo().isIoProcessor());
      }
   }

   MpiFileIo& fileIo()
   { return fileIo_; }

private:

   MpiFileIo       fileIo_;

};

TEST_BEGIN(MpiFileIoTest)
TEST_ADD(MpiFileIoTest, testSetCommunicator)
TEST_ADD(MpiFileIoTest, testIsIoProcessor1)
TEST_ADD(MpiFileIoTest, testIsIoProcessor2)
TEST_END(MpiFileIoTest)

#endif

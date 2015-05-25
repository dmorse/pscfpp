#include "MpiSendRecvTest.h"
#include "MpiFileIoTest.h"
#include "MpiLoaderTest.h"
#include "MpiLoggerTest.h"
#include "MpiTestComposite.h"

using namespace Util;

int main()
{
   MPI::Init();
   Vector::commitMpiType();
   IntVector::commitMpiType();

   
   MpiTestComposite test;
   test.run();

   //TEST_RUNNER(MpiSendRecvTest) test1;
   //test1.run();

   //TEST_RUNNER(MpiFileIoTest) test2;
   //test2.run();

   //TEST_RUNNER(MpiLoaderTest) test3;
   //test3.run();

   //TEST_RUNNER(MpiLoggerTest) test4;
   //test3.run();

   MPI::Finalize();
}

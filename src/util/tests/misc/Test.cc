#include "ExceptionTest.h"
#include "MiscTestComposite.h"
//#include "ioUtilTest.h"
//#include "SetableTest.h"
//#include "BitTest.h"

int main() 
{

   #ifdef UTIL_MPI
   MPI::Init();
   #endif

   #ifndef UTIL_MPI
   TEST_RUNNER(ExceptionTest) test1;
   test1.run();
   #endif

   MiscTestComposite test;
   test.run();

   #if 0
   TEST_RUNNER(ioUtilTest) test2;
   test2.run();

   TEST_RUNNER(SetableTest) test3;
   test3.run();

   TEST_RUNNER(BitTest) test4;
   test4.run();
   #endif

   #ifdef UTIL_MPI
   MPI::Finalize();
   #endif

}

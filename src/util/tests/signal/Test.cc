#include "MethodFunctorTest.h"
#include "SignalTest.h"

int main() 
{

   #ifdef UTIL_MPI
   MPI::Init();
   #endif

   TEST_RUNNER(MethodFunctorTest) test1;
   test1.run();

   TEST_RUNNER(SignalTest) test2;
   test2.run();

   #ifdef UTIL_MPI
   MPI::Finalize();
   #endif

}

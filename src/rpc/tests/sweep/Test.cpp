/*
* This program runs all unit tests in the rpc/tests/sweep directory.
*/ 

#include <util/global.h>
#include "SweepTestComposite.h"

#include <test/TestRunner.h>
#include <test/CompositeTestRunner.h>

int main(int argc, char* argv[])
{
   SweepTestComposite runner;
   runner.run();
}

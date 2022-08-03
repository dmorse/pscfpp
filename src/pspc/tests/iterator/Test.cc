/*
* This program runs all unit tests in the pspc/tests/iterator directory.
*/ 

#include <util/global.h>
#include "IteratorTestComposite.h"

#include <test/TestRunner.h>
#include <test/CompositeTestRunner.h>

int main(int argc, char* argv[])
{
   IteratorTestComposite runner;
   runner.run();
}
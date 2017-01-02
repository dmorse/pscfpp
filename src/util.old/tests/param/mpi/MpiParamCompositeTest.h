#ifndef MPI_PARAM_COMPOSITE_TEST_H
#define MPI_PARAM_COMPOSITE_TEST_H

#include <util/global.h>

#ifdef UTIL_MPI

#include <util/param/ParamComposite.h>
#include <util/param/Factory.h>
#include <util/param/Manager.h>
#include <util/space/Vector.h>
#include <util/space/IntVector.h>
#include <util/misc/FileMaster.cpp>

#ifndef TEST_MPI
#define TEST_MPI
#endif

#include <test/ParamFileTest.h>
#include <test/UnitTestRunner.h>

#include <iostream>
#include <fstream>

using namespace Util;

#include "../ParamTestClasses.h"

class MpiParamCompositeTest : public ParamFileTest
{

   AComposite object;

public:

   void setUp()
   {  
      Label::clear(); 
      ParamComponent::setEcho(false);
   }

   void breakDown()
   {  
      Label::clear(); 
      ParamComponent::setEcho(false);
      // Log::close(); 
   }

   void testReadWrite1() 
   {
      printMethod(TEST_FUNC);
      int     value0;
      int     optInt;
      long    value1;
      double  value2;
      std::string  str;
      int     value3[3];
      double  value4[3];
      double  value5[3][3];
      DArray<double> value6;
      value6.allocate(4);
      Vector    value7;
      IntVector value8;
      DMatrix<double> value9;
      value9.allocate(2,2);
      E         e;
      AManager  manager;

      object.setIoCommunicator(communicator()); 
      //ParamComponent::setEcho(true);

      openFile("in/ParamComposite");
      if (ParamComponent::echo()) Log::file() << std::endl;

      object.readBegin(file(), "AComposite");
      object.read<int>(file(), "value0", value0);
      object.read<long>(file(), "value1", value1);
      object.read<int>(file(), "optInt", optInt, false);
      object.read<double>(file(), "value2", value2);
      object.read<std::string>(file(), "str", str);
      object.readCArray<int>(file(), "value3", value3, 3, false);
      object.readCArray<double>(file(), "value4", value4, 3);
      object.readCArray2D<double>(file(), "value5", value5[0], 2, 2, 3);
      object.readDArray<double>(file(), "value6", value6, 4);
      object.readBlank(file());
      object.read<Vector>(file(), "value7", value7);
      object.read<IntVector>(file(), "value8", value8);
      object.readDMatrix<double>(file(), "value9", value9, 2, 2);
      object.readParamComposite(file(), e);
      object.readParamComposite(file(), manager);
      object.readEnd(file());

      if (mpiRank() == 1) {
         std::cout << std::endl;
         object.writeParam(std::cout);
      }

   }

   void testReadWrite2() 
   {
      printMethod(TEST_FUNC);

      object.setIoCommunicator(communicator()); 
      //ParamComponent::setEcho(true);

      openFile("in/ParamComposite");
      if (ParamComponent::echo()) Log::file() << std::endl;
      object.readParam(file());
      file().close();

      if (mpiRank() == 1) {
         std::cout << std::endl;
         object.writeParam(std::cout);
      }

   }

   void testReadWrite3() 
   {
      printMethod(TEST_FUNC);

      FileMaster    fileMaster;
      std::ofstream logFile;

      fileMaster.setDirectoryId(mpiRank());
      fileMaster.openInputFile("ParamComposite", file());
      object.readParam(file());
      file().close();

      fileMaster.openOutputFile("log", logFile);
      Log::setFile(logFile);
      Log::file() << std::endl;
      object.writeParam(Log::file());
      Log::close();
      
      //logFile << std::endl;
      //object.writeParam(logFile);

   }

   void testReadWrite4() 
   {
      printMethod(TEST_FUNC);

      BComposite optional;

      object.setIoCommunicator(communicator()); 
      optional.setIoCommunicator(communicator()); 
      ParamComponent::setEcho(true);

      openFile("in/ParamComposite");
      if (ParamComponent::echo()) {
         Log::file() << std::endl;
      }
      optional.readParamOptional(file());
      object.readParam(file());
      file().close();

      TEST_ASSERT(object.isRequired());
      TEST_ASSERT(object.isActive());
      TEST_ASSERT(!optional.isRequired());
      TEST_ASSERT(!optional.isActive());

      if (mpiRank() == 1) {
         std::cout << std::endl;
         optional.writeParam(std::cout);
         object.writeParam(std::cout);
      }

   }

};


TEST_BEGIN(MpiParamCompositeTest)
TEST_ADD(MpiParamCompositeTest, testReadWrite1)
TEST_ADD(MpiParamCompositeTest, testReadWrite2)
TEST_ADD(MpiParamCompositeTest, testReadWrite3)
TEST_ADD(MpiParamCompositeTest, testReadWrite4)
TEST_END(MpiParamCompositeTest)

#endif // ifdef UTIL_MPI
#endif // ifndef MPI_PARAM_COMPOSITE_TEST_H

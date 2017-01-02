#ifndef MPI_MANAGER_TEST_H
#define MPI_MANAGER_TEST_H

#ifdef UTIL_MPI

#include <util/global.h>
#include <util/param/ParamComponent.h>
#include <util/param/ParamComposite.h>
#include <util/param/Manager.h>
#include <util/param/Factory.h>
#include <util/misc/FileMaster.cpp>

#ifndef TEST_MPI
#define TEST_MPI
#endif

#include <test/UnitTest.h>
#include <test/UnitTestRunner.h>

using namespace Util;

#include "../ParamTestClasses.h"

class MpiManagerTest : public UnitTest 
{

public:

   void setUp() {
      Label::clear();
   }

   void breakDown() 
   {
      Label::clear();
      ParamComponent::setEcho(false);
   }

   /**
   * Test factory using string.
   */
   void testFactory() 
   {
      printMethod(TEST_FUNC);

      AFactory    factory;
      A*          ptr;
      std::string name("B");

      ptr = factory.factory(name);

      TEST_ASSERT(ptr->className() == "B");
 
      // if (mpiRank() == 1) {
      //   std::cout <<  ptr->className() << std::endl;
      //   std::cout.flush();
      //}
   }

   /**
   * Read single file in/Factory
   */
   void testFactoryReadObject() 
   {
      printMethod(TEST_FUNC);

      AFactory    factory;
      AManager    manager;
      std::string className;
      A*          ptr = 0;
      bool        isEnd = false;

      manager.setIoCommunicator(communicator());

      std::ifstream in;
      if (manager.isIoProcessor()) {
         in.open("in/Factory");
      }
      ptr = factory.readObject(in, manager, className, isEnd);
      TEST_ASSERT(ptr->className() == "B");

      //if (mpiRank() == 1) {
      //   //std::cout <<  ptr->className() << std::endl;
      //   ptr->writeParam(std::cout);
      //}
 
   }

   /**
   * Read single file in/Manager.
   */
   void testManager1() 
   {
      printMethod(TEST_FUNC);
      AManager       manager;
      CustomAFactory factory;
      manager.setIoCommunicator(communicator());
      manager.addSubfactory(factory);

      std::ifstream in;
      if (mpiRank() == 0) {
         in.open("in/Manager");
         if (ParamComponent::echo()) printEndl();
      }
      manager.readParam(in);
      TEST_ASSERT(manager.isRequired());
      TEST_ASSERT(manager.isActive());
      if (mpiRank() == 1) {
         std::cout << std::endl;
         manager.writeParam(std::cout);
         std::cout.flush();
      }

   }

   /**
   * Read optional manager from single file in/Manager.
   */
   void testManager2() 
   {
      printMethod(TEST_FUNC);
      AManager       manager;
      CustomAFactory factory;
      manager.setIoCommunicator(communicator());
      manager.addSubfactory(factory);
      TEST_ASSERT(manager.isRequired());
      TEST_ASSERT(manager.isActive());

      std::ifstream in;
      if (mpiRank() == 0) {
         in.open("in/Manager");
         if (ParamComponent::echo()) printEndl();
      }
      manager.readParamOptional(in);
      TEST_ASSERT(!manager.isRequired());
      TEST_ASSERT(manager.isActive());
      if (mpiRank() == 1) {
         std::cout << std::endl;
         manager.writeParam(std::cout);
         std::cout.flush();
      }

   }

   /**
   * Read multiple files n/Manager, where n=mpiRank().
   */
   void testManager3() 
   {
      printMethod(TEST_FUNC);

      AManager      manager;
      FileMaster    fileMaster;
      std::ifstream in;
      std::ofstream log;

      fileMaster.setDirectoryId(mpiRank());
      fileMaster.openInputFile("Manager", in);
      fileMaster.openOutputFile("Manager.log", log);
      if (ParamComponent::echo()) printEndl();

      manager.readParam(in);

      Log::setFile(log);
      Log::file() << std::endl;
      manager.writeParam(Log::file());
      Log::file().flush();
      Log::close();
      
      //log << std::endl;
      //manager.writeParam(log);
      //log.flush();

   }

};

TEST_BEGIN(MpiManagerTest)
TEST_ADD(MpiManagerTest, testFactory)
TEST_ADD(MpiManagerTest, testFactoryReadObject)
TEST_ADD(MpiManagerTest, testManager1)
TEST_ADD(MpiManagerTest, testManager2)
TEST_ADD(MpiManagerTest, testManager3)
TEST_END(MpiManagerTest)

#endif
#endif

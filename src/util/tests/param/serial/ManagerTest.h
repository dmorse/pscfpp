#ifndef MANAGER_TEST_H
#define MANAGER_TEST_H

#include <util/param/ParamComposite.h>
#include <util/param/Manager.h>
#include <util/param/Factory.h>

#include <test/UnitTest.h>
#include <test/UnitTestRunner.h>

using namespace Util;

#include "../ParamTestClasses.h"

class ManagerTest : public UnitTest 
{

public:

   void testFactory() 
   {
      printMethod(TEST_FUNC);

      AFactory    factory;
      A*          ptr;
      std::string name("B");

      ptr = factory.factory(name);
      TEST_ASSERT(ptr->className() == name);
 
      printEndl();
      std::cout <<  "Classname = " << ptr->className() << std::endl;

   }

   void testCustomFactory() 
   {
      printMethod(TEST_FUNC);

      CustomAFactory factory;
      A*             ptr;
      std::string name("D");

      ptr = factory.factory(name);
      TEST_ASSERT(ptr != 0);
      TEST_ASSERT(ptr->className() == name);
 
      printEndl();
      std::cout <<  "Classname = " << ptr->className() << std::endl;

   }

   void testManager() 
   {
      printMethod(TEST_FUNC);
      AManager manager;
      std::ifstream in;
      openInputFile("in/Manager", in);
      manager.readParam(in);

      printEndl();
      manager.writeParam(std::cout);

   }

   void testCustomManager() 
   {
      printMethod(TEST_FUNC);
      AManager manager;
      std::ifstream in;
      openInputFile("in/CustomManager", in);

      CustomAFactory factory;
      manager.addSubfactory(factory);
      manager.readParam(in);

      printEndl();
      manager.writeParam(std::cout);

   }

};

TEST_BEGIN(ManagerTest)
TEST_ADD(ManagerTest, testFactory)
TEST_ADD(ManagerTest, testCustomFactory)
TEST_ADD(ManagerTest, testManager)
TEST_ADD(ManagerTest, testCustomManager)
TEST_END(ManagerTest)

#endif

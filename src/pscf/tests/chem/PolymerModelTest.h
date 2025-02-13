#ifndef POLYMER_MODEL_TEST_H
#define POLYMER_MODEL_TEST_H

#include <test/UnitTest.h>
#include <test/UnitTestRunner.h>

#include <pscf/chem/PolymerModel.h>
#include <util/param/ScalarParam.h>

#include <fstream>

using namespace Pscf;
using namespace Util;

class PolymerModelTest : public UnitTest 
{

public:

   void setUp()
   {
      //setVerbose(1);
   }

   void tearDown()
   {
      setVerbose(0);
      PolymerModel::setModel(PolymerModel::Thread);
   }

   void testReadWrite() {
      printMethod(TEST_FUNC);

      std::ifstream in;
      openInputFile("in/PolymerModel", in);

      PolymerModel::Type b;
      in >> b;
      TEST_ASSERT(b == PolymerModel::Thread);
      if (verbose() > 0) {
         printEndl();
         std::cout << b << "  ";
      }
      in >> b;
      TEST_ASSERT(b == PolymerModel::Bead);
      if (verbose() > 0) {
         std::cout << b << std::endl ;
      }

      // If uncommented out, this one fails to read "Thingy"
      //in >> b;
   }

   void testGlobal() {
      printMethod(TEST_FUNC);

      TEST_ASSERT(PolymerModel::model() == PolymerModel::Thread);
      TEST_ASSERT(PolymerModel::isThread());
      int nSetInit = PolymerModel::nSet();

      PolymerModel::setModel(PolymerModel::Thread);
      TEST_ASSERT(PolymerModel::model() == PolymerModel::Thread);
      TEST_ASSERT(PolymerModel::isThread());
      TEST_ASSERT(!PolymerModel::isBead());
      TEST_ASSERT(nSetInit + 1 == PolymerModel::nSet());

      PolymerModel::setModel(PolymerModel::Bead);
      TEST_ASSERT(PolymerModel::model() == PolymerModel::Bead);
      TEST_ASSERT(!PolymerModel::isThread());
      TEST_ASSERT(PolymerModel::isBead());
      TEST_ASSERT(nSetInit + 2 == PolymerModel::nSet());
   }

   void testParam() {
      printMethod(TEST_FUNC);

      PolymerModel::Type model;
      bool isRequired = true;
      Util::ScalarParam<PolymerModel::Type> param("model", model, isRequired);
   }

};

TEST_BEGIN(PolymerModelTest)
TEST_ADD(PolymerModelTest, testReadWrite)
TEST_ADD(PolymerModelTest, testGlobal)
TEST_ADD(PolymerModelTest, testParam)
TEST_END(PolymerModelTest)

#endif

#ifndef PARAM_COMPOSITE_TEST_H
#define PARAM_COMPOSITE_TEST_H

#include <util/global.h>

#include <util/param/ParamComposite.h>
#include <util/param/Factory.h>
#include <util/param/Manager.h>
#include <util/space/Vector.h>
#include <util/space/IntVector.h>

#include <util/archives/MemoryCounter.h>
#include <util/archives/MemoryOArchive.h>
#include <util/archives/MemoryIArchive.h>

#include <util/archives/Serializable.h>
#include <util/archives/Serializable_includes.h>

#include <test/UnitTest.h>
#include <test/UnitTestRunner.h>

#include <iostream>
#include <fstream>
#include <string>

using namespace Util;

#include "../ParamTestClasses.h"

class ParamCompositeTest : public UnitTest
{

  ParamComposite paramComposite_;
  std::ifstream  file_;

public:

   void setUp()
   {
      Label::clear();
      ParamComponent::setEcho(false);
   }

   void testConstructor() 
   {}

   #if 0
   void testAddWrite() 
   {
      printMethod(TEST_FUNC);

      int     value0 = 4;
      long    value1 = 25;
      double  value2 = 3.0;
      std::string str("Thunk");
      int     value3[3];
      value3[0] = 4;
      value3[1] = 5;
      value3[2] = 6;
      double value4[3];
      value4[0] = 4.0;
      value4[1] = 5.0;
      value4[2] = 6.0;
      double value5[3][3]; // Physically 3 x 3, logicall 2 x 2
      value5[0][0] = 4.0;
      value5[0][1] = 5.0;
      value5[1][0] = 6.0;
      value5[1][1] = 7.0;
      DArray<double> value6;
      value6.allocate(4);
      value6[0] = 14.5;
      value6[1] = 15.4;
      value6[2] = 16.3;
      value6[3] = 16.2;
      FArray<double, 4> value7;
      value7[0] = 13.5;
      value7[1] = 14.4;
      value7[2] = 15.3;
      value7[3] = 15.2;

      paramComposite_.addBegin("AComposite");
      paramComposite_.add<int>("value0", value0);
      paramComposite_.add<long>("value1", value1);
      paramComposite_.add<double>("value2", value2);
      paramComposite_.add<std::string>("str", str);
      paramComposite_.addCArray<int>("value3", value3, 3);
      paramComposite_.addCArray<double>("value4", value4, 3);
      paramComposite_.addCArray2D<double>("value5", &value5[0][0], 2, 2, 3);
      paramComposite_.addDArray<double>("value6", value6, 4);
      paramComposite_.addFArray<double, 4>("value7", value7);
      paramComposite_.addEnd();

      printEndl();
      paramComposite_.writeParam(std::cout);
   }
   #endif

   void testReadWrite1() 
   {
      printMethod(TEST_FUNC);
      int     value0;
      int     optInt;
      long    value1;
      double  value2;
      std::string str;
      int     value3[3];
      double  value4[3];
      double  value5[3][3]; // Physically 3 x 3, logically 2 x 2
      DArray<double> value6;
      value6.allocate(4);
      Vector    value7;
      IntVector value8;
      DMatrix<double>  value9;
      value9.allocate(2, 2);
      E e;
      AManager  manager;

      openInputFile("in/ParamComposite", file_);

      paramComposite_.setEcho();
      if (ParamComponent::echo()) std::cout << std::endl;
      paramComposite_.readBegin(file_, "AComposite");
      paramComposite_.read<int>(file_, "value0", value0);
      paramComposite_.readOptional<int>(file_, "optInt", optInt); // Optional parameter
      paramComposite_.read<long>(file_, "value1", value1);
      paramComposite_.read<double>(file_, "value2", value2);
      paramComposite_.read<std::string>(file_, "str", str);
      paramComposite_.readBlank(file_);
      paramComposite_.readCArray<int>(file_, "value3", value3, 3);
      paramComposite_.readCArray<double>(file_, "value4", value4, 3);
      paramComposite_.readCArray2D<double>(file_, "value5", &value5[0][0], 2, 2, 3);
      paramComposite_.readDArray<double>(file_, "value6", value6, 4);
      paramComposite_.read<Vector>(file_, "value7", value7);
      paramComposite_.read<IntVector>(file_, "value8", value8);
      paramComposite_.readDMatrix<double>(file_, "value9", value9, 2, 2);
      paramComposite_.readParamComposite(file_, e);
      paramComposite_.readParamComposite(file_, manager);
      paramComposite_.readEnd(file_);
      paramComposite_.setEcho(false);

      printEndl();
      paramComposite_.writeParam(std::cout);
   }

   void testReadWrite2() 
   {
      printMethod(TEST_FUNC);

      ParamComponent::setEcho(true);
      if (ParamComponent::echo()) std::cout << std::endl;

      BComposite bcomp;
      AComposite acomp;

      openInputFile("in/ParamComposite", file_);
      bcomp.readParamOptional(file_);
      acomp.readParam(file_);

      printEndl();
      bcomp.writeParam(std::cout);
      acomp.writeParam(std::cout);
   }

   void testSaveLoadWrite() 
   {
      printMethod(TEST_FUNC);

      int     value0 = 4;
      long    value1 = 25;
      double  value2 = 3.0;
      int     value3[3];
      value3[0] = 4;
      value3[1] = 5;
      value3[2] = 6;
      double value4[3];
      value4[0] = 4.0;
      value4[1] = 5.0;
      value4[2] = 6.0;
      double value5[3][3]; // Physically 3 x 3, logically 2 x 2
      value5[0][0] = 4.0;
      value5[0][1] = 5.0;
      value5[1][0] = 6.0;
      value5[1][1] = 7.0;
      DArray<double> value6;
      value6.allocate(4);
      value6[0] = 14.5;
      value6[1] = 15.4;
      value6[2] = 16.3;
      value6[3] = 16.2;
      FArray<double, 4> value7;
      value7[0] = 13.5;
      value7[1] = 14.4;
      value7[2] = 15.3;
      value7[3] = 15.2;

      printEndl();

      Serializable::OArchive oar;
      openOutputFile("out/save.bin", oar.file());
      oar << value0;
      oar << value1;
      oar << value2;
      oar.pack(value3, 3);
      oar.pack(value4, 3);
      oar.pack(value5[0], 2, 2, 3);
      oar << value6;
      oar << value7;
      oar.file().close();

      Serializable::IArchive iar;
      openInputFile("out/save.bin", iar.file());

      int     cValue0;
      long    cValue1;
      double  cValue2;
      int     cValue3[3];
      double  cValue4[3];
      double  cValue5[3][3];
      DArray<double> cValue6;
      cValue6.allocate(4);
      FArray<double, 4> cValue7;

      ParamComposite clone;
      clone.addBegin("AComposite");
      clone.loadParameter<int>(iar, "value0", cValue0);
      clone.loadParameter<long>(iar, "value1", cValue1);
      clone.loadParameter<double>(iar, "value2", cValue2);
      clone.loadCArray<int>(iar, "value3", cValue3, 3);
      clone.loadCArray<double>(iar, "value4", cValue4, 3);
      clone.loadCArray2D<double>(iar, "value5", &cValue5[0][0], 2, 2, 3);
      clone.loadDArray<double>(iar, "value6", cValue6, 4);
      clone.loadFArray<double, 4>(iar, "value7", cValue7);
      clone.addEnd();

      TEST_ASSERT(cValue0 == 4);
      TEST_ASSERT(cValue1 == 25);
      TEST_ASSERT(eq(cValue2, 3.0));
      TEST_ASSERT(eq(cValue3[0], 4));
      TEST_ASSERT(eq(cValue3[1], 5));
      TEST_ASSERT(eq(cValue3[2], 6));
      TEST_ASSERT(eq(cValue4[0], 4.0));
      TEST_ASSERT(eq(cValue4[1], 5.0));
      TEST_ASSERT(eq(cValue4[2], 6.0));
      TEST_ASSERT(eq(cValue5[0][0], 4.0));
      TEST_ASSERT(eq(cValue5[0][1], 5.0));
      TEST_ASSERT(eq(cValue5[1][0], 6.0));
      TEST_ASSERT(eq(cValue5[1][1], 7.0));
      TEST_ASSERT(eq(cValue6[0], 14.5));
      TEST_ASSERT(eq(cValue6[1], 15.4));
      TEST_ASSERT(eq(cValue6[2], 16.3));
      TEST_ASSERT(eq(cValue6[3], 16.2));
      TEST_ASSERT(eq(cValue7[0], 13.5));
      TEST_ASSERT(eq(cValue7[1], 14.4));
      TEST_ASSERT(eq(cValue7[2], 15.3));
      TEST_ASSERT(eq(cValue7[3], 15.2));
      
      printEndl();
      clone.writeParam(std::cout);
   }

   void testReadSaveLoadWrite1() 
   {
      printMethod(TEST_FUNC);
      int     value0;
      long    value1;
      double  value2;
      std::string str;
      int     value3[3];
      double  value4[3];
      double  value5[3][3]; // Physically 3 x 3, logicall 2 x 2
      DArray<double> value6;
      value6.allocate(4);
      Vector    value7;
      IntVector value8;
      DMatrix<double>  value9;
      value9.allocate(2, 2);
      E e;
      AManager  manager;

      openInputFile("in/ParamComposite", file_);

      //paramComposite_.setEcho();
      paramComposite_.readBegin(file_, "AComposite");
      paramComposite_.read<int>(file_, "value0", value0);
      paramComposite_.read<long>(file_, "value1", value1);
      paramComposite_.read<double>(file_, "value2", value2);
      paramComposite_.read<std::string>(file_, "str", str);
      paramComposite_.readBlank(file_);
      paramComposite_.readCArray<int>(file_, "value3", value3, 3);
      paramComposite_.readCArray<double>(file_, "value4", value4, 3);
      paramComposite_.readCArray2D<double>(file_, "value5", &value5[0][0], 2, 2, 3);
      paramComposite_.readDArray<double>(file_, "value6", value6, 4);
      paramComposite_.read<Vector>(file_, "value7", value7);
      paramComposite_.read<IntVector>(file_, "value8", value8);
      paramComposite_.readDMatrix<double>(file_, "value9", value9, 2, 2);
      paramComposite_.readParamComposite(file_, e);
      paramComposite_.readParamComposite(file_, manager);
      paramComposite_.readEnd(file_);
      paramComposite_.setEcho(false);

      Serializable::OArchive oar;
      openOutputFile("out/save2.bin", oar.file());
      //paramComposite_.save(oar);
      oar << value0;
      oar << value1;
      oar << value2;
      oar << str;
      oar.pack(value3, 3);
      oar.pack(value4, 3);
      oar.pack(value5[0], 2, 2, 3);
      oar << value6;
      oar << value7;
      oar << value8;
      oar << value9;
      oar.file().close();

      int     cValue0;
      long    cValue1;
      double  cValue2;
      std::string cStr;
      int     cValue3[3];
      double  cValue4[3];
      double  cValue5[2][2];
      DArray<double> cValue6;
      cValue6.allocate(4);
      Vector    cValue7;
      IntVector cValue8;
      DMatrix<double>  cValue9;
      cValue9.allocate(2, 2);
      E cE;
      AManager cManager;

      Serializable::IArchive iar;
      openInputFile("out/save2.bin", iar.file());

      ParamComposite clone;
      clone.addBegin("AComposite");
      clone.loadParameter<int>(iar, "value0", cValue0);
      clone.loadParameter<long>(iar, "value1", cValue1);
      clone.loadParameter<double>(iar, "value2", cValue2);
      clone.loadParameter<std::string>(iar, "str", cStr);
      clone.addBlank();
      clone.loadCArray<int>(iar, "value3", cValue3, 3);
      clone.loadCArray<double>(iar, "value4", cValue4, 3);
      clone.loadCArray2D<double>(iar, "value5", &cValue5[0][0], 2, 2, 3);
      clone.loadDArray<double>(iar, "value6", cValue6, 4);
      clone.loadParameter<Vector>(iar, "value7", cValue7);
      clone.loadParameter<IntVector>(iar, "value8", cValue8);
      clone.loadDMatrix<double>(iar, "value9", cValue9, 2, 2);
      //clone.loadParamComposite(iar, cE);
      //clone.loadParamComposite(iar, cManager);
      clone.addEnd();

      printEndl();
      clone.writeParam(std::cout);
   }

   void testReadSaveLoadWrite2() 
   {
      printMethod(TEST_FUNC);

      openInputFile("in/ParamComposite", file_);

      ParamComponent::setEcho(true);
      printEndl();

      BComposite absent;
      AComposite original;
      absent.readParamOptional(file_);
      original.readParam(file_);

      // printEndl();
      // original.writeParam(std::cout);

      Serializable::OArchive oar;
      openOutputFile("out/save1.bin", oar.file());
      absent.saveOptional(oar);
      original.save(oar);
      oar.file().close();

      BComposite absentClone;
      AComposite clone;
      Serializable::IArchive iar;
      openInputFile("out/save1.bin", iar.file());
      absentClone.loadOptional(iar);
      clone.load(iar);

      printEndl();
      absentClone.writeParam(std::cout);
      clone.writeParam(std::cout);
   }

   void testMemoryArchiveSerialize() 
   {
      printMethod(TEST_FUNC);

      openInputFile("in/ParamComposite", file_);

      AComposite original;
      original.readParam(file_);

      MemoryCounter  car;
      car << original;
      int size = car.size();

      MemoryOArchive oar;
      oar.allocate(size);
      oar << original;

      MemoryIArchive iar;
      iar = oar;

      AComposite clone;
      iar >> clone;

      printEndl();
      clone.writeParam(std::cout);
   }


};

TEST_BEGIN(ParamCompositeTest)
TEST_ADD(ParamCompositeTest, testConstructor)
//TEST_ADD(ParamCompositeTest, testAddWrite)
TEST_ADD(ParamCompositeTest, testReadWrite1)
TEST_ADD(ParamCompositeTest, testReadWrite2)
TEST_ADD(ParamCompositeTest, testSaveLoadWrite)
TEST_ADD(ParamCompositeTest, testReadSaveLoadWrite1)
TEST_ADD(ParamCompositeTest, testReadSaveLoadWrite2)
//TEST_ADD(ParamCompositeTest, testMemoryArchiveSerialize)
TEST_END(ParamCompositeTest)

#endif // ifndef PARAM_COMPOSITE_TEST_H

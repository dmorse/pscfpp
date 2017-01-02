#ifndef PARAMETER_TEST_H
#define PARAMETER_TEST_H

#include <test/UnitTest.h>
#include <test/UnitTestRunner.h>

#include <util/param/Parameter.h>
#include <util/param/ScalarParam.h>
#include <util/param/CArrayParam.h>
#include <util/param/DArrayParam.h>
#include <util/param/FArrayParam.h>
#include <util/param/CArray2DParam.h>
#include <util/param/DMatrixParam.h>
#include <util/param/DSymmMatrixParam.h>

using namespace Util;

class ParameterTest : public UnitTest 
{

public:

   void setUp()
   { 
      setVerbose(2); 
      Label::clear();
      // ParamComponent::setEcho(true);
   }

   void tearDown()
   {
      Label::clear();
      ParamComponent::setEcho(false);
   }

   void testParamIntConstructor1() {
      printMethod(TEST_FUNC);
      int        requiredVal = 4;
      Parameter* requiredPrm;
      requiredPrm = new ScalarParam<int>("Required", requiredVal);
      delete requiredPrm;
   }

   void testParamIntConstructor2() {
      printMethod(TEST_FUNC);
      int        requiredVal = 4;
      Parameter* requiredPrm;
      requiredPrm = new ScalarParam<int>("Required", requiredVal, false);
      delete requiredPrm;
   }

   void testParamIntWrite() {
      printMethod(TEST_FUNC);
      int        requiredVal = 4;
      Parameter *requiredPrm;
      requiredPrm = new ScalarParam<int>("Required", requiredVal);

      if (verbose() > 0) {
         printEndl();
         requiredPrm->writeParam(std::cout);
      }
      delete requiredPrm;
   }

   void testParamIntRead1() 
   {
      printMethod(TEST_FUNC);
      int        requiredVal;
      Parameter *requiredPrm;
      requiredPrm  = new ScalarParam<int>("Required", requiredVal);
      std::ifstream in;
      openInputFile("in/ScalarParamInt", in);
      if (ParamComponent::echo()) std::cout << std::endl;
      requiredPrm->readParam(in);
      in.close();
      TEST_ASSERT(Label::isClear());
      TEST_ASSERT(requiredPrm->label() == "Required");
      TEST_ASSERT(requiredVal == 36);
      if (verbose() > 0) {
         printEndl();
         requiredPrm->writeParam(std::cout);
      }
      delete requiredPrm;
   }

   void testParamIntRead2() 
   {
      /// Test absentPrm optional parameter
      printMethod(TEST_FUNC);
      int absentVal;
      int requiredVal;
      Parameter* absentPrm = new ScalarParam<int>("Absent", absentVal, false);
      Parameter* requiredPrm  = new ScalarParam<int>("Required", requiredVal);
      std::ifstream in;
      openInputFile("in/ScalarParamInt", in);
      if (ParamComponent::echo()) std::cout << std::endl;
      absentPrm->readParam(in);
      TEST_ASSERT(!Label::isClear());
      requiredPrm->readParam(in);
      TEST_ASSERT(Label::isClear());
      TEST_ASSERT(requiredPrm->label() == "Required");
      //TEST_ASSERT(requiredVal == 36);
      if (verbose() > 0) {
         printEndl();
         requiredPrm->writeParam(std::cout);
      }
      in.close();
      delete absentPrm;
      delete requiredPrm;
   }

   void testParamIntReadSaveLoad()
   {
      printMethod(TEST_FUNC);
      int absentVal;
      int requiredVal;
      int presentVal;
      Parameter *absentPrm;
      Parameter *requiredPrm;
      Parameter *presentPrm;
      absentPrm = new ScalarParam<int>("Absent", absentVal, false);
      requiredPrm = new ScalarParam<int>("Required", requiredVal);
      presentPrm = new ScalarParam<int>("Present", presentVal, false);

      //Parameter::setEcho(true);
      if (ParamComponent::echo()) printEndl();

      std::ifstream in;
      openInputFile("in/ScalarParamInt", in);
      absentPrm->readParam(in);
      TEST_ASSERT(!Label::isClear());
      requiredPrm->readParam(in);
      TEST_ASSERT(Label::isClear());
      presentPrm->readParam(in);
      TEST_ASSERT(Label::isClear());

      TEST_ASSERT(!absentPrm->isRequired());
      TEST_ASSERT(!absentPrm->isActive());
      TEST_ASSERT(requiredPrm->isRequired());
      TEST_ASSERT(requiredPrm->isActive());
      TEST_ASSERT(!presentPrm->isRequired());
      TEST_ASSERT(presentPrm->isActive());

      #if 0
      if (verbose() > 0) {
         printEndl();
         absentPrm->writeParam(std::cout);
         requiredPrm->writeParam(std::cout);
         presentPrm->writeParam(std::cout);
      }
      #endif

      // Save to archive
      Serializable::OArchive oar;
      openOutputFile("out/binary", oar.file());
      Parameter::saveOptional(oar, absentVal, false);
      oar << requiredVal;
      Parameter::saveOptional(oar, presentVal, true);
      oar.file().close();

      // Load from archive
      int absentVal2;
      int requiredVal2;
      int presentVal2;
      Parameter* absentPrm2;
      Parameter* requiredPrm2;
      Parameter* presentPrm2;
      absentPrm2 = 
           new ScalarParam<int>("Absent", absentVal2, false);
      requiredPrm2 = 
           new ScalarParam<int>("Required", requiredVal2);
      presentPrm2 = 
           new ScalarParam<int>("Present", presentVal2, false);

      Serializable::IArchive iar;
      openInputFile("out/binary", iar.file());
      absentPrm2->load(iar);
      requiredPrm2->load(iar);
      presentPrm2->load(iar);
      iar.file().close();
      TEST_ASSERT(!absentPrm2->isRequired());
      TEST_ASSERT(!absentPrm2->isActive());
      TEST_ASSERT(requiredPrm2->isRequired());
      TEST_ASSERT(!presentPrm2->isRequired());
      TEST_ASSERT(presentPrm2->isActive());
      TEST_ASSERT(requiredVal2 == requiredVal);
      TEST_ASSERT(presentVal2 == presentVal);

      if (verbose() > 0) {
         printEndl();
         absentPrm2->writeParam(std::cout);
         requiredPrm2->writeParam(std::cout);
         presentPrm2->writeParam(std::cout);
      }

      delete absentPrm;
      delete requiredPrm;
      delete presentPrm;
      delete absentPrm2;
      delete requiredPrm2;
      delete presentPrm2;
   }

   void testParamDoubleWrite() {
      printMethod(TEST_FUNC);
      double requiredVal = 4.0;
      Parameter     *requiredPrm;
      requiredPrm = new ScalarParam<double>("Required", requiredVal);
      if (verbose() > 0) {
         printEndl();
         requiredPrm->writeParam(std::cout);
      }
      delete requiredPrm;
   }

   void testParamDoubleRead() {
      printMethod(TEST_FUNC);
      double requiredVal;
      Parameter *requiredPrm;
      requiredPrm = new ScalarParam<double>("Required", requiredVal);
      std::ifstream in;
      openInputFile("in/ScalarParamDouble", in);
      if (ParamComponent::echo()) std::cout << std::endl;
      requiredPrm->readParam(in);
      in.close();
      if (verbose() > 0) {
         printEndl();
         requiredPrm->writeParam(std::cout);
      }
      delete requiredPrm;
   }

   void testParamDoubleReadSaveLoad()
   {
      printMethod(TEST_FUNC);
      double absentVal;
      double requiredVal;
      double presentVal;
      Parameter *absentPrm;
      Parameter *requiredPrm;
      Parameter *presentPrm;
      absentPrm = new ScalarParam<double>("Absent", absentVal, false);
      requiredPrm = new ScalarParam<double>("Required", requiredVal);
      presentPrm = new ScalarParam<double>("Present", presentVal, false);

      //Parameter::setEcho(true);
      if (ParamComponent::echo()) printEndl();

      std::ifstream in;
      openInputFile("in/ScalarParamDouble", in);
      absentPrm->readParam(in);
      TEST_ASSERT(!Label::isClear());
      requiredPrm->readParam(in);
      TEST_ASSERT(Label::isClear());
      presentPrm->readParam(in);
      TEST_ASSERT(Label::isClear());

      TEST_ASSERT(!absentPrm->isRequired());
      TEST_ASSERT(!absentPrm->isActive());
      TEST_ASSERT(requiredPrm->isRequired());
      TEST_ASSERT(requiredPrm->isActive());
      TEST_ASSERT(!presentPrm->isRequired());
      TEST_ASSERT(presentPrm->isActive());

      #if 0
      if (verbose() > 0) {
         printEndl();
         absentPrm->writeParam(std::cout);
         requiredPrm->writeParam(std::cout);
         presentPrm->writeParam(std::cout);
      }
      #endif

      // Save to archive
      Serializable::OArchive oar;
      openOutputFile("out/binary", oar.file());
      Parameter::saveOptional(oar, absentVal, false);
      oar << requiredVal;
      Parameter::saveOptional(oar, presentVal, true);
      oar.file().close();

      // Load from archive
      double absentVal2;
      double requiredVal2;
      double presentVal2;
      Parameter* absentPrm2;
      Parameter* requiredPrm2;
      Parameter* presentPrm2;
      absentPrm2 = 
           new ScalarParam<double>("Absent", absentVal2, false);
      requiredPrm2 = 
           new ScalarParam<double>("Required", requiredVal2);
      presentPrm2 = 
           new ScalarParam<double>("Present", presentVal2, false);

      Serializable::IArchive iar;
      openInputFile("out/binary", iar.file());
      absentPrm2->load(iar);
      requiredPrm2->load(iar);
      presentPrm2->load(iar);
      iar.file().close();
      TEST_ASSERT(!absentPrm2->isRequired());
      TEST_ASSERT(!absentPrm2->isActive());
      TEST_ASSERT(requiredPrm2->isRequired());
      TEST_ASSERT(!presentPrm2->isRequired());
      TEST_ASSERT(presentPrm2->isActive());
      TEST_ASSERT(requiredVal2 == requiredVal);
      TEST_ASSERT(presentVal2 == presentVal);

      if (verbose() > 0) {
         printEndl();
         absentPrm2->writeParam(std::cout);
         requiredPrm2->writeParam(std::cout);
         presentPrm2->writeParam(std::cout);
      }

      delete absentPrm;
      delete requiredPrm;
      delete presentPrm;
      delete absentPrm2;
      delete requiredPrm2;
      delete presentPrm2;
   }

   void testParamStringWrite() {
      printMethod(TEST_FUNC);
      std::string requiredVal = "stringy";
      Parameter *requiredPrm;
      requiredPrm = new ScalarParam<std::string>("Required", requiredVal);
      if (verbose() > 0) {
         printEndl();
         requiredPrm->writeParam(std::cout);
      }
      delete requiredPrm;
   }

   void testParamStringRead1() {
      printMethod(TEST_FUNC);
      std::string requiredVal;
      Parameter *requiredPrm;
      requiredPrm = new ScalarParam<std::string>("Required", requiredVal);
      std::ifstream in;
      openInputFile("in/ScalarParamString", in);
      if (ParamComponent::echo()) std::cout << std::endl;
      requiredPrm->readParam(in);
      if (verbose() > 0) {
         printEndl();
         requiredPrm->writeParam(std::cout);
      }
      delete requiredPrm;
   }

   void testParamStringRead2() {
      printMethod(TEST_FUNC);
      int absentVal;
      std::string requiredVal;
      Parameter* absentPrm = new ScalarParam<int>("Absent", absentVal, false);
      Parameter* requiredPrm = new ScalarParam<std::string>("Required", requiredVal);
      std::ifstream in;
      openInputFile("in/ScalarParamString", in);
      if (ParamComponent::echo()) std::cout << std::endl;
      absentPrm->readParam(in);
      requiredPrm->readParam(in);
      if (verbose() > 0) {
         printEndl();
         requiredPrm->writeParam(std::cout);
      }
      delete absentPrm;
      delete requiredPrm;
   }







   void testParamStringReadSaveLoad()
   {
      printMethod(TEST_FUNC);
      std::string absentVal;
      std::string requiredVal;
      std::string presentVal;
      Parameter *absentPrm;
      Parameter *requiredPrm;
      Parameter *presentPrm;
      absentPrm = new ScalarParam<std::string>("Absent", absentVal, false);
      requiredPrm = new ScalarParam<std::string>("Required", requiredVal);
      presentPrm = new ScalarParam<std::string>("Present", presentVal, false);

      //Parameter::setEcho(true);
      if (ParamComponent::echo()) printEndl();

      std::ifstream in;
      openInputFile("in/ScalarParamString", in);
      absentPrm->readParam(in);
      TEST_ASSERT(!Label::isClear());
      requiredPrm->readParam(in);
      TEST_ASSERT(Label::isClear());
      presentPrm->readParam(in);
      TEST_ASSERT(Label::isClear());

      TEST_ASSERT(!absentPrm->isRequired());
      TEST_ASSERT(!absentPrm->isActive());
      TEST_ASSERT(requiredPrm->isRequired());
      TEST_ASSERT(requiredPrm->isActive());
      TEST_ASSERT(!presentPrm->isRequired());
      TEST_ASSERT(presentPrm->isActive());

      #if 0
      if (verbose() > 0) {
         printEndl();
         absentPrm->writeParam(std::cout);
         requiredPrm->writeParam(std::cout);
         presentPrm->writeParam(std::cout);
      }
      #endif

      // Save to archive
      Serializable::OArchive oar;
      openOutputFile("out/binary", oar.file());
      Parameter::saveOptional(oar, absentVal, false);
      oar << requiredVal;
      Parameter::saveOptional(oar, presentVal, true);
      oar.file().close();

      // Load from archive
      std::string absentVal2;
      std::string requiredVal2;
      std::string presentVal2;
      Parameter* absentPrm2;
      Parameter* requiredPrm2;
      Parameter* presentPrm2;
      absentPrm2 = 
           new ScalarParam<std::string>("Absent", absentVal2, false);
      requiredPrm2 = 
           new ScalarParam<std::string>("Required", requiredVal2);
      presentPrm2 = 
           new ScalarParam<std::string>("Present", presentVal2, false);

      Serializable::IArchive iar;
      openInputFile("out/binary", iar.file());
      absentPrm2->load(iar);
      requiredPrm2->load(iar);
      presentPrm2->load(iar);
      iar.file().close();
      TEST_ASSERT(!absentPrm2->isRequired());
      TEST_ASSERT(!absentPrm2->isActive());
      TEST_ASSERT(requiredPrm2->isRequired());
      TEST_ASSERT(!presentPrm2->isRequired());
      TEST_ASSERT(presentPrm2->isActive());

      TEST_ASSERT(requiredVal2 == requiredVal);
      TEST_ASSERT(presentVal2 == presentVal);
      TEST_ASSERT(requiredVal2 == "MyValue");
      TEST_ASSERT(presentVal2 == "HerValue");

      if (verbose() > 0) {
         printEndl();
         absentPrm2->writeParam(std::cout);
         requiredPrm2->writeParam(std::cout);
         presentPrm2->writeParam(std::cout);
      }

      delete absentPrm;
      delete requiredPrm;
      delete presentPrm;
      delete absentPrm2;
      delete requiredPrm2;
      delete presentPrm2;
   }

   void testCArrayParamIntWrite() 
   {
      printMethod(TEST_FUNC);
      int requiredVal[3];
      requiredVal[0] = 3;
      requiredVal[1] = 34;
      requiredVal[2] = 8;
      Parameter *requiredPrm;
      requiredPrm = new CArrayParam<int>("Required", requiredVal, 3);
      if (verbose() > 0) {
         printEndl();
         requiredPrm->writeParam(std::cout);
      }
      delete requiredPrm;
   }

   void testCArrayParamIntRead() {
      printMethod(TEST_FUNC);
      int requiredVal[3];
      Parameter *requiredPrm;
      requiredPrm = new CArrayParam<int>("Required", requiredVal, 3);
      std::ifstream in;
      openInputFile("in/ArrayParamInt", in);
      if (ParamComponent::echo()) std::cout << std::endl;
      requiredPrm->readParam(in);
      if (verbose() > 0) {
         printEndl();
         requiredPrm->writeParam(std::cout);
      }
      delete requiredPrm;
   }

   void testCArrayParamDoubleWrite() {
      printMethod(TEST_FUNC);
      double requiredVal[3];
      requiredVal[0] = 3.0;
      requiredVal[1] = 34.7;
      requiredVal[2] = 8.97296;
      Parameter *requiredPrm;
      requiredPrm = new CArrayParam<double>("Required", requiredVal, 3);
      if (verbose() > 0) {
         printEndl();
         requiredPrm->writeParam(std::cout);
      }
      delete requiredPrm;
   }

   void testCArrayParamDoubleRead() {
      printMethod(TEST_FUNC);
      double requiredVal[3];
      Parameter *requiredPrm;
      requiredPrm = new CArrayParam<double>("Required", requiredVal, 3);
      std::ifstream in;
      openInputFile("in/ArrayParamDouble", in);
      if (ParamComponent::echo()) std::cout << std::endl;
      requiredPrm->readParam(in);
      if (verbose() > 0) {
         printEndl();
         requiredPrm->writeParam(std::cout);
      }
      delete requiredPrm;
   }

   void testCArrayParamDoubleReadSaveLoad() {
      printMethod(TEST_FUNC);
      double absentVal[3];
      double requiredVal[3];
      double presentVal[3];
      Parameter *absentPrm;
      Parameter *requiredPrm;
      Parameter *presentPrm;
      absentPrm = new CArrayParam<double>("Absent", absentVal, 3, false);
      requiredPrm  = new CArrayParam<double>("Required", requiredVal, 3);
      presentPrm = new CArrayParam<double>("Present", presentVal, 3, false);
      std::ifstream in;
      openInputFile("in/ArrayParamDouble", in);
      if (ParamComponent::echo()) std::cout << std::endl;
      absentPrm->readParam(in);
      requiredPrm->readParam(in);
      presentPrm->readParam(in);
      TEST_ASSERT(!absentPrm->isRequired());
      TEST_ASSERT(!absentPrm->isActive());
      TEST_ASSERT(requiredPrm->isRequired());
      TEST_ASSERT(requiredPrm->isActive());
      TEST_ASSERT(!presentPrm->isRequired());
      TEST_ASSERT(presentPrm->isActive());

      // Save to archive
      Serializable::OArchive oar;
      openOutputFile("out/binary", oar.file());
      Parameter::saveOptionalCArray(oar, absentVal, 3, false);
      //requiredPrm->save(oar);
      oar.pack(requiredVal, 3);
      Parameter::saveOptionalCArray(oar, presentVal, 3, true);
      oar.file().close();

      // Load from archive
      double absentVal2[3];
      double requiredVal2[3];
      double presentVal2[3];
      Serializable::IArchive iar;
      openInputFile("out/binary", iar.file());
      Parameter* absentPrm2 = new CArrayParam<double>("Absent", absentVal2, 3, false);
      Parameter* requiredPrm2 = new CArrayParam<double>("Required", requiredVal2, 3);
      Parameter* presentPrm2 = new CArrayParam<double>("Present", presentVal2, 3, false);
      absentPrm2->load(iar);
      requiredPrm2->load(iar);
      presentPrm2->load(iar);
      iar.file().close();
      TEST_ASSERT(!absentPrm2->isRequired());
      TEST_ASSERT(!absentPrm2->isActive());
      TEST_ASSERT(requiredPrm2->label() == "Required");
      TEST_ASSERT(requiredPrm2->isRequired());
      TEST_ASSERT(presentPrm2->label() == "Present");
      TEST_ASSERT(!presentPrm2->isRequired());
      TEST_ASSERT(presentPrm2->isActive());

      if (verbose() > 0) {
         printEndl();
         absentPrm2->writeParam(std::cout);
         requiredPrm2->writeParam(std::cout);
         presentPrm2->writeParam(std::cout);
      }

      delete absentPrm;
      delete requiredPrm;
      delete presentPrm;
      delete absentPrm2;
      delete requiredPrm2;
      delete presentPrm2;
   }

   void testDArrayParamIntWrite() 
   {
      printMethod(TEST_FUNC);
      DArray<int> requiredVal;
      requiredVal.allocate(3);
      requiredVal[0] = 3;
      requiredVal[1] = 34;
      requiredVal[2] = 8;
      Parameter* requiredPrm = new DArrayParam<int>("Required", requiredVal, 3);
      if (verbose() > 0) {
         printEndl();
         requiredPrm->writeParam(std::cout);
      }
      delete requiredPrm;
   }

   void testDArrayParamIntRead() {
      printMethod(TEST_FUNC);
      DArray<int> requiredVal;
      requiredVal.allocate(3);
      Parameter* requiredPrm = new DArrayParam<int>("Required", requiredVal, 3);
      std::ifstream in;
      openInputFile("in/ArrayParamInt", in);
      if (ParamComponent::echo()) std::cout << std::endl;
      requiredPrm->readParam(in);
      if (verbose() > 0) {
         printEndl();
         requiredPrm->writeParam(std::cout);
      }
      delete requiredPrm;
   }

   void testDArrayParamDoubleWrite() {
      printMethod(TEST_FUNC);
      DArray<double> requiredVal;
      requiredVal.allocate(3);
      requiredVal[0] = 3.0;
      requiredVal[1] = 34.7;
      requiredVal[2] = 8.97296;
      Parameter *requiredPrm;
      requiredPrm = new DArrayParam<double>("Required", requiredVal, 3);
      if (verbose() > 0) {
         printEndl();
         requiredPrm->writeParam(std::cout);
      }
      delete requiredPrm;
   }

   void testDArrayParamDoubleRead() {
      printMethod(TEST_FUNC);
      DArray<double> requiredVal;
      requiredVal.allocate(3);
      Parameter *requiredPrm;
      requiredPrm = new DArrayParam<double>("Required", requiredVal, 3);
      std::ifstream in;
      openInputFile("in/ArrayParamDouble", in);
      if (ParamComponent::echo()) std::cout << std::endl;
      requiredPrm->readParam(in);
      in.close();
      if (verbose() > 0) {
         printEndl();
         requiredPrm->writeParam(std::cout);
      }
      delete requiredPrm;
   }

   void testDArrayParamDoubleReadSaveLoad()
   {
      printMethod(TEST_FUNC);
      DArray<double> absentVal;
      DArray<double> requiredVal;
      DArray<double> presentVal;
      absentVal.allocate(3);
      requiredVal.allocate(3);
      presentVal.allocate(3);
      Parameter *absentPrm;
      Parameter *requiredPrm;
      Parameter *presentPrm;
      absentPrm = new DArrayParam<double>("Absent", absentVal, 3, false);
      requiredPrm = new DArrayParam<double>("Required", requiredVal, 3);
      presentPrm = new DArrayParam<double>("Present", presentVal, 3, false);

      //Parameter::setEcho(true);
      if (ParamComponent::echo()) printEndl();

      std::ifstream in;
      openInputFile("in/ArrayParamDouble", in);
      absentPrm->readParam(in);
      TEST_ASSERT(!Label::isClear());
      requiredPrm->readParam(in);
      TEST_ASSERT(Label::isClear());
      presentPrm->readParam(in);
      TEST_ASSERT(Label::isClear());

      TEST_ASSERT(!absentPrm->isRequired());
      TEST_ASSERT(!absentPrm->isActive());
      TEST_ASSERT(requiredPrm->isRequired());
      TEST_ASSERT(requiredPrm->isActive());
      TEST_ASSERT(!presentPrm->isRequired());
      TEST_ASSERT(presentPrm->isActive());

      #if 0
      if (verbose() > 0) {
         printEndl();
         absentPrm->writeParam(std::cout);
         requiredPrm->writeParam(std::cout);
         presentPrm->writeParam(std::cout);
      }
      #endif

      // Save to archive
      Serializable::OArchive oar;
      openOutputFile("out/binary", oar.file());
      Parameter::saveOptional(oar, absentVal, false);
      oar << requiredVal;
      Parameter::saveOptional(oar, presentVal, true);
      oar.file().close();

      // Load from archive
      DArray<double> absentVal2;
      DArray<double> requiredVal2;
      DArray<double> presentVal2;
      absentVal2.allocate(3);
      requiredVal2.allocate(3);
      presentVal2.allocate(3);
      Parameter* absentPrm2;
      Parameter* requiredPrm2;
      Parameter* presentPrm2;
      absentPrm2 = 
           new DArrayParam<double>("Absent", absentVal2, 3, false);
      requiredPrm2 = 
           new DArrayParam<double>("Required", requiredVal2, 3);
      presentPrm2 = 
           new DArrayParam<double>("Present", presentVal2, 3, false);

      Serializable::IArchive iar;
      openInputFile("out/binary", iar.file());
      absentPrm2->load(iar);
      requiredPrm2->load(iar);
      presentPrm2->load(iar);
      iar.file().close();
      TEST_ASSERT(!absentPrm2->isRequired());
      TEST_ASSERT(!absentPrm2->isActive());
      TEST_ASSERT(requiredPrm2->isRequired());
      TEST_ASSERT(!presentPrm2->isRequired());
      TEST_ASSERT(presentPrm2->isActive());

      if (verbose() > 0) {
         printEndl();
         absentPrm2->writeParam(std::cout);
         requiredPrm2->writeParam(std::cout);
         presentPrm2->writeParam(std::cout);
      }

      delete absentPrm;
      delete requiredPrm;
      delete presentPrm;
      delete absentPrm2;
      delete requiredPrm2;
      delete presentPrm2;
   }

   void testFArrayParamIntWrite() 
   {
      printMethod(TEST_FUNC);
      FArray<int, 3> requiredVal;
      requiredVal[0] = 3;
      requiredVal[1] = 34;
      requiredVal[2] = 8;
      Parameter *requiredPrm;
      requiredPrm = new FArrayParam<int, 3>("Required", requiredVal);
      if (verbose() > 0) {
         printEndl();
         requiredPrm->writeParam(std::cout);
      }
      delete requiredPrm;
   }

   void testFArrayParamIntRead() {
      printMethod(TEST_FUNC);
      FArray<int, 3> requiredVal;
      Parameter *requiredPrm;
      requiredPrm = new FArrayParam<int, 3>("Required", requiredVal);
      std::ifstream in;
      openInputFile("in/ArrayParamInt", in);
      if (ParamComponent::echo()) std::cout << std::endl;
      requiredPrm->readParam(in);
      in.close();
      if (verbose() > 0) {
         printEndl();
         requiredPrm->writeParam(std::cout);
      }
      delete requiredPrm;
   }

   void testFArrayParamDoubleWrite() {
      printMethod(TEST_FUNC);
      FArray<double,3> requiredVal;
      requiredVal[0] = 3.0;
      requiredVal[1] = 34.7;
      requiredVal[2] = 8.97296;
      Parameter *requiredPrm;
      requiredPrm = new FArrayParam<double, 3>("Required", requiredVal);
      if (verbose() > 0) {
         printEndl();
         requiredPrm->writeParam(std::cout);
      }
      delete requiredPrm;
   }

   void testFArrayParamDoubleRead() {
      printMethod(TEST_FUNC);
      FArray<double,3> requiredVal;
      Parameter *requiredPrm;
      requiredPrm = new FArrayParam<double, 3>("Required", requiredVal);
      std::ifstream in;
      openInputFile("in/ArrayParamDouble", in);
      if (ParamComponent::echo()) std::cout << std::endl;
      requiredPrm->readParam(in);
      in.close();
      if (verbose() > 0) {
         printEndl();
         requiredPrm->writeParam(std::cout);
      }
      delete requiredPrm;
   }

   void testFArrayParamDoubleReadSaveLoad()
   {
      printMethod(TEST_FUNC);
      FArray<double, 3> absentVal;
      FArray<double, 3> requiredVal;
      FArray<double, 3> presentVal;
      // absentVal.allocate(3);
      // requiredVal.allocate(3);
      // presentVal.allocate(3);
      Parameter *absentPrm;
      Parameter *requiredPrm;
      Parameter *presentPrm;
      absentPrm = new FArrayParam<double, 3>("Absent", absentVal, false);
      requiredPrm = new FArrayParam<double, 3>("Required", requiredVal);
      presentPrm = new FArrayParam<double, 3>("Present", presentVal, false);

      //Parameter::setEcho(true);
      if (ParamComponent::echo()) printEndl();

      std::ifstream in;
      openInputFile("in/ArrayParamDouble", in);
      absentPrm->readParam(in);
      TEST_ASSERT(!Label::isClear());
      requiredPrm->readParam(in);
      TEST_ASSERT(Label::isClear());
      presentPrm->readParam(in);
      TEST_ASSERT(Label::isClear());

      TEST_ASSERT(!absentPrm->isRequired());
      TEST_ASSERT(!absentPrm->isActive());
      TEST_ASSERT(requiredPrm->isRequired());
      TEST_ASSERT(requiredPrm->isActive());
      TEST_ASSERT(!presentPrm->isRequired());
      TEST_ASSERT(presentPrm->isActive());

      #if 0
      if (verbose() > 0) {
         printEndl();
         absentPrm->writeParam(std::cout);
         requiredPrm->writeParam(std::cout);
         presentPrm->writeParam(std::cout);
      }
      #endif

      // Save to archive
      Serializable::OArchive oar;
      openOutputFile("out/binary", oar.file());
      Parameter::saveOptional(oar, absentVal, false);
      oar << requiredVal;
      Parameter::saveOptional(oar, presentVal, true);
      oar.file().close();

      // Load from archive
      FArray<double, 3> absentVal2;
      FArray<double, 3> requiredVal2;
      FArray<double, 3> presentVal2;
      //absentVal2.allocate(3);
      //requiredVal2.allocate(3);
      //presentVal2.allocate(3);
      Parameter* absentPrm2;
      Parameter* requiredPrm2;
      Parameter* presentPrm2;
      absentPrm2 = 
           new FArrayParam<double, 3>("Absent", absentVal2, false);
      requiredPrm2 = 
           new FArrayParam<double, 3>("Required", requiredVal2);
      presentPrm2 = 
           new FArrayParam<double, 3>("Present", presentVal2, false);

      Serializable::IArchive iar;
      openInputFile("out/binary", iar.file());
      absentPrm2->load(iar);
      requiredPrm2->load(iar);
      presentPrm2->load(iar);
      iar.file().close();
      TEST_ASSERT(!absentPrm2->isRequired());
      TEST_ASSERT(!absentPrm2->isActive());
      TEST_ASSERT(requiredPrm2->isRequired());
      TEST_ASSERT(!presentPrm2->isRequired());
      TEST_ASSERT(presentPrm2->isActive());

      if (verbose() > 0) {
         printEndl();
         absentPrm2->writeParam(std::cout);
         requiredPrm2->writeParam(std::cout);
         presentPrm2->writeParam(std::cout);
      }

      delete absentPrm;
      delete requiredPrm;
      delete presentPrm;
      delete absentPrm2;
      delete requiredPrm2;
      delete presentPrm2;
   }

   void testCArray2DParamDoubleWrite() {
      printMethod(TEST_FUNC);
      double requiredVal[2][2];
      requiredVal[0][0] = 3.0;
      requiredVal[0][1] = 34.7;
      requiredVal[1][0] = 8.97296;
      requiredVal[1][1] = 27.54;
      Parameter *requiredPrm;
      requiredPrm = 
           new CArray2DParam<double>("Required", &requiredVal[0][0], 2, 2, 2);
      if (verbose() > 0) {
         printEndl();
         requiredPrm->writeParam(std::cout);
      }
      delete requiredPrm;
   }

   void testCArray2DParamDoubleReadSaveLoad() 
   {
      printMethod(TEST_FUNC);
      double absentVal[2][2];
      double requiredVal[2][2];
      double presentVal[2][2];
      Parameter *absentPrm;
      Parameter *requiredPrm;
      Parameter *presentPrm;
      absentPrm = 
            new CArray2DParam<double>("Absent", absentVal[0], 2, 2, 2, false);
      requiredPrm  = 
            new CArray2DParam<double>("Required", requiredVal[0], 2, 2, 2);
      presentPrm = 
            new CArray2DParam<double>("Present", presentVal[0], 2, 2, 2, false);

      //Parameter::setEcho(true);
      if (ParamComponent::echo()) printEndl();

      std::ifstream in;
      openInputFile("in/MatrixParamDouble", in);
      absentPrm->readParam(in);
      TEST_ASSERT(!Label::isClear());
      requiredPrm->readParam(in);
      TEST_ASSERT(Label::isClear());
      presentPrm->readParam(in);
      TEST_ASSERT(Label::isClear());

      TEST_ASSERT(!absentPrm->isRequired());
      TEST_ASSERT(!absentPrm->isActive());
      TEST_ASSERT(requiredPrm->isRequired());
      TEST_ASSERT(requiredPrm->isActive());
      TEST_ASSERT(!presentPrm->isRequired());
      TEST_ASSERT(presentPrm->isActive());

      #if 0
      if (verbose() > 0) {
         printEndl();
         absentPrm->writeParam(std::cout);
         requiredPrm->writeParam(std::cout);
         presentPrm->writeParam(std::cout);
      }
      #endif

      // Save to archive
      Serializable::OArchive oar;
      openOutputFile("out/binary", oar.file());
      Parameter::saveOptionalCArray2D(oar, absentVal[0], 2, 2, 2, false);
      oar.pack(requiredVal[0], 2, 2, 2);
      Parameter::saveOptionalCArray2D(oar, presentVal[0], 2, 2, 2, true);
      oar.file().close();

      // Load from archive
      double absentVal2[2][2];
      double requiredVal2[2][2];
      double presentVal2[2][2];
      Parameter* absentPrm2;
      Parameter* requiredPrm2;
      Parameter* presentPrm2;
      absentPrm2 = 
           new CArray2DParam<double>("Absent", absentVal2[0], 2, 2, 2, false);
      requiredPrm2 = 
           new CArray2DParam<double>("Required", requiredVal2[0], 2, 2, 2);
      presentPrm2 = 
           new CArray2DParam<double>("Present", presentVal2[0], 2, 2, 2, false);

      Serializable::IArchive iar;
      openInputFile("out/binary", iar.file());
      absentPrm2->load(iar);
      requiredPrm2->load(iar);
      presentPrm2->load(iar);
      iar.file().close();
      TEST_ASSERT(!absentPrm2->isRequired());
      TEST_ASSERT(!absentPrm2->isActive());
      TEST_ASSERT(requiredPrm2->isRequired());
      TEST_ASSERT(!presentPrm2->isRequired());
      TEST_ASSERT(presentPrm2->isActive());

      if (verbose() > 0) {
         printEndl();
         absentPrm2->writeParam(std::cout);
         requiredPrm2->writeParam(std::cout);
         presentPrm2->writeParam(std::cout);
      }

      delete absentPrm;
      delete requiredPrm;
      delete presentPrm;
      delete absentPrm2;
      delete requiredPrm2;
      delete presentPrm2;
   }

   void testDMatrixParamDoubleWrite() {
      printMethod(TEST_FUNC);
      DMatrix<double> requiredVal;
      requiredVal.allocate(2, 2);
      requiredVal(0, 0) = 3.0;
      requiredVal(0, 1) = 34.7;
      requiredVal(1, 0) = 8.97296;
      requiredVal(1, 1) = 27.54;
      Parameter *requiredPrm;
      requiredPrm = new DMatrixParam<double>("Required", requiredVal, 2, 2);
      if (verbose() > -3) {
         printEndl();
         requiredPrm->writeParam(std::cout);
      }
      delete requiredPrm;
   }

   void testDMatrixParamDoubleRead() {
      printMethod(TEST_FUNC);
      DMatrix<double> requiredVal;
      requiredVal.allocate(2, 2);
      Parameter *requiredPrm;
      requiredPrm = new DMatrixParam<double>("Required", requiredVal, 2, 2);
      std::ifstream in;
      openInputFile("in/MatrixParamDouble", in);
      if (ParamComponent::echo()) std::cout << std::endl;
      requiredPrm->readParam(in);
      if (verbose() > 0) {
         printEndl();
         requiredPrm->writeParam(std::cout);
      }
      delete requiredPrm;
   }


   void testDMatrixParamDoubleReadSaveLoad() 
   {
      printMethod(TEST_FUNC);
      DMatrix<double> absentVal;
      DMatrix<double> requiredVal;
      DMatrix<double> presentVal;
      absentVal.allocate(2, 2);
      requiredVal.allocate(2, 2);
      presentVal.allocate(2, 2);
      Parameter *absentPrm;
      Parameter *requiredPrm;
      Parameter *presentPrm;
      absentPrm = 
            new DMatrixParam<double>("Absent", absentVal, 2, 2, false);
      requiredPrm  = 
            new DMatrixParam<double>("Required", requiredVal, 2, 2);
      presentPrm = 
            new DMatrixParam<double>("Present", presentVal, 2, 2, false);

      //Parameter::setEcho(true);
      if (ParamComponent::echo()) printEndl();

      std::ifstream in;
      openInputFile("in/MatrixParamDouble", in);
      absentPrm->readParam(in);
      TEST_ASSERT(!Label::isClear());
      requiredPrm->readParam(in);
      TEST_ASSERT(Label::isClear());
      presentPrm->readParam(in);
      TEST_ASSERT(Label::isClear());

      TEST_ASSERT(!absentPrm->isRequired());
      TEST_ASSERT(!absentPrm->isActive());
      TEST_ASSERT(requiredPrm->isRequired());
      TEST_ASSERT(requiredPrm->isActive());
      TEST_ASSERT(!presentPrm->isRequired());
      TEST_ASSERT(presentPrm->isActive());

      #if 0
      if (verbose() > 0) {
         printEndl();
         absentPrm->writeParam(std::cout);
         requiredPrm->writeParam(std::cout);
         presentPrm->writeParam(std::cout);
      }
      #endif

      // Save to archive
      Serializable::OArchive oar;
      openOutputFile("out/binary", oar.file());
      Parameter::saveOptional(oar, absentVal, false);
      oar << requiredVal;
      Parameter::saveOptional(oar, presentVal, true);
      oar.file().close();

      // Load from archive
      DMatrix<double> absentVal2;
      DMatrix<double> requiredVal2;
      DMatrix<double> presentVal2;
      absentVal2.allocate(2, 2);
      requiredVal2.allocate(2, 2);
      presentVal2.allocate(2, 2);
      Parameter* absentPrm2;
      Parameter* requiredPrm2;
      Parameter* presentPrm2;
      absentPrm2 = 
           new DMatrixParam<double>("Absent", absentVal2, 2, 2, false);
      requiredPrm2 = 
           new DMatrixParam<double>("Required", requiredVal2, 2, 2);
      presentPrm2 = 
           new DMatrixParam<double>("Present", presentVal2, 2, 2, false);

      Serializable::IArchive iar;
      openInputFile("out/binary", iar.file());
      absentPrm2->load(iar);
      requiredPrm2->load(iar);
      presentPrm2->load(iar);
      iar.file().close();
      TEST_ASSERT(!absentPrm2->isRequired());
      TEST_ASSERT(!absentPrm2->isActive());
      TEST_ASSERT(requiredPrm2->isRequired());
      TEST_ASSERT(!presentPrm2->isRequired());
      TEST_ASSERT(presentPrm2->isActive());

      if (verbose() > 0) {
         printEndl();
         absentPrm2->writeParam(std::cout);
         requiredPrm2->writeParam(std::cout);
         presentPrm2->writeParam(std::cout);
      }

      delete absentPrm;
      delete requiredPrm;
      delete presentPrm;
      delete absentPrm2;
      delete requiredPrm2;
      delete presentPrm2;
   }

   void testDSymmMatrixParamDoubleWrite() {
      printMethod(TEST_FUNC);
      DMatrix<double> requiredVal;
      requiredVal.allocate(2, 2);
      requiredVal(0, 0) = 3.0;
      requiredVal(0, 1) = 34.7;
      requiredVal(1, 0) = 34.7;
      requiredVal(1, 1) = 27.54;
      Parameter *requiredPrm;
      requiredPrm = new DSymmMatrixParam<double>("Required", requiredVal, 2);
      if (verbose() > -3) {
         printEndl();
         requiredPrm->writeParam(std::cout);
      }
      delete requiredPrm;
   }

   void testDSymmMatrixParamDoubleRead() {
      printMethod(TEST_FUNC);
      DMatrix<double> requiredVal;
      requiredVal.allocate(2, 2);
      Parameter *requiredPrm;
      requiredPrm = new DSymmMatrixParam<double>("Required", requiredVal, 2);
      std::ifstream in;
      openInputFile("in/SymmMatrixParamDouble", in);
      if (ParamComponent::echo()) std::cout << std::endl;
      requiredPrm->readParam(in);
      if (verbose() > -1) {
         printEndl();
         requiredPrm->writeParam(std::cout);
      }
      delete requiredPrm;
   }

   void testDSymmMatrixParamDoubleReadSaveLoad() 
   {
      printMethod(TEST_FUNC);
      DMatrix<double> absentVal;
      DMatrix<double> requiredVal;
      DMatrix<double> presentVal;
      absentVal.allocate(2, 2);
      requiredVal.allocate(2, 2);
      presentVal.allocate(2, 2);
      Parameter *absentPrm;
      Parameter *requiredPrm;
      Parameter *presentPrm;
      absentPrm = 
            new DSymmMatrixParam<double>("Absent", absentVal, 2, false);
      requiredPrm  = 
            new DSymmMatrixParam<double>("Required", requiredVal, 2);
      presentPrm = 
            new DSymmMatrixParam<double>("Present", presentVal, 2, false);

      //Parameter::setEcho(true);
      if (ParamComponent::echo()) printEndl();

      std::ifstream in;
      openInputFile("in/SymmMatrixParamDouble", in);
      absentPrm->readParam(in);
      TEST_ASSERT(!Label::isClear());
      requiredPrm->readParam(in);
      TEST_ASSERT(Label::isClear());
      presentPrm->readParam(in);
      TEST_ASSERT(Label::isClear());

      TEST_ASSERT(!absentPrm->isRequired());
      TEST_ASSERT(!absentPrm->isActive());
      TEST_ASSERT(requiredPrm->isRequired());
      TEST_ASSERT(requiredPrm->isActive());
      TEST_ASSERT(!presentPrm->isRequired());
      TEST_ASSERT(presentPrm->isActive());

      #if 0
      if (verbose() > 0) {
         printEndl();
         absentPrm->writeParam(std::cout);
         requiredPrm->writeParam(std::cout);
         presentPrm->writeParam(std::cout);
      }
      #endif

      // Save to archive
      Serializable::OArchive oar;
      openOutputFile("out/binary", oar.file());
      Parameter::saveOptional(oar, absentVal, false);
      oar << requiredVal;
      Parameter::saveOptional(oar, presentVal, true);
      oar.file().close();

      // Load from archive
      DMatrix<double> absentVal2;
      DMatrix<double> requiredVal2;
      DMatrix<double> presentVal2;
      absentVal2.allocate(2, 2);
      requiredVal2.allocate(2, 2);
      presentVal2.allocate(2, 2);
      Parameter* absentPrm2;
      Parameter* requiredPrm2;
      Parameter* presentPrm2;
      absentPrm2 = 
           new DSymmMatrixParam<double>("Absent", absentVal2, 2, false);
      requiredPrm2 = 
           new DSymmMatrixParam<double>("Required", requiredVal2, 2);
      presentPrm2 = 
           new DSymmMatrixParam<double>("Present", presentVal2, 2, false);

      Serializable::IArchive iar;
      openInputFile("out/binary", iar.file());
      absentPrm2->load(iar);
      requiredPrm2->load(iar);
      presentPrm2->load(iar);
      iar.file().close();
      TEST_ASSERT(!absentPrm2->isRequired());
      TEST_ASSERT(!absentPrm2->isActive());
      TEST_ASSERT(requiredPrm2->isRequired());
      TEST_ASSERT(!presentPrm2->isRequired());
      TEST_ASSERT(presentPrm2->isActive());

      if (verbose() > 0) {
         printEndl();
         absentPrm2->writeParam(std::cout);
         requiredPrm2->writeParam(std::cout);
         presentPrm2->writeParam(std::cout);
      }

      delete absentPrm;
      delete requiredPrm;
      delete presentPrm;
      delete absentPrm2;
      delete requiredPrm2;
      delete presentPrm2;
   }
   #if 0
   #endif

};

TEST_BEGIN(ParameterTest)
TEST_ADD(ParameterTest, testParamIntConstructor1)
TEST_ADD(ParameterTest, testParamIntConstructor2)
TEST_ADD(ParameterTest, testParamIntWrite)
TEST_ADD(ParameterTest, testParamIntRead1)
TEST_ADD(ParameterTest, testParamIntRead2)
TEST_ADD(ParameterTest, testParamIntReadSaveLoad)
TEST_ADD(ParameterTest, testParamDoubleWrite)
TEST_ADD(ParameterTest, testParamDoubleRead)
TEST_ADD(ParameterTest, testParamDoubleReadSaveLoad)
TEST_ADD(ParameterTest, testParamStringWrite)
TEST_ADD(ParameterTest, testParamStringRead1)
TEST_ADD(ParameterTest, testParamStringRead2)
TEST_ADD(ParameterTest, testParamStringReadSaveLoad)
TEST_ADD(ParameterTest, testCArrayParamIntWrite)
TEST_ADD(ParameterTest, testCArrayParamIntRead)
TEST_ADD(ParameterTest, testCArrayParamDoubleWrite)
TEST_ADD(ParameterTest, testCArrayParamDoubleRead)
TEST_ADD(ParameterTest, testCArrayParamDoubleReadSaveLoad)
TEST_ADD(ParameterTest, testDArrayParamIntWrite)
TEST_ADD(ParameterTest, testDArrayParamIntRead)
TEST_ADD(ParameterTest, testDArrayParamDoubleWrite)
TEST_ADD(ParameterTest, testDArrayParamDoubleRead)
TEST_ADD(ParameterTest, testDArrayParamDoubleReadSaveLoad)
TEST_ADD(ParameterTest, testFArrayParamIntWrite)
TEST_ADD(ParameterTest, testFArrayParamIntRead)
TEST_ADD(ParameterTest, testFArrayParamDoubleWrite)
TEST_ADD(ParameterTest, testFArrayParamDoubleRead)
TEST_ADD(ParameterTest, testFArrayParamDoubleReadSaveLoad)
TEST_ADD(ParameterTest, testCArray2DParamDoubleWrite)
TEST_ADD(ParameterTest, testCArray2DParamDoubleReadSaveLoad)
TEST_ADD(ParameterTest, testDMatrixParamDoubleWrite)
TEST_ADD(ParameterTest, testDMatrixParamDoubleRead)
TEST_ADD(ParameterTest, testDMatrixParamDoubleReadSaveLoad)
TEST_ADD(ParameterTest, testDSymmMatrixParamDoubleWrite)
TEST_ADD(ParameterTest, testDSymmMatrixParamDoubleRead)
TEST_ADD(ParameterTest, testDSymmMatrixParamDoubleReadSaveLoad)
TEST_END(ParameterTest)

#endif

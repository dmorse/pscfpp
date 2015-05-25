#ifndef PARAM_FILE_TEST_H
#define PARAM_FILE_TEST_H

#include "UnitTest.h"

#include <string>
#include <iostream>
#include <fstream>

/**
* A UnitTest with a built-in input file.
*
* \ingroup Test_Module
**/
class ParamFileTest : public UnitTest 
{

public:

   /**
   * Constructor.
   */
   ParamFileTest()
    : UnitTest()
   {}

   /**
   * Destructor.
   */
   ~ParamFileTest()
   { closeFile(); }

   /**
   * Close the input file.
   */
   virtual void tearDown()
   {  closeFile(); }

   /**
   * Open the input file.
   */
   void openFile(const char *fileName)
   {
      if (isIoProcessor()) {
         openInputFile(std::string(fileName), file_);
      }
   }

   /**
   * Close the input file.
   */
   void closeFile()
   {  if (file_.is_open()) file_.close(); }

   /**
   * Returns input file by reference.
   */
   std::ifstream& file()
   {  return file_; }

private:

   std::ifstream  file_;

};

#endif

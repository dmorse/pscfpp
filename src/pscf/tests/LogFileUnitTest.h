#ifndef PSCF_LOG_FILE_UNIT_TEST_H
#define PSCF_LOG_FILE_UNIT_TEST_H

#include <test/UnitTest.h>
#include <util/misc/Log.h>

#include <fstream>

class LogFileUnitTest : public UnitTest
{

public:

   void tearDown()
   {  closeLogFile(); }

protected:

   std::ofstream logFile_;

   void openLogFile(char const * filename)
   {
      openOutputFile(filename, logFile_);
      Util::Log::setFile(logFile_);
   }

   void closeLogFile()
   {
      if (logFile_.is_open()) {
         Util::Log::close();
      }
   }


};
#endif

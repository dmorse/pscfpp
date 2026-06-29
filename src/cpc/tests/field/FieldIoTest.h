#ifndef CPC_FIELD_IO_TEST_H
#define CPC_FIELD_IO_TEST_H

#include <test/UnitTest.h>
#include <test/UnitTestRunner.h>

#include <cpc/field/Domain.h>
#include <cpc/field/FieldIo.h>

#include <prdc/field/cpu/RField.h>
#include <prdc/field/cpu/RFieldDft.h>
#include <prdc/field/cpu/RFieldComparison.h>
#include <prdc/field/cpu/RFieldDftComparison.h>
#include <prdc/field/cpu/FFT.h>
#include <prdc/fieldIo/fieldCheck.h>
#include <prdc/crystal/BFieldComparison.h>
#include <prdc/crystal/Basis.h>
#include <prdc/crystal/UnitCell.h>

#include <pscf/mesh/Mesh.h>
#include <pscf/mesh/MeshIterator.h>

#include <util/containers/DArray.h>
#include <util/misc/FileMaster.h>
#include <util/format/Dbl.h>

#include <iostream>
#include <fstream>

using namespace Util;
using namespace Pscf;
using namespace Pscf::Cpc;
using namespace Pscf::Prdc;
using namespace Pscf::Prdc::Cpu;

class FieldIoTest : public UnitTest 
{

   std::ofstream logFile_;
   FileMaster fileMaster_;
   int nMonomer_;

public:

   void setUp()
   {
      setVerbose(0);
      nMonomer_ = 2;
      openLogFile("out/fieldIoTest.log");
      fileMaster_.setRootPrefix(UnitTest::filePrefix());
   }

   void tearDown()
   {
      if (logFile_.is_open()) {
         logFile_.close();
      }
   }

   void openLogFile(char const * filename)
   {
      openOutputFile(filename, logFile_);
      Log::setFile(logFile_);
   }

   /*
   * Open and read parameter header to initialize Domain<D> system.
   */
   template <int D>
   void readParam(std::string filename, Domain<D>& domain)
   {
      std::ifstream in;
      openInputFile(filename, in);
      domain.readParam(in);
      domain.setFileMaster(fileMaster_);
      in.close();
   }

   /*
   * Open and read file header to initialize Domain<D> system.
   */
   template <int D>
   void readHeader(std::string filename, Domain<D>& domain)
   {
      std::ifstream in;
      openInputFile(filename, in);
      domain.readFieldHeader(in, nMonomer_);
      domain.setFileMaster(fileMaster_);
      domain.fieldIo().setNMonomer(nMonomer_);
      in.close();
   }

   template <int D>
   void readFields(std::string filename, 
                   Domain<D>& domain,
                   DArray< CField<D> >& fields)
   {
      std::ifstream in;
      openInputFile(filename, in);
      domain.fieldIo().readFields(in, fields, domain.unitCell());
      in.close();
   }

   template <int D>
   void writeFields(std::string filename, Domain<D>& domain,
                   DArray< CField<D> > const & fields)
   {
      std::ofstream out;
      openOutputFile(filename, out);
      domain.fieldIo().writeFields(out, fields, domain.unitCell());
      out.close();
   }

   void testReadHeader() 
   {
      printMethod(TEST_FUNC);

      Domain<3> domain;
      domain.setFileMaster(fileMaster_);
      readHeader("in/w_bcc.rf", domain);

      TEST_ASSERT(domain.mesh().dimension(0) == 32);
      TEST_ASSERT(domain.mesh().dimension(1) == 32);
      TEST_ASSERT(domain.mesh().dimension(2) == 32);
      TEST_ASSERT(domain.unitCell().lattice() == UnitCell<3>::Cubic);
      //TEST_ASSERT(nMonomer_ == 2);

      if (verbose() > 0) {
         std::cout << "\n";
         std::cout << "Cell  = " << domain.unitCell() << "\n";
         std::cout << "Ngrid = " << domain.mesh().dimensions() << "\n";
      }

      DArray< CField<3> >  fc;
      allocateFields(fc, nMonomer_, domain.mesh().dimensions());

   }

   void testReadWrite() 
   {
      printMethod(TEST_FUNC);

      Domain<3> domain;
      domain.setFileMaster(fileMaster_);
      readHeader("in/w_bcc.rf", domain);

      TEST_ASSERT(domain.mesh().dimension(0) == 32);
      TEST_ASSERT(domain.mesh().dimension(1) == 32);
      TEST_ASSERT(domain.mesh().dimension(2) == 32);
      TEST_ASSERT(domain.unitCell().lattice() == UnitCell<3>::Cubic);
      //TEST_ASSERT(nMonomer_ == 2);

      if (verbose() > 0) {
         std::cout << "\n";
         std::cout << "Cell  = " << domain.unitCell() << "\n";
         std::cout << "Ngrid = " << domain.mesh().dimensions() << "\n";
      }

      // Write arbitrary values to file
      DArray< CField<3> >  f0;
      allocateFields(f0, nMonomer_, domain.mesh().dimensions());
      double x, y;
      int rank;
      double meshSize = (double) domain.mesh().size();
      MeshIterator<3> iter(domain.mesh().dimensions());
      for (int i=0; i < nMonomer_; ++i) {
         CField<3>& field = f0[i];
         for (iter.begin(); !iter.atEnd(); ++iter) {
            rank = iter.rank();
            x = (double) rank/ meshSize;
            x = x - 1.0;
            y = 0.5 + x - x*x*x;
            x += (double) (i+1);
            y -= (double) (i+1);
            assign(field[rank], x, y);
         }
      }
      writeFields("out/random3D.cf", domain, f0);

      DArray< CField<3> >  f1;
      allocateFields(f1, nMonomer_, domain.mesh().dimensions());
      readFields("out/random3D.cf", domain, f1);
      for (int i=0; i < nMonomer_; ++i) {
         CField<3>& field0= f0[i];
         CField<3>& field1= f1[i];
         for (iter.begin(); !iter.atEnd(); ++iter) {
            rank = iter.rank();
            TEST_ASSERT(eq(field0[rank][0], field1[rank][0]));
            TEST_ASSERT(eq(field0[rank][1], field1[rank][1]));
         }
      }

   }

};

TEST_BEGIN(FieldIoTest)
TEST_ADD(FieldIoTest, testReadHeader)
TEST_ADD(FieldIoTest, testReadWrite)
TEST_END(FieldIoTest)

#endif

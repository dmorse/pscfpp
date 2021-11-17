#ifndef PSPC_FIELD_IO_TEST_H
#define PSPC_FIELD_IO_TEST_H

#include <test/UnitTest.h>
#include <test/UnitTestRunner.h>

#include <pspc/field/FieldIo.h>
#include <pspc/field/RField.h>
#include <pspc/field/RFieldDft.h>
#include <pspc/field/FFT.h>

#include <pscf/crystal/Basis.h>
#include <pscf/crystal/UnitCell.h>
#include <pscf/mesh/Mesh.h>
#include <pscf/mesh/MeshIterator.h>

#include <util/containers/DArray.h>
#include <util/math/Constants.h>
#include <util/misc/FileMaster.h>
#include <util/format/Dbl.h>

#include <iostream>
#include <fstream>

using namespace Util;
using namespace Pscf;
using namespace Pscf::Pspc;

class FieldIoTest : public UnitTest 
{

   std::ofstream logFile_;

public:

   /*
   * A subsystem containing a FieldIo and associated components.
   */
   template <int D>
   class FieldIoSystem
   {
   public:

      UnitCell<D> unitCell;
      Mesh<D> mesh;
      IntVec<D> nGrid;
      FFT<D> fft;
      Basis<D> basis;
      FileMaster fileMaster;
      FieldIo<D> fieldIo;
      std::string groupName;
      int nMonomer;

      // Constructor
      FieldIoSystem()
      {  fieldIo.associate(mesh, fft, groupName, basis, fileMaster); }
     
      // Read field file header 
      void initialize(std::istream& in)
      {
         // Read standard field header
         int ver1, ver2;
         Pscf::readFieldHeader(in, ver1, ver2, 
                         unitCell, groupName, nMonomer);
 
         // Read grid dimensions 
         std::string label;
         in >> label;
         UTIL_CHECK(label == "ngrid");
         IntVec<D> nGrid;
         in >> nGrid;

         // Initialize mesh, fft and basis
         mesh.setDimensions(nGrid);
         fft.setup(mesh.dimensions());
         basis.makeBasis(mesh, unitCell, groupName);
      }

   };

   void setUp()
   {
      setVerbose(0);
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
   * Open and read file header to initialize FieldIoSystem<D> system.
   */
   template <int D>
   void readHeader(std::string filename, FieldIoSystem<D>& system)
   {
      std::ifstream in;
      openInputFile(filename, in);
      system.initialize(in);
      in.close();
   }

   // Allocate an array of fields in symmetry adapated format
   void allocateFieldsBasis(int nMonomer, int nStar,
                            DArray< DArray<double> >& fields)
   {
      fields.allocate(nMonomer);
      for (int i = 0; i < nMonomer; ++i) {   
         fields[i].allocate(nStar);
      }
   }

   // Allocate an array of r-grid fields
   template <int D>
   void allocateFieldsRGrid(int nMonomer, IntVec<D> dimensions,
                            DArray< RField<D> >& fields)
   {
      fields.allocate(nMonomer);
      for (int i = 0; i < nMonomer; ++i) {   
         fields[i].allocate(dimensions);
      }
   }

   // Allocate an array of k-grid fields
   template <int D>
   void allocateFieldsKGrid(int nMonomer, IntVec<D> dimensions,
                            DArray< RFieldDft<D> >& fields)
   {
      fields.allocate(nMonomer);
      for (int i = 0; i < nMonomer; ++i) {
         fields[i].allocate(dimensions);
      }
   }

   void testFieldIoSystem3D() 
   {
      printMethod(TEST_FUNC);

      FieldIoSystem<3> sys;
      readHeader<3>("in/w_bcc.rf", sys);

      TEST_ASSERT(sys.mesh.dimension(0) == 32);
      TEST_ASSERT(sys.mesh.dimension(1) == 32);
      TEST_ASSERT(sys.mesh.dimension(2) == 32);
      TEST_ASSERT(sys.unitCell.lattice() == UnitCell<3>::Cubic);
      TEST_ASSERT(sys.basis.nBasis() == 489);

      if (verbose() > 0) {
         std::cout << "\n";
         std::cout << "Cell  = " << sys.unitCell << "\n";
         std::cout << "Ngrid = " << sys.mesh.dimensions() << "\n";
         if (verbose() > 1) {
            sys.basis.outputStars(std::cout);
         }
      }

      DArray< DArray<double> > fieldsB;
      allocateFieldsBasis(sys.nMonomer, sys.basis.nStar(), fieldsB);

      DArray< RField<3> >  fieldsR;
      allocateFieldsRGrid(sys.nMonomer, sys.mesh.dimensions(), fieldsR);

      DArray< RFieldDft<3> > fieldsK;
      allocateFieldsKGrid(sys.nMonomer, sys.mesh.dimensions(), fieldsK);

   }

};

TEST_BEGIN(FieldIoTest)
TEST_ADD(FieldIoTest, testFieldIoSystem3D)
TEST_END(FieldIoTest)

#endif

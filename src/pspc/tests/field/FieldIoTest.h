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

public:

   void setUp()
   {
      setVerbose(2);
   }

   void tearDown()
   {}


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

   template <int D>
   void readHeader(std::string filename, FieldIoSystem<D>& system)
   {
      std::ifstream in;
      openInputFile(filename, in);
      system.initialize(in);
      in.close();
   }

   void testSystem3D() 
   {
      printMethod(TEST_FUNC);

      FieldIoSystem<3> system;
      readHeader<3>("in/w_bcc.rf", system);

      TEST_ASSERT(system.mesh.dimension(0) == 32);
      TEST_ASSERT(system.mesh.dimension(1) == 32);
      TEST_ASSERT(system.mesh.dimension(2) == 32);
      TEST_ASSERT(system.unitCell.lattice() == UnitCell<3>::Cubic);
      TEST_ASSERT(system.basis.nBasis() == 489);

      //std::cout << std::endl;
      //std::cout << "Cell  = " << system.unitCell << std::endl;
      //std::cout << "Ngrid = " << system.mesh.dimensions() << std::endl;
      //system.basis.outputStars(std::cout);
   }

};

TEST_BEGIN(FieldIoTest)
TEST_ADD(FieldIoTest, testSystem3D)
TEST_END(FieldIoTest)

#endif

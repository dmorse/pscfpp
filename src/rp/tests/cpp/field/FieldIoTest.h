#ifndef RPC_FIELD_IO_TEST_H
#define RPC_FIELD_IO_TEST_H

#include <test/UnitTest.h>
#include <test/UnitTestRunner.h>

#include <rp/field/Domain.h>
#include <rp/field/FieldIo.h>

#include <pscf/field/cpp/RField.h>
#include <pscf/field/cpp/RFieldDft.h>
#include <pscf/field/cpp/RFieldComparison.h>
#include <pscf/field/cpp/RFieldDftComparison.h>
#include <pscf/field/cpp/FFT.h>
#include <pscf/crystal/BFieldComparison.h>
#include <pscf/crystal/Basis.h>
#include <pscf/crystal/UnitCell.h>

#include <pscf/mesh/Mesh.h>
#include <pscf/mesh/MeshIterator.h>
#include <pscf/interaction/Interaction.h>

#include <util/containers/DArray.h>
#include <util/containers/FSArray.h>
#include <util/misc/FileMaster.h>
#include <util/misc/Exception.h>
#include <util/format/Dbl.h>

#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
#include <map>
#include <cmath>

using namespace Util;
using namespace Pscf;

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
   * Open and read parameter header to initialize Rp::Domain<D,CPT> system.
   */
   template <int D>
   void readParam(std::string filename, Rp::Domain<D,CPT>& domain)
   {
      std::ifstream in;
      openInputFile(filename, in);
      domain.readParam(in);
      in.close();
   }

   /*
   * Open and read file header to initialize Rp::Domain<D,CPT> system.
   */
   template <int D>
   void readHeader(std::string filename, Rp::Domain<D,CPT>& domain)
   {
      std::ifstream in;
      openInputFile(filename, in);
      domain.readRGridFieldHeader(in, nMonomer_);
      in.close();
   }

   /*
   * Initialize Rp::Domain<D,CPT> system from the text of a file header.
   */
   template <int D>
   void readHeaderString(std::string header, Rp::Domain<D,CPT>& domain)
   {
      std::istringstream in(header);
      domain.readRGridFieldHeader(in, nMonomer_);
   }

   // Allocate an array of fields in symmetry adapated format
   void allocateFields(int nMonomer, int nStar,
                            DArray< DArray<double> >& fields)
   {
      fields.allocate(nMonomer);
      for (int i = 0; i < nMonomer; ++i) {   
         fields[i].allocate(nStar);
      }
   }

   // Allocate an array of r-grid fields
   template <int D>
   void allocateFields(int nMonomer, IntVec<D> dimensions,
                            DArray< RField<D,CPT> >& fields)
   {
      fields.allocate(nMonomer);
      for (int i = 0; i < nMonomer; ++i) {   
         fields[i].allocate(dimensions);
      }
   }

   // Allocate an array of k-grid fields
   template <int D>
   void allocateFields(int nMonomer, IntVec<D> dimensions,
                            DArray< RFieldDft<D,CPT> >& fields)
   {
      fields.allocate(nMonomer);
      for (int i = 0; i < nMonomer; ++i) {
         fields[i].allocate(dimensions);
      }
   }

   template <int D>
   void readFields(std::string filename, Rp::Domain<D,CPT>& domain,
                   DArray< DArray<double> >& fields)
   {
      std::ifstream in;
      openInputFile(filename, in);
      domain.fieldIo().readFieldsBasis(in, fields, domain.unitCell());
      in.close();
   }

   template <int D>
   void readFields(std::string filename, Rp::Domain<D,CPT>& domain,
                   DArray< RField<D,CPT> >& fields)
   {
      std::ifstream in;
      openInputFile(filename, in);
      domain.fieldIo().readFieldsRGrid(in, fields, domain.unitCell());
      in.close();
   }

   template <int D>
   void readFields(std::string filename, Rp::Domain<D,CPT>& domain,
                   DArray< RFieldDft<D,CPT> >& fields)
   {
      std::ifstream in;
      openInputFile(filename, in);
      domain.fieldIo().readFieldsKGrid(in, fields, domain.unitCell());
      in.close();
   }

   template <int D>
   void writeFields(std::string filename, Rp::Domain<D,CPT>& domain,
                   DArray< DArray<double> > const & fields)
   {
      std::ofstream out;
      openOutputFile(filename, out);
      domain.fieldIo().writeFieldsBasis(out, fields, domain.unitCell());
      out.close();
   }

   template <int D>
   void writeFields(std::string filename, Rp::Domain<D,CPT>& domain,
                   DArray< RField<D,CPT> > const & fields)
   {
      std::ofstream out;
      openOutputFile(filename, out);
      domain.fieldIo().writeFieldsRGrid(out, fields, domain.unitCell());
      out.close();
   }

   template <int D>
   void writeFields(std::string filename, Rp::Domain<D,CPT>& domain,
                   DArray< RFieldDft<D,CPT> > const & fields)
   {
      std::ofstream out;
      openOutputFile(filename, out);
      domain.fieldIo().writeFieldsKGrid(out, fields, domain.unitCell());
      out.close();
   }

   void testReadHeader() 
   {
      printMethod(TEST_FUNC);

      Rp::Domain<3,CPT> domain;
      domain.setFileMaster(fileMaster_);
      readHeader("in/w_bcc.rf", domain);

      TEST_ASSERT(domain.mesh().dimension(0) == 32);
      TEST_ASSERT(domain.mesh().dimension(1) == 32);
      TEST_ASSERT(domain.mesh().dimension(2) == 32);
      TEST_ASSERT(domain.unitCell().lattice() == UnitCell<3>::Cubic);
      TEST_ASSERT(domain.basis().nBasis() == 489);
      //TEST_ASSERT(nMonomer_ == 2);

      if (verbose() > 0) {
         std::cout << "\n";
         std::cout << "Cell  = " << domain.unitCell() << "\n";
         std::cout << "Ngrid = " << domain.mesh().dimensions() << "\n";
         if (verbose() > 1) {
            domain.basis().outputStars(std::cout);
         }
      }

      DArray< DArray<double> > fb;
      allocateFields(nMonomer_, domain.basis().nBasis(), fb);

      DArray< RField<3,CPT> >  fr;
      allocateFields(nMonomer_, domain.mesh().dimensions(), fr);

      DArray< RFieldDft<3,CPT> > fk;
      allocateFields(nMonomer_, domain.mesh().dimensions(), fk);

   }

   void testBasisIo_lam() 
   {
      printMethod(TEST_FUNC);

      // Create and initialize domain
      Rp::Domain<1,CPT> domain;
      domain.setFileMaster(fileMaster_);
      readHeader("in/w_lam.rf", domain);

      // Create required containers and stream
      DArray< DArray<double> > bf_0;
      allocateFields(nMonomer_, domain.basis().nBasis(), bf_0);
      DArray< DArray<double> > bf_1;
      allocateFields(nMonomer_, domain.basis().nBasis(), bf_1);
      std::ifstream in;
      std::ofstream out;

      // Read in/w_lam.bf into bf_0 in basis format
      openInputFile("in/w_lam.bf", in);
      domain.fieldIo().readFieldsBasis(in, bf_0, domain.unitCell());
      in.close();

      // Write bf_0 to out/w_lam.bf in basis format
      openOutputFile("out/w_lam.bf", out);
      domain.fieldIo().writeFieldsBasis(out, bf_0, domain.unitCell());
      out.close();

      // Read out/w_lam.bf to bf_1 in basis format
      openInputFile("out/w_lam.bf", in);
      domain.fieldIo().readFieldsBasis(in, bf_1, domain.unitCell());
      in.close();

      // Compare initial and final fields in basis format
      BFieldComparison comparison;
      comparison.compare(bf_0, bf_1);
      if (verbose() > 0) {
         std::cout  << std::endl;
         std::cout  << Dbl(comparison.maxDiff(), 21, 13) << std::endl;
         std::cout  << Dbl(comparison.rmsDiff(), 21, 13) << std::endl;
      }
      TEST_ASSERT(comparison.maxDiff() < 1.0E-12);

   }

   void testBasisIo_hex() 
   {
      printMethod(TEST_FUNC);

      // Create and initialize domain
      Rp::Domain<2,CPT> domain;
      domain.setFileMaster(fileMaster_);
      readHeader("in/w_hex.rf", domain);

      // Create required containers and stream
      DArray< DArray<double> > bf_0;
      allocateFields(nMonomer_, domain.basis().nBasis(), bf_0);
      DArray< DArray<double> > bf_1;
      allocateFields(nMonomer_, domain.basis().nBasis(), bf_1);
      std::ifstream in;
      std::ofstream out;

      // Read in/w_hex.bf into bf_0 in basis format
      openInputFile("in/w_hex.bf", in);
      domain.fieldIo().readFieldsBasis(in, bf_0, domain.unitCell());
      in.close();

      // Write bf_0 to out/w_hex.bf in basis format
      openOutputFile("out/w_hex.bf", out);
      domain.fieldIo().writeFieldsBasis(out, bf_0, domain.unitCell());
      out.close();

      // Read out/w_hex.bf to bf_1 in basis format
      openInputFile("out/w_hex.bf", in);
      domain.fieldIo().readFieldsBasis(in, bf_1, domain.unitCell());
      in.close();

      // Compare initial and final fields in basis format
      BFieldComparison comparison;
      comparison.compare(bf_0, bf_1);
      if (verbose() > 0) {
         std::cout  << std::endl;
         std::cout  << Dbl(comparison.maxDiff(), 21, 13) << std::endl;
         std::cout  << Dbl(comparison.rmsDiff(), 21, 13) << std::endl;
      }
      TEST_ASSERT(comparison.maxDiff() < 1.0E-12);

   }

   void testBasisIo_bcc() 
   {
      printMethod(TEST_FUNC);

      Rp::Domain<3,CPT> domain;
      domain.setFileMaster(fileMaster_);
      readHeader("in/w_bcc.rf", domain);

      DArray< DArray<double> > bf_0;
      allocateFields(nMonomer_, domain.basis().nBasis(), bf_0);
      DArray< DArray<double> > bf_1;
      allocateFields(nMonomer_, domain.basis().nBasis(), bf_1);
      std::ifstream in;
      std::ofstream out;

      openInputFile("in/w_bcc.bf", in);
      domain.fieldIo().readFieldsBasis(in, bf_0, domain.unitCell());
      in.close();

      openOutputFile("out/w_bcc.bf", out);
      domain.fieldIo().writeFieldsBasis(out, bf_0, domain.unitCell());
      out.close();

      openInputFile("out/w_bcc.bf", in);
      domain.fieldIo().readFieldsBasis(in, bf_1, domain.unitCell());
      in.close();

      BFieldComparison comparison;
      comparison.compare(bf_0, bf_1);
      if (verbose() > 0) {
         std::cout  << std::endl;
         std::cout  << Dbl(comparison.maxDiff(), 21, 13) << std::endl;
         std::cout  << Dbl(comparison.rmsDiff(), 21, 13) << std::endl;
      }
      TEST_ASSERT(comparison.maxDiff() < 1.0E-12);

   }

   void testBasisIo_c15_1() 
   {
      printMethod(TEST_FUNC);

      // Initialize domain by reading header of r-grid file
      Rp::Domain<3,CPT> domain;
      domain.setFileMaster(fileMaster_);
      readHeader("in/c_c15_1.rf", domain);

      // Create containers and streams
      DArray< DArray<double> > bf_0;
      allocateFields(nMonomer_, domain.basis().nBasis(), bf_0);
      DArray< DArray<double> > bf_1;
      allocateFields(nMonomer_, domain.basis().nBasis(), bf_1);
      std::ifstream in;
      std::ofstream out;

      // Read to bf_0 from input file in basis format
      openInputFile("in/w_c15_1.bf", in);
      domain.fieldIo().readFieldsBasis(in, bf_0, domain.unitCell());
      in.close();

      // Write from bf_0 to output file in basis format
      openOutputFile("out/w_c15_1.bf", out);
      domain.fieldIo().writeFieldsBasis(out, bf_0, domain.unitCell());
      out.close();

      // Read to bf_1 from output file in basis format
      openInputFile("out/w_c15_1.bf", in);
      domain.fieldIo().readFieldsBasis(in, bf_1, domain.unitCell());
      in.close();

      // Compare bf_0 and bf_1
      BFieldComparison comparison;
      comparison.compare(bf_0, bf_1);
      if (verbose() > 0) {
         std::cout  << std::endl;
         std::cout  << Dbl(comparison.maxDiff(),21,13) << std::endl;
         std::cout  << Dbl(comparison.rmsDiff(),21,13) << std::endl;
      }
      TEST_ASSERT(comparison.maxDiff() < 1.0E-12);

   }

   void testBasisIo_altG() 
   {
      printMethod(TEST_FUNC);

      Rp::Domain<3,CPT> domain;
      domain.setFileMaster(fileMaster_);
      readHeader("in/w_altG.rf", domain);

      DArray< DArray<double> > bf_0;
      allocateFields(nMonomer_, domain.basis().nBasis(), bf_0);

      DArray< DArray<double> > bf_1;
      allocateFields(nMonomer_, domain.basis().nBasis(), bf_1);

      std::ifstream in;
      openInputFile("in/w_altG.bf", in);
      domain.fieldIo().readFieldsBasis(in, bf_0, domain.unitCell());
      in.close();

      std::ofstream out;
      openOutputFile("out/w_altG.bf", out);
      domain.fieldIo().writeFieldsBasis(out, bf_0, domain.unitCell());
      out.close();

      openInputFile("out/w_altG.bf", in);
      domain.fieldIo().readFieldsBasis(in, bf_1, domain.unitCell());
      in.close();

      BFieldComparison comparison;
      comparison.compare(bf_0, bf_1);
      if (verbose() > 0) {
         std::cout  << std::endl;
         std::cout  << Dbl(comparison.maxDiff(),21,13) << std::endl;
         std::cout  << Dbl(comparison.rmsDiff(),21,13) << std::endl;
      }
      TEST_ASSERT(comparison.maxDiff() < 1.0E-12);

   }

   void testBasisIo_altG_fort() 
   {
      printMethod(TEST_FUNC);

      Rp::Domain<3,CPT> domain;
      domain.setFileMaster(fileMaster_);
      readHeader("in/w_altG.rf", domain);

      DArray< DArray<double> > bf_0;
      allocateFields(nMonomer_, domain.basis().nBasis(), bf_0);

      DArray< DArray<double> > bf_1;
      allocateFields(nMonomer_, domain.basis().nBasis(), bf_1);

      std::ifstream in;
      openInputFile("in/w_altG_fort.bf", in);
      domain.fieldIo().readFieldsBasis(in, bf_0, domain.unitCell());
      in.close();

      std::ofstream out;
      openOutputFile("out/w_altG_fort.bf", out);
      domain.fieldIo().writeFieldsBasis(out, bf_0, domain.unitCell());
      out.close();

      openInputFile("out/w_altG_fort.bf", in);
      domain.fieldIo().readFieldsBasis(in, bf_1, domain.unitCell());
      in.close();

      //setVerbose(1);
      BFieldComparison comparison;
      comparison.compare(bf_0, bf_1);
      if (verbose() > 0) {
         std::cout  << std::endl;
         std::cout  << Dbl(comparison.maxDiff(), 21, 12) << std::endl;
         std::cout  << Dbl(comparison.rmsDiff(), 21, 12) << std::endl;
      }
      TEST_ASSERT(comparison.rmsDiff() < 1.0E-11);

   }

   /*
   * Read a basis file whose header gives a unit cell other than the one
   * used to construct the basis, as a file written after a flexible-cell
   * iteration does. The basis is that of a rhombic cell with gamma =
   * 1.9, the file header gives gamma = 2.1, and the minimum images of
   * some waves differ between these cells.
   */
   void testBasisIo_cellChange()
   {
      printMethod(TEST_FUNC);

      Rp::Domain<2,CPT> domain;
      domain.setFileMaster(fileMaster_);
      readHeaderString("format 1 0 dim 2 crystal_system rhombic "
                       "N_cell_param 2 cell_param 1.7 1.9 "
                       "group_name p_2 N_monomer 2 mesh 12 12", domain);
      int nBasis = domain.basis().nBasis();

      DArray< DArray<double> > bf_0;
      allocateFields(nMonomer_, nBasis, bf_0);
      for (int i = 0; i < nMonomer_; ++i) {
         for (int j = 0; j < nBasis; ++j) {
            bf_0[i][j] = 1.0/double(1 + i + j);
         }
      }

      // Write bf_0 with the header cell gamma = 2.1, read it back
      UnitCell<2> unitCell;
      FSArray<double, 6> parameters;
      parameters.append(1.7);
      parameters.append(2.1);
      unitCell.set(UnitCell<2>::Rhombic, parameters);
      std::stringstream file;
      domain.fieldIo().writeFieldsBasis(file, bf_0, unitCell);
      DArray< DArray<double> > bf_1;
      allocateFields(nMonomer_, nBasis, bf_1);
      bool thrown = false;
      try {
         domain.fieldIo().readFieldsBasis(file, bf_1, unitCell);
      } catch (Exception& e) {
         thrown = true;
      }
      TEST_ASSERT(!thrown);

      BFieldComparison comparison;
      comparison.compare(bf_0, bf_1);
      if (verbose() > 0) {
         std::cout  << std::endl;
         std::cout  << Dbl(comparison.maxDiff(),21,13) << std::endl;
      }
      TEST_ASSERT(comparison.maxDiff() < 1.0E-10);
   }

   /*
   * Write basis fields for a rectangular cell and space group "group"
   * on a mesh of n0 x n0 points, read the file on a mesh of n1 x n1
   * points, and compare the Fourier components of all waves that both
   * meshes represent without aliasing.
   */
   void checkBasisIoMesh(std::string group, int n0, int n1)
   {
      std::string header = "format 1 0 dim 2 crystal_system rectangular "
                           "N_cell_param 2 cell_param 1.7 2.1 group_name "
                           + group + " N_monomer 2 mesh ";
      Rp::Domain<2,CPT> domain0, domain1;
      domain0.setFileMaster(fileMaster_);
      domain1.setFileMaster(fileMaster_);
      readHeaderString(header + std::to_string(n0) + " "
                       + std::to_string(n0), domain0);
      readHeaderString(header + std::to_string(n1) + " "
                       + std::to_string(n1), domain1);

      DArray< DArray<double> > bf_0, bf_1;
      allocateFields(nMonomer_, domain0.basis().nBasis(), bf_0);
      allocateFields(nMonomer_, domain1.basis().nBasis(), bf_1);
      for (int i = 0; i < nMonomer_; ++i) {
         for (int j = 0; j < domain0.basis().nBasis(); ++j) {
            bf_0[i][j] = 1.0/double(1 + i + j);
         }
      }

      // Write bf_0 on mesh n0, read it into bf_1 on mesh n1
      std::stringstream file;
      domain0.fieldIo().writeFieldsBasis(file, bf_0, domain0.unitCell());
      bool thrown = false;
      try {
         domain1.fieldIo().readFieldsBasis(file, bf_1, domain1.unitCell());
      } catch (Exception& e) {
         thrown = true;
      }
      TEST_ASSERT(!thrown);

      // Compare components of waves k with |k_i| < n/2 for n = min(n0,n1)
      DArray< RFieldDft<2,CPT> > kf_0, kf_1;
      allocateFields(nMonomer_, domain0.mesh().dimensions(), kf_0);
      allocateFields(nMonomer_, domain1.mesh().dimensions(), kf_1);
      domain0.fieldIo().convertBasisToKGrid(bf_0, kf_0);
      domain1.fieldIo().convertBasisToKGrid(bf_1, kf_1);
      Mesh<2> dftMesh0(kf_0[0].dftDimensions());
      Mesh<2> dftMesh1(kf_1[0].dftDimensions());
      int n = (n0 < n1) ? n0 : n1;
      IntVec<2> k, p0, p1;
      int r0, r1;
      double maxDiff = 0.0;
      for (k[0] = 1 - n/2; k[0] < n/2; ++k[0]) {
         for (k[1] = 0; k[1] < n/2; ++k[1]) {
            p0 = k;
            p1 = k;
            if (k[0] < 0) {
               p0[0] += n0;
               p1[0] += n1;
            }
            r0 = dftMesh0.rank(p0);
            r1 = dftMesh1.rank(p1);
            for (int i = 0; i < nMonomer_; ++i) {
               for (int j = 0; j < 2; ++j) {
                  double diff = std::abs(kf_0[i][r0][j] - kf_1[i][r1][j]);
                  if (diff > maxDiff) maxDiff = diff;
               }
            }
         }
      }
      if (verbose() > 0) {
         std::cout << std::endl << group << " " << n0 << " -> " << n1
                   << " : " << Dbl(maxDiff, 21, 13) << std::endl;
      }
      TEST_ASSERT(maxDiff < 1.0E-10);
   }

   /*
   * A line of a basis file for star, with field components 0.5.
   */
   std::string starLine(Basis<2>::Star const & star)
   {
      std::ostringstream out;
      out << "0.5 0.5 " << star.waveMin[0] << " " << star.waveMin[1]
          << " " << star.size;
      return out.str();
   }

   /*
   * Write basis fields of a domain, replace some of the data lines,
   * read the result, and return true iff the read throws.
   *
   * \param lines  new data lines, indexed by line number
   */
   bool readReplacedBasisThrows(Rp::Domain<2,CPT>& domain,
                                std::map<int, std::string> const & lines)
   {
      DArray< DArray<double> > bf;
      allocateFields(nMonomer_, domain.basis().nBasis(), bf);
      for (int i = 0; i < nMonomer_; ++i) {
         for (int j = 0; j < domain.basis().nBasis(); ++j) {
            bf[i][j] = 1.0/double(1 + i + j);
         }
      }
      std::stringstream file;
      domain.fieldIo().writeFieldsBasis(file, bf, domain.unitCell());
      std::string line, text;
      int iData = -3;  // index of a data line, -3 in the header
      while (std::getline(file, line)) {
         if (line.find("N_basis") != std::string::npos) iData = -2;
         if (lines.count(iData)) line = lines.at(iData);
         text += line + "\n";
         if (iData > -3) ++iData;
      }
      std::istringstream in(text);
      try {
         domain.fieldIo().readFieldsBasis(in, bf, domain.unitCell());
      } catch (Exception& e) {
         return true;
      }
      return false;
   }

   /*
   * Read basis files written for another mesh.
   */
   void testBasisIo_mesh()
   {
      printMethod(TEST_FUNC);

      // Group without inversion, file for a coarser mesh
      checkBasisIoMesh("p_m", 12, 16);

      // File for a finer mesh, with more lines than the basis has stars
      checkBasisIoMesh("p_m", 16, 12);

      // File for a finer mesh, with stars cancelled on the coarser mesh
      checkBasisIoMesh("p_g", 16, 12);

      // A file for the same mesh must still be rejected if it lists a
      // cancelled star, or an open star that is not followed by its
      // partner (in the next line, or at the end of the file)
      Rp::Domain<2,CPT> domain;
      domain.setFileMaster(fileMaster_);
      readHeaderString("format 1 0 dim 2 crystal_system rectangular "
                       "N_cell_param 2 cell_param 1.7 2.1 group_name p_g "
                       "N_monomer 2 mesh 12 12", domain);
      Basis<2> const & basis = domain.basis();
      int nb = basis.nBasis();
      int ic = -1, io = -1, il = -1;  // cancelled, first open, last star
      for (int is = 0; is < basis.nStar(); ++is) {
         if (basis.star(is).cancel) {
            if (ic < 0) ic = is;
         } else {
            if (io < 0 && basis.star(is).invertFlag == 1) io = is;
            il = is;
         }
      }
      TEST_ASSERT(ic >= 0 && io >= 0);
      std::map<int, std::string> lines;
      lines[0] = starLine(basis.star(ic));
      TEST_ASSERT(readReplacedBasisThrows(domain, lines));
      lines.clear();
      lines[basis.star(io).basisId + 1] = starLine(basis.star(0));
      TEST_ASSERT(readReplacedBasisThrows(domain, lines));
      lines.clear();
      if (basis.star(il).invertFlag == -1) {
         lines[nb - 2] = starLine(basis.star(0));
      }
      lines[nb - 1] = starLine(basis.star(io));
      TEST_ASSERT(readReplacedBasisThrows(domain, lines));
   }

   void testRGridIo_lam() 
   {
      printMethod(TEST_FUNC);

      Rp::Domain<1,CPT> domain;
      domain.setFileMaster(fileMaster_);
      readHeader("in/w_lam.rf", domain);

      DArray< RField<1,CPT> > rf_0;
      allocateFields(nMonomer_, domain.mesh().dimensions(), rf_0);
      DArray< RField<1,CPT> > rf_1;
      allocateFields(nMonomer_, domain.mesh().dimensions(), rf_1);

      readFields("in/w_lam.rf", domain, rf_0);
      writeFields("out/w_lam.rf", domain, rf_0);
      readFields("out/w_lam.rf", domain, rf_1);

      RFieldComparison<1,CPT> comparison;
      comparison.compare(rf_0, rf_1);
      if (verbose() > 0) {
         std::cout  << "\n";
         std::cout  << Dbl(comparison.maxDiff(), 21, 13) << "\n";
         std::cout  << Dbl(comparison.rmsDiff(), 21, 13) << "\n";
      }
      TEST_ASSERT(comparison.maxDiff() < 1.0E-12);
   }

   void testRGridIo_hex() 
   {
      printMethod(TEST_FUNC);

      Rp::Domain<2,CPT> domain;
      domain.setFileMaster(fileMaster_);
      readHeader("in/w_hex.rf", domain);

      DArray< RField<2,CPT> > rf_0;
      allocateFields(nMonomer_, domain.mesh().dimensions(), rf_0);
      DArray< RField<2,CPT> > rf_1;
      allocateFields(nMonomer_, domain.mesh().dimensions(), rf_1);

      readFields("in/w_hex.rf", domain, rf_0);
      writeFields("out/w_hex.rf", domain, rf_0);
      readFields("out/w_hex.rf", domain, rf_1);

      RFieldComparison<2,CPT> comparison;
      comparison.compare(rf_0, rf_1);
      if (verbose() > 0) {
         std::cout  << "\n";
         std::cout  << Dbl(comparison.maxDiff(), 21, 13) << "\n";
         std::cout  << Dbl(comparison.rmsDiff(), 21, 13) << "\n";
      }
      TEST_ASSERT(comparison.maxDiff() < 1.0E-12);
   }

   /*
   * Conversions between r-grid and k-grid files, and rescaling of an
   * r-grid file, without a space group in the domain or in the file.
   */
   void testConvertRGridKGrid_noGroup()
   {
      printMethod(TEST_FUNC);

      // Domain of in/w_hex.rf, but without a space group
      Rp::Domain<2,CPT> domain;
      domain.setFileMaster(fileMaster_);
      readHeaderString("format 1 0 dim 2 crystal_system hexagonal "
                       "N_cell_param 1 cell_param 1.6651597308 "
                       "N_monomer 2 mesh 32 32", domain);
      TEST_ASSERT(!domain.hasGroup());
      domain.fieldIo().setNMonomer(nMonomer_);

      // Convert in/w_hex.rf (whose header gives group p_6_m_m)
      bool thrown = false;
      try {
         domain.fieldIo().convertRGridToKGrid("in/w_hex.rf",
                                              "out/w_hex_ng.kf");
         domain.fieldIo().convertKGridToRGrid("out/w_hex_ng.kf",
                                              "out/w_hex_ng.rf");
         domain.fieldIo().scaleFieldsRGrid("in/w_hex.rf",
                                           "out/w_hex_ng_x2.rf", 2.0);
      } catch (Exception& e) {
         thrown = true;
      }
      TEST_ASSERT(!thrown);

      DArray< RField<2,CPT> > rf_0, rf_1;
      allocateFields(nMonomer_, domain.mesh().dimensions(), rf_0);
      allocateFields(nMonomer_, domain.mesh().dimensions(), rf_1);
      readFields("in/w_hex.rf", domain, rf_0);
      readFields("out/w_hex_ng.rf", domain, rf_1);
      RFieldComparison<2,CPT> comparison;
      comparison.compare(rf_0, rf_1);
      TEST_ASSERT(comparison.maxDiff() < 1.0E-10);

      // Domain with group p_6_m_m, and an r-grid file without a group
      Rp::Domain<2,CPT> domainG;
      domainG.setFileMaster(fileMaster_);
      readHeader("in/w_hex.rf", domainG);
      domainG.fieldIo().setNMonomer(nMonomer_);
      rf_0[0][0] += 0.1;
      domainG.fieldIo().writeFieldsRGrid("out/w_hex_asym.rf", rf_0,
                                         domainG.unitCell(), false);

      // The output files must not declare the group
      domainG.fieldIo().convertRGridToKGrid("out/w_hex_asym.rf",
                                            "out/w_hex_asym.kf");
      domainG.fieldIo().convertKGridToRGrid("out/w_hex_asym.kf",
                                            "out/w_hex_asym_2.rf");
      bool isSymmetric;
      isSymmetric = domainG.fieldIo().readFieldsRGrid("out/w_hex_asym_2.rf",
                                               rf_1, domainG.unitCell());
      TEST_ASSERT(!isSymmetric);
   }

   /*
   * Write a single field without symmetry to a named r-grid file (as
   * WRITE_Q_SLICE does in a system without a group), and read it back.
   */
   void testRGridIo_singleField()
   {
      printMethod(TEST_FUNC);

      Rp::Domain<1,CPT> domain;
      domain.setFileMaster(fileMaster_);
      readHeader("in/w_lam.rf", domain);
      DArray< RField<1,CPT> > rf;
      allocateFields(nMonomer_, domain.mesh().dimensions(), rf);
      readFields("in/w_lam.rf", domain, rf);

      bool isSymmetric = false;
      domain.fieldIo().writeFieldRGrid("out/w_lam_0.rf", rf[0],
                                       domain.unitCell(), isSymmetric);
      RField<1,CPT> field;
      bool thrown = false;
      try {
         isSymmetric = domain.fieldIo().readFieldRGrid("out/w_lam_0.rf",
                                               field, domain.unitCell());
      } catch (Exception& e) {
         thrown = true;
      }
      TEST_ASSERT(!thrown);
      TEST_ASSERT(!isSymmetric);

      RFieldComparison<1,CPT> comparison;
      comparison.compare(rf[0], field);
      TEST_ASSERT(comparison.maxDiff() < 1.0E-10);
   }

   void testRGridIo_bcc() 
   {
      printMethod(TEST_FUNC);

      Rp::Domain<3,CPT> domain;
      domain.setFileMaster(fileMaster_);
      readHeader("in/w_bcc.rf", domain);

      DArray< RField<3,CPT> > rf_0;
      allocateFields(nMonomer_, domain.mesh().dimensions(), rf_0);
      DArray< RField<3,CPT> > rf_1;
      allocateFields(nMonomer_, domain.mesh().dimensions(), rf_1);

      readFields("in/w_bcc.rf", domain, rf_0);
      writeFields("out/w_bcc.rf", domain, rf_0);
      readFields("out/w_bcc.rf", domain, rf_1);

      RFieldComparison<3,CPT> comparison;
      comparison.compare(rf_0, rf_1);
      if (verbose() > 0) {
         std::cout  << "\n";
         std::cout  << Dbl(comparison.maxDiff(),21,13) << "\n";
         std::cout  << Dbl(comparison.rmsDiff(),21,13) << "\n";
      }
      TEST_ASSERT(comparison.maxDiff() < 1.0E-12);
   }

   void testConvertBasisKGridBasis_bcc() 
   {
      printMethod(TEST_FUNC);

      Rp::Domain<3,CPT> domain;
      domain.setFileMaster(fileMaster_);
      readHeader("in/w_bcc.rf", domain);

      DArray< DArray<double> > bf_0;
      allocateFields(nMonomer_, domain.basis().nBasis(), bf_0);
      DArray< DArray<double> > bf_1;
      allocateFields(nMonomer_, domain.basis().nBasis(), bf_1);
      DArray< RFieldDft<3,CPT> > kf_0;
      allocateFields(nMonomer_, domain.mesh().dimensions(), kf_0);

      readFields("in/w_bcc.bf", domain, bf_0);
      domain.fieldIo().convertBasisToKGrid(bf_0, kf_0);
      domain.fieldIo().convertKGridToBasis(kf_0, bf_1);

      BFieldComparison comparison;
      comparison.compare(bf_0, bf_1);
      TEST_ASSERT(comparison.maxDiff() < 1.0E-12);

      if (verbose() > 0) {
         std::cout  << "\n";
         std::cout  << Dbl(comparison.maxDiff(),21,13) << "\n";
         std::cout  << Dbl(comparison.rmsDiff(),21,13) << "\n";
      }

   }

   void testConvertBasisRGridBasis_bcc() 
   {
      printMethod(TEST_FUNC);

      Rp::Domain<3,CPT> domain;
      domain.setFileMaster(fileMaster_);
      readHeader("in/w_bcc.rf", domain);

      DArray< DArray<double> > bf_0;
      allocateFields(nMonomer_, domain.basis().nBasis(), bf_0);
      DArray< DArray<double> > bf_1;
      allocateFields(nMonomer_, domain.basis().nBasis(), bf_1);
      DArray< RField<3,CPT> > rf_0;
      allocateFields(nMonomer_, domain.mesh().dimensions(), rf_0);

      readFields("in/w_bcc.bf", domain, bf_0);
      domain.fieldIo().convertBasisToRGrid(bf_0, rf_0);
      domain.fieldIo().convertRGridToBasis(rf_0, bf_1);

      BFieldComparison comparison;
      comparison.compare(bf_0, bf_1);
      TEST_ASSERT(comparison.maxDiff() < 1.0E-12);

      if (verbose() > 0) {
         std::cout  << "\n";
         std::cout  << Dbl(comparison.maxDiff(),21,13) << "\n";
         std::cout  << Dbl(comparison.rmsDiff(),21,13) << "\n";
      }
   }

   void testConvertBasisKGridBasis_altG() 
   {
      printMethod(TEST_FUNC);
      nMonomer_ = 3;

      Rp::Domain<3,CPT> domain;
      domain.setFileMaster(fileMaster_);
      readHeader("in/w_altG.rf", domain);

      std::ofstream  out;
      openOutputFile("out/altG.stars", out);
      domain.basis().outputStars(out);
      out.close();

      #if 1
      DArray< DArray<double> > bf_0;
      allocateFields(nMonomer_, domain.basis().nBasis(), bf_0);
      DArray< DArray<double> > bf_1;
      allocateFields(nMonomer_, domain.basis().nBasis(), bf_1);
      DArray< RFieldDft<3,CPT> > kf_0;
      allocateFields(nMonomer_, domain.mesh().dimensions(), kf_0);

      readFields("in/w_altG.bf", domain, bf_0);
      domain.fieldIo().convertBasisToKGrid(bf_0, kf_0);
      domain.fieldIo().convertKGridToBasis(kf_0, bf_1);

      BFieldComparison comparison;
      comparison.compare(bf_0, bf_1);

      if (verbose() > 0) {
         std::cout  << "\n";
         std::cout  << Dbl(comparison.maxDiff(),21,13) << "\n";
         std::cout  << Dbl(comparison.rmsDiff(),21,13) << "\n";
      }
      #endif

      TEST_ASSERT(comparison.maxDiff() < 1.0E-10);
   }

   void testConvertBasisKGridBasis_c15_1() 
   {
      printMethod(TEST_FUNC);
      nMonomer_ = 2;

      Rp::Domain<3,CPT> domain;
      domain.setFileMaster(fileMaster_);
      readHeader("in/c_c15_1.rf", domain);

      #if 0
      std::ofstream  out;
      openOutputFile("out/waves_c15_1", out);
      domain.basis().outputWaves(out);
      out.close();
      #endif

      DArray< DArray<double> > bf_0;
      allocateFields(nMonomer_, domain.basis().nBasis(), bf_0);
      DArray< DArray<double> > bf_1;
      allocateFields(nMonomer_, domain.basis().nBasis(), bf_1);
      DArray< RFieldDft<3,CPT> > kf_0;
      allocateFields(nMonomer_, domain.mesh().dimensions(), kf_0);

      readFields("in/w_c15_1.bf", domain, bf_0);
      domain.fieldIo().convertBasisToKGrid(bf_0, kf_0);
      domain.fieldIo().convertKGridToBasis(kf_0, bf_1);

      BFieldComparison comparison;
      comparison.compare(bf_0, bf_1);

      // setVerbose(1);
      if (verbose() > 0) {
         std::cout  << "\n";
         std::cout  << Dbl(comparison.maxDiff(),21,13) << "\n";
         std::cout  << Dbl(comparison.rmsDiff(),21,13) << "\n";
      }
      TEST_ASSERT(comparison.rmsDiff() < 1.0E-10);
   }

   void testEstimateWBasis() 
   {
      printMethod(TEST_FUNC);
      nMonomer_ = 2;

      Rp::Domain<3,CPT> domain;
      domain.setFileMaster(fileMaster_);
      domain.fieldIo().setNMonomer(nMonomer_);
      readHeader("in/c_c15_1.rf", domain);

      Interaction v;
      v.setNMonomer(nMonomer_);
      std::ifstream in;
      openInputFile("in/Interaction", in);
      v.readParam(in);
      in.close();

      domain.fieldIo().estimateWBasis("in/c_c15_1.bf","out/w_c15_est.bf",
                                      v.chi());

      DArray< DArray<double> > bf_0, bf_1;
      allocateFields(nMonomer_, domain.basis().nBasis(), bf_0);
      allocateFields(nMonomer_, domain.basis().nBasis(), bf_1);
      readFields("out/w_c15_est.bf", domain, bf_0);
      readFields("in/w_c15_est.bf", domain, bf_1);

      BFieldComparison comparison;
      comparison.compare(bf_0, bf_1);

      if (verbose() > 0) {
         std::cout  << "\n";
         std::cout  << Dbl(comparison.maxDiff(),21,13) << "\n";
         std::cout  << Dbl(comparison.rmsDiff(),21,13) << "\n";
      }
      TEST_ASSERT(comparison.rmsDiff() < 1.0E-8);
   }

   void testKGridIo_lam() 
   {
      printMethod(TEST_FUNC);

      Rp::Domain<1,CPT> domain;
      domain.setFileMaster(fileMaster_);
      readHeader("in/w_lam.rf", domain);

      DArray< DArray<double> > bf_0;
      allocateFields(nMonomer_, domain.basis().nBasis(), bf_0);

      DArray< RFieldDft<1,CPT> > kf_0, kf_1;
      allocateFields(nMonomer_, domain.mesh().dimensions(), kf_0);
      allocateFields(nMonomer_, domain.mesh().dimensions(), kf_1);

      // Read basis form to bf_0, convert to k-grid form kf_0
      readFields("in/w_lam.bf", domain, bf_0);
      domain.fieldIo().convertBasisToKGrid(bf_0, kf_0);

      // Write kf_0 and read back into kf_1
      writeFields("out/w_lam.kf", domain, kf_0);
      readFields("out/w_lam.kf", domain, kf_1);

      // Compare kf_0 and kf_1 in k-grid format
      RFieldDftComparison<1,CPT> comparison;
      comparison.compare(kf_0, kf_1);
      if (verbose() > 0) {
         std::cout  << std::endl;
         std::cout  << Dbl(comparison.maxDiff(), 21, 13) << "\n";
         std::cout  << Dbl(comparison.rmsDiff(), 21, 13) << "\n";
      }
      TEST_ASSERT(comparison.maxDiff() < 1.0E-12);
   }

   /*
   * Read field files in which the last number is not followed by a
   * newline.
   */
   void testReadWithoutFinalNewline_lam()
   {
      printMethod(TEST_FUNC);

      Rp::Domain<1,CPT> domain;
      domain.setFileMaster(fileMaster_);
      readHeader("in/w_lam.rf", domain);
      Rp::FieldIo<1,CPT> const & fieldIo = domain.fieldIo();

      DArray< DArray<double> > bf_0, bf_1;
      allocateFields(nMonomer_, domain.basis().nBasis(), bf_0);
      allocateFields(nMonomer_, domain.basis().nBasis(), bf_1);
      DArray< RField<1,CPT> > rf_0, rf_1;
      allocateFields(nMonomer_, domain.mesh().dimensions(), rf_0);
      allocateFields(nMonomer_, domain.mesh().dimensions(), rf_1);
      DArray< RFieldDft<1,CPT> > kf_0, kf_1;
      allocateFields(nMonomer_, domain.mesh().dimensions(), kf_0);
      allocateFields(nMonomer_, domain.mesh().dimensions(), kf_1);
      readFields("in/w_lam.bf", domain, bf_0);
      fieldIo.convertBasisToRGrid(bf_0, rf_0);
      fieldIo.convertBasisToKGrid(bf_0, kf_0);

      // Write the fields in each format, remove the final newline
      std::ostringstream outB, outR, outK;
      fieldIo.writeFieldsBasis(outB, bf_0, domain.unitCell());
      fieldIo.writeFieldsRGrid(outR, rf_0, domain.unitCell());
      fieldIo.writeFieldsKGrid(outK, kf_0, domain.unitCell());
      std::string textB = outB.str();
      std::string textR = outR.str();
      std::string textK = outK.str();
      textB.erase(textB.find_last_not_of(" \n") + 1);
      textR.erase(textR.find_last_not_of(" \n") + 1);
      textK.erase(textK.find_last_not_of(" \n") + 1);
      std::istringstream inB(textB), inR(textR), inK(textK);

      bool thrown = false;
      try {
         fieldIo.readFieldsBasis(inB, bf_1, domain.unitCell());
         fieldIo.readFieldsRGrid(inR, rf_1, domain.unitCell());
         fieldIo.readFieldsKGrid(inK, kf_1, domain.unitCell());
      } catch (Exception& e) {
         thrown = true;
      }
      TEST_ASSERT(!thrown);

      BFieldComparison comparisonB;
      comparisonB.compare(bf_0, bf_1);
      TEST_ASSERT(comparisonB.maxDiff() < 1.0E-10);
      RFieldComparison<1,CPT> comparisonR;
      comparisonR.compare(rf_0, rf_1);
      TEST_ASSERT(comparisonR.maxDiff() < 1.0E-10);
      RFieldDftComparison<1,CPT> comparisonK;
      comparisonK.compare(kf_0, kf_1);
      TEST_ASSERT(comparisonK.maxDiff() < 1.0E-10);
   }

   /*
   * Write k-grid fields with negative components that have three-digit
   * exponents. Every data line must have 1 + 2*nMonomer columns.
   */
   void testKGridIo_columns()
   {
      printMethod(TEST_FUNC);

      Rp::Domain<1,CPT> domain;
      domain.setFileMaster(fileMaster_);
      readHeader("in/w_lam.rf", domain);

      DArray< RFieldDft<1,CPT> > kf;
      allocateFields(nMonomer_, domain.mesh().dimensions(), kf);
      for (int i = 0; i < nMonomer_; ++i) {
         for (int j = 0; j < kf[i].capacity(); ++j) {
            kf[i][j][0] = -1.5E-100;
            kf[i][j][1] = -2.5E-100;
         }
      }
      std::stringstream file;
      domain.fieldIo().writeFieldsKGrid(file, kf, domain.unitCell());

      // Skip the header, which ends with the mesh dimensions
      std::string line;
      while (std::getline(file, line)) {
         if (line.find("mesh") != std::string::npos) break;
      }
      std::getline(file, line);

      // Count the columns of every data line
      int nLine = 0;
      while (std::getline(file, line)) {
         std::istringstream lineStream(line);
         std::string word;
         int nColumn = 0;
         while (lineStream >> word) ++nColumn;
         TEST_ASSERT(nColumn == 1 + 2*nMonomer_);
         ++nLine;
      }
      TEST_ASSERT(nLine == kf[0].capacity());
   }

   void testKGridIo_bcc() 
   {
      printMethod(TEST_FUNC);

      Rp::Domain<3,CPT> domain;
      domain.setFileMaster(fileMaster_);
      readHeader("in/w_bcc.rf", domain);

      DArray< DArray<double> > bf_0;
      allocateFields(nMonomer_, domain.basis().nBasis(), bf_0);
      DArray< RFieldDft<3,CPT> > kf_0;
      allocateFields(nMonomer_, domain.mesh().dimensions(), kf_0);
      DArray< RFieldDft<3,CPT> > kf_1;
      allocateFields(nMonomer_, domain.mesh().dimensions(), kf_1);

      readFields("in/w_bcc.bf", domain, bf_0);
      domain.fieldIo().convertBasisToKGrid(bf_0, kf_0);

      writeFields("out/w_bcc.kf", domain, kf_0);
      readFields("out/w_bcc.kf", domain, kf_1);

      RFieldDftComparison<3,CPT> comparison;
      comparison.compare(kf_0, kf_1);
      if (verbose() > 0) {
         std::cout  << "\n";
         std::cout  << Dbl(comparison.maxDiff(),21,13) << "\n";
         std::cout  << Dbl(comparison.rmsDiff(),21,13) << "\n";
      }
      TEST_ASSERT(comparison.maxDiff() < 1.0E-11);

   }

   void testKGridIo_altG() 
   {
      printMethod(TEST_FUNC);

      Rp::Domain<3,CPT> domain;
      domain.setFileMaster(fileMaster_);
      readHeader("in/w_altG.rf", domain);

      DArray< DArray<double> > bf_0;
      allocateFields(nMonomer_, domain.basis().nBasis(), bf_0);
      DArray< RFieldDft<3,CPT> > kf_0;
      allocateFields(nMonomer_, domain.mesh().dimensions(), kf_0);
      DArray< RFieldDft<3,CPT> > kf_1;
      allocateFields(nMonomer_, domain.mesh().dimensions(), kf_1);

      readFields("in/w_altG.bf", domain, bf_0);
      domain.fieldIo().convertBasisToKGrid(bf_0, kf_0);

      writeFields("out/w_altG.kf", domain, kf_0);
      readFields("out/w_altG.kf", domain, kf_1);

      RFieldDftComparison<3,CPT> comparison;
      comparison.compare(kf_0, kf_1);
      if (verbose() > 0) {
         std::cout  << "\n";
         std::cout  << Dbl(comparison.maxDiff(),21,13) << "\n";
         std::cout  << Dbl(comparison.rmsDiff(),21,13) << "\n";
      }
      TEST_ASSERT(comparison.maxDiff() < 1.0E-11);
   }

   void testConvertBasisKGridRGridKGrid_bcc() 
   {
      printMethod(TEST_FUNC);

      Rp::Domain<3,CPT> domain;
      domain.setFileMaster(fileMaster_);
      readHeader("in/w_bcc.rf", domain);

      DArray< DArray<double> > bf_0;
      allocateFields(nMonomer_, domain.basis().nBasis(), bf_0);

      DArray< DArray<double> > bf_1;
      allocateFields(nMonomer_, domain.basis().nBasis(), bf_1);

      DArray< RFieldDft<3,CPT> > kf_0;
      allocateFields(nMonomer_, domain.mesh().dimensions(), kf_0);

      DArray< RFieldDft<3,CPT> > kf_1;
      allocateFields(nMonomer_, domain.mesh().dimensions(), kf_1);

      DArray< RFieldDft<3,CPT> > kf_2;
      allocateFields(nMonomer_, domain.mesh().dimensions(), kf_2);

      DArray< RField<3,CPT> > rf_0;
      allocateFields(nMonomer_, domain.mesh().dimensions(), rf_0);

      readFields("in/w_bcc.bf", domain, bf_0);
      domain.fieldIo().convertBasisToKGrid(bf_0, kf_0);

      kf_2 = kf_0;
      domain.fieldIo().convertKGridToRGrid(kf_0, rf_0);

      #if 0
      // Demonstrate that input kf_0 is modified by above (it is)
      RFieldDftComparison<3,CPT> check;
      check.compare(kf_2, kf_0);
      std::cout  << std::endl;
      std::cout  << Dbl(check.maxDiff(), 21, 13) << "\n";
      std::cout  << Dbl(check.rmsDiff(), 21, 13) << "\n";
      #endif

      domain.fieldIo().convertRGridToKGrid(rf_0, kf_1);

      RFieldDftComparison<3,CPT> comparison;
      comparison.compare(kf_2, kf_1);
      TEST_ASSERT(comparison.maxDiff() < 1.0E-10);
      if (verbose() > 0) {
        std::cout  << std::endl;
        std::cout  << Dbl(comparison.maxDiff(), 21, 13) << "\n";
        std::cout  << Dbl(comparison.rmsDiff(), 21, 13) << "\n";
      }
   }

   void testConvertBasisKGridRGridKGrid_c15_1() 
   {
      printMethod(TEST_FUNC);

      // Read header
      Rp::Domain<3,CPT> domain;
      domain.setFileMaster(fileMaster_);
      readHeader("in/c_c15_1.rf", domain);

      // Allocate required fields
      DArray< DArray<double> > bf_0;
      allocateFields(nMonomer_, domain.basis().nBasis(), bf_0);

      DArray< DArray<double> > bf_1;
      allocateFields(nMonomer_, domain.basis().nBasis(), bf_1);

      DArray< RFieldDft<3,CPT> > kf_0;
      allocateFields(nMonomer_, domain.mesh().dimensions(), kf_0);

      DArray< RFieldDft<3,CPT> > kf_1;
      allocateFields(nMonomer_, domain.mesh().dimensions(), kf_1);

      DArray< RFieldDft<3,CPT> > kf_2;
      allocateFields(nMonomer_, domain.mesh().dimensions(), kf_2);

      DArray< RField<3,CPT> > rf_0;
      allocateFields(nMonomer_, domain.mesh().dimensions(), rf_0);

      // Read fields in basis format 
      readFields("in/w_c15_1.bf", domain, bf_0);

      // Convert basis -> kgrid -> rgrid -> kgrid -> basis
      domain.fieldIo().convertBasisToKGrid(bf_0, kf_0);

      kf_2 = kf_0;
      domain.fieldIo().convertKGridToRGrid(kf_0, rf_0);

      #if 0
      // Demonstrate that input kf_0 is modified by above (it is)
      RFieldDftComparison<3,CPT> check;
      check.compare(kf_2, kf_0);
      std::cout  << std::endl;
      std::cout  << Dbl(check.maxDiff(), 21, 13) << "\n";
      std::cout  << Dbl(check.rmsDiff(), 21, 13) << "\n";
      #endif

      domain.fieldIo().convertRGridToKGrid(rf_0, kf_1);

      RFieldDftComparison<3,CPT> Bcomparison;
      Bcomparison.compare(kf_2, kf_1);
      TEST_ASSERT(Bcomparison.maxDiff() < 1.0E-10);
      if (verbose() > 0) {
        std::cout  << std::endl;
        std::cout  << Dbl(Bcomparison.maxDiff(), 21, 13) << "\n";
        std::cout  << Dbl(Bcomparison.rmsDiff(), 21, 13) << "\n";
      }

      domain.fieldIo().convertKGridToBasis(kf_1, bf_1);
      BFieldComparison comparison;
      comparison.compare(bf_0, bf_1);
      TEST_ASSERT(comparison.maxDiff() < 1.0E-10);

   }

   void testReplicate_bcc() 
   {
      printMethod(TEST_FUNC);

      Rp::Domain<3,CPT> domain;
      domain.setFileMaster(fileMaster_);
      readHeader("in/w_bcc.rf", domain);
      IntVec<3> dimensions = domain.mesh().dimensions();

      DArray< RField<3,CPT> > rf_0;
      allocateFields(nMonomer_, dimensions, rf_0);
      readFields("in/w_bcc.rf", domain, rf_0);

      IntVec<3> replicas;
      replicas[0] = 2;
      replicas[1] = 2;
      replicas[2] = 2;

      std::ofstream  out;
      openOutputFile("out/w_bcc_replica.rf", out);
      domain.fieldIo().replicateUnitCell(out, rf_0, 
                                         domain.unitCell(),
                                         replicas);

      // Read replicated field header 
      Rp::Domain<3,CPT> domain_rep;
      domain_rep.setFileMaster(fileMaster_);
      readHeader("out/w_bcc_replica.rf", domain_rep);
      IntVec<3> dimensions_rep = domain_rep.mesh().dimensions();
      for (int i=0; i < 3; ++i) {
         //TEST_ASSERT(domain_rep.mesh().dimension(i) 
         //            == 2*domain.mesh().dimension(i));
         TEST_ASSERT(dimensions_rep[i] == 2*dimensions[i]);
         TEST_ASSERT(dimensions[i] == 32);
      }

      // Read replicated field
      DArray< RField<3,CPT> > rf_1;
      allocateFields(nMonomer_, dimensions_rep, rf_1);
      readFields("out/w_bcc_replica.rf", domain_rep, rf_1);

      // Compare fields
      Mesh<3> mesh(dimensions);
      Mesh<3> mesh_rep(dimensions_rep);
      IntVec<3> p;
      IntVec<3> p_rep;
      int rank, rank_rep;
      for (p[0] = 0; p[0] < dimensions[0]; ++p[0]) {
         p_rep[0] = p[0] + dimensions[0];
         for (p[1] = 0; p[1] < dimensions[1]; ++p[1]) {
            p_rep[1] = p[1] + dimensions[1];
            for (p[2] = 0; p[2] < dimensions[2]; ++p[2]) {
               p_rep[2] = p[2] + dimensions[2];
               rank = mesh.rank(p);
               rank_rep = mesh_rep.rank(p_rep);
               //std::cout << Int(p[0]) << Int(p[1]) << Int(p[2])
               //          << Dbl(rf_0[0][rank])
               //          << Dbl(rf_1[0][rank_rep]) 
               //          << std::endl;
               for (int i=0; i < nMonomer_; ++i) {
                  TEST_ASSERT(eq( rf_0[i][rank], rf_1[i][rank_rep]));
               }
            }
         }
      }

   }

   void testExpand_lam_13() 
   {
      printMethod(TEST_FUNC);

      Rp::Domain<1,CPT> domain;
      domain.setFileMaster(fileMaster_);
      readHeader("in/w_lam.rf", domain);
      IntVec<1> dimensions = domain.mesh().dimensions();

      DArray< RField<1,CPT> > rf_0;
      allocateFields(nMonomer_, dimensions, rf_0);
      readFields("in/w_lam.rf", domain, rf_0);

      int d = 3;
      DArray<int> newGridDimensions;
      newGridDimensions.allocate(2);
      newGridDimensions[0] = dimensions[0];
      newGridDimensions[1] = dimensions[0];

      std::ofstream  out;
      openOutputFile("out/w_lam_exp.rf", out);
      domain.fieldIo().expandRGridDimension(out, rf_0, 
                                  domain.unitCell(), d,
                                  newGridDimensions);

      // Read header for expanded field (d=3)
      Rp::Domain<3,CPT> domain_exp;
      domain_exp.setFileMaster(fileMaster_);
      readHeader("out/w_lam_exp.rf", domain_exp);
      IntVec<3> dimensions_exp = domain_exp.mesh().dimensions();
      dimensions_exp[0] = dimensions[0];
      dimensions_exp[1] = newGridDimensions[0];
      dimensions_exp[2] = newGridDimensions[1];

      // Read data for expanded field (d=3)
      DArray< RField<3,CPT> > rf_1;
      allocateFields(nMonomer_, dimensions_exp, rf_1);
      readFields("out/w_lam_exp.rf", domain_exp, rf_1);

      // Compare fields
      Mesh<3> mesh_exp(dimensions_exp);
      IntVec<3> p;
      int rank, rank_exp;
      double value;
      for (int i=0; i < nMonomer_; ++i) {
         for (p[0] = 0; p[0] < dimensions_exp[0]; ++p[0]) {
            rank = p[0];
            value = rf_0[i][rank];
            for (p[1] = 0; p[1] < dimensions_exp[1]; ++p[1]) {
               for (p[2] = 0; p[2] < dimensions_exp[2]; ++p[2]) {
                  rank_exp = mesh_exp.rank(p);
                  TEST_ASSERT(eq(value, rf_1[i][rank_exp]));
               }
            }
         }
      }

   }

   void testExpand_hex_23() 
   {
      printMethod(TEST_FUNC);

      // Create domain and initialize by reading header
      Rp::Domain<2,CPT> domain;
      domain.setFileMaster(fileMaster_);
      readHeader("in/w_hex.rf", domain);
      IntVec<2> dimensions = domain.mesh().dimensions();

      // Read fields into container rf_0
      DArray< RField<2,CPT> > rf_0;
      allocateFields(nMonomer_, dimensions, rf_0);
      readFields("in/w_hex.rf", domain, rf_0);

      // Setup extra dimensions
      int d = 3;
      DArray<int> newGridDimensions;
      newGridDimensions.allocate(1);
      newGridDimensions[0] = dimensions[0];

      // Create and write expanded fields
      std::ofstream  out;
      openOutputFile("out/w_hex_exp.rf", out);
      domain.fieldIo().expandRGridDimension(out, rf_0, 
                                 domain.unitCell(), d,
                                 newGridDimensions);

      // Read header for expanded field (d=3)
      Rp::Domain<3,CPT> domain_exp;
      domain_exp.setFileMaster(fileMaster_);
      readHeader("out/w_hex_exp.rf", domain_exp);
      IntVec<3> dimensions_exp = domain_exp.mesh().dimensions();
      TEST_ASSERT(dimensions_exp[0] == dimensions[0]);
      TEST_ASSERT(dimensions_exp[1] == dimensions[1]);
      TEST_ASSERT(dimensions_exp[2] == newGridDimensions[0]);

      // Read data for expanded field (d=3)
      DArray< RField<3,CPT> > rf_1;
      allocateFields(nMonomer_, dimensions_exp, rf_1);
      readFields("out/w_hex_exp.rf", domain_exp, rf_1);

      // Compare fields
      Mesh<2> mesh(dimensions);
      Mesh<3> mesh_exp(dimensions_exp);
      IntVec<2> p;
      IntVec<3> p_exp;
      int rank, rank_exp;
      double value;
      for (int i=0; i < nMonomer_; ++i) {
         for (p_exp[0] = 0; p_exp[0] < dimensions_exp[0]; ++p_exp[0]) {
            p[0] = p_exp[0];
            for (p_exp[1] = 0; p_exp[1] < dimensions_exp[1]; ++p_exp[1]) {
               p[1] = p_exp[1];
               rank = mesh.rank(p);
               value = rf_0[i][rank];
               for (p_exp[2] = 0; p_exp[2] < dimensions_exp[2]; ++p_exp[2]){
                  rank_exp = mesh_exp.rank(p_exp);
                  TEST_ASSERT(eq(value, rf_1[i][rank_exp]));
               }
            }
         }
      }

   }

};

TEST_BEGIN(FieldIoTest)
TEST_ADD(FieldIoTest, testReadHeader)
TEST_ADD(FieldIoTest, testBasisIo_lam)
TEST_ADD(FieldIoTest, testBasisIo_hex)
TEST_ADD(FieldIoTest, testBasisIo_hex)
TEST_ADD(FieldIoTest, testBasisIo_bcc)
TEST_ADD(FieldIoTest, testBasisIo_c15_1)
TEST_ADD(FieldIoTest, testBasisIo_altG)
TEST_ADD(FieldIoTest, testBasisIo_altG_fort)
TEST_ADD(FieldIoTest, testBasisIo_cellChange)
TEST_ADD(FieldIoTest, testBasisIo_mesh)
TEST_ADD(FieldIoTest, testRGridIo_lam)
TEST_ADD(FieldIoTest, testRGridIo_hex)
TEST_ADD(FieldIoTest, testRGridIo_bcc)
TEST_ADD(FieldIoTest, testConvertRGridKGrid_noGroup)
TEST_ADD(FieldIoTest, testRGridIo_singleField)
TEST_ADD(FieldIoTest, testConvertBasisKGridBasis_bcc)
TEST_ADD(FieldIoTest, testConvertBasisRGridBasis_bcc)
TEST_ADD(FieldIoTest, testConvertBasisKGridBasis_altG)
TEST_ADD(FieldIoTest, testConvertBasisKGridBasis_c15_1)
TEST_ADD(FieldIoTest, testEstimateWBasis)
TEST_ADD(FieldIoTest, testKGridIo_bcc)
TEST_ADD(FieldIoTest, testKGridIo_altG)
TEST_ADD(FieldIoTest, testKGridIo_lam)
TEST_ADD(FieldIoTest, testKGridIo_columns)
TEST_ADD(FieldIoTest, testReadWithoutFinalNewline_lam)
TEST_ADD(FieldIoTest, testConvertBasisKGridRGridKGrid_bcc)
TEST_ADD(FieldIoTest, testConvertBasisKGridRGridKGrid_c15_1)
TEST_ADD(FieldIoTest, testReplicate_bcc)
TEST_ADD(FieldIoTest, testExpand_lam_13)
TEST_ADD(FieldIoTest, testExpand_hex_23)
TEST_END(FieldIoTest)

#endif

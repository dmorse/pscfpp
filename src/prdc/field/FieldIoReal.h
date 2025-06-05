#ifndef RPC_FIELD_IO_REAL_H
#define RPC_FIELD_IO_REAL_H

/*
* PSCF - Polymer Self-Consistent Field Theory
*
* Copyright 2015 - 2025, The Regents of the University of Minnesota
* Distributed under the terms of the GNU General Public License.
*/

#include <prdc/crystal/UnitCell.h>     // nested LatticeSystem enum
#include <pscf/math/IntVec.h>          // template with default
#include <util/containers/DArray.h>    // member

// Forward declarations for classes used in interfaces
namespace Util {
   class FileMaster;
}
namespace Pscf {
   template <int D> class Mesh;
   namespace Prdc {
      template <int D> class Basis;
      template <int D> class SpaceGroup;
   }
}

namespace Pscf {
namespace Prdc {

   using namespace Util;
   using namespace Pscf;
   using namespace Pscf::Prdc;

   /**
   * File input/output, format conversions and other utilities for fields.
   *
   * This class provides functions to read and write arrays that contain
   * fields in any of three representations (symmetry-adapted basis,
   * r-space grid, or Fourier k-space grid), functions to convert among 
   * these representations, and other utilities for manipulating fields. 
   * The member functions that implement field IO operations define the 
   * file formats for these field representations.
   *
   * <b>Template parameters:</b>
   *
   *    - D    : dimension of space, i.e., 1, 2, or 3
   *    - RFT  : real field r-grid type, e.g., RField<D> 
   *    - KFT  : real field k-grid type, e.g., RFieldDft<D> 
   *    - FFT  : fast Fourier transform type, e.g., FFT<D> 
   *
   * The FieldIoReal template is a base class for two class templates
   * named FieldIo<int D> defined in namespaces Pscf::Rpc and Pscf::Rpg. 
   * The Pscf::Rpc::FieldIo<int D> template is derived from a partial 
   * specialization of FieldIoReal with parameters RFT = Cpu::RField<D>, 
   * KFT = Cpu::RFieldDft<D>, and FFT = Cpu::FFT<D> that are all defined 
   * in the Pscf::Cpu namespace, and that all use standard CPU hardware. 
   * The analogous template Rpg::Field<D> in the Pscf::Rpg namespace is 
   * derived from a partial specialization of FieldIoReal in which these
   * three parameters are class templates with the same names that are 
   * defined in the Pscf::Cuda namespace, and that all use a GPU. 
   *
   * <b>Basis construction as side effect of reading field files:</b> 
   * Each member functions that read fields from a file may construct a
   * symmetry adapted basis in an associated Basis<D> object as a side 
   * effect of reading a field header. All of these functions call member 
   * function readFieldHeader to read the field file header. If an 
   * associated space group is known but the symmetry adapted basis has 
   * not been initialized upon entry, then readFieldHeader will initialize 
   * the associated Basis<D> object using the unit cell parameters found
   * in the field file header. 
   *
   * <b> Pure virtual member functions </b>: This class template defines 
   * several pure virtual functions for which different implementations 
   * are required for Cpu and Cuda code. Cpu and Cuda implementations 
   * of such functions usually differ because the Cuda versions must 
   * explicitly transfer data between Cpu and Gpu memory spaces. 
   *
   * \ingroup Prdc_Field_Module
   */
   template <int D, class RFT, class KFT, class FFT>
   class FieldIoReal
   {

   public:

      /// \name Construction, Initialization and Destruction
      ///@{

      /**
      * Constructor.
      */
      FieldIoReal();

      /**
      * Destructor.
      */
      virtual ~FieldIoReal();

      /**
      * Create association with other objects in the parent Domain.
      *
      * \param mesh  associated spatial discretization Mesh<D>
      * \param fft   associated FFT object for fast transforms
      * \param lattice  lattice system type (enumeration value)
      * \param hasGroup true if a space group is declared
      * \param groupName space group name string
      * \param group  associated SpaceGroup object
      * \param basis  associated Basis object
      */
      void associate(Mesh<D> const & mesh,
                     FFT const & fft,
                     typename UnitCell<D>::LatticeSystem const & lattice,
                     bool const & hasGroup,
                     std::string const & groupName,
                     SpaceGroup<D> const & group,
                     Basis<D> & basis);

      /**
      * Create an association with a FileMaster.
      *
      * \param fileMaster  associated FileMaster (for file paths)
      */
      void setFileMaster(FileMaster const & fileMaster);

      /**
      * Set the number of monomers types.
      * 
      * This is used to allocate arrays of fields used to as workspace
      * by some functions that read field files.
      *
      * \param nMonomer  number of monomer types
      */
      void setNMonomer(int nMonomer);

      ///@}
      /// \name Field File IO - Symmetry Adapted Basis Format
      ///@{

      /**
      * Read an array of fields in basis format from an input stream.
      *
      * This function reads fields in a symmetry adapted basis format from
      * input stream in. The header of the field file must declare a group
      * name, and this name must agree with that declared in the parameter
      * file. 
      *
      * The capacity of DArray fields is equal to nMonomer, and element
      * fields[i] is a DArray<double> containing components of the field
      * associated with monomer type i, defined as coefficients of 
      * symmetry-adapted basis functions.
      *
      * \param in  input stream (e.g., input file stream)
      * \param fields  array of fields (symmetry adapted basis components)
      * \param unitCell  associated unit cell object
      */
      void
      readFieldsBasis(std::istream& in, 
		      DArray< DArray<double> >& fields,
                      UnitCell<D> & unitCell) const;

      /**
      * Read an array of fields in basis format from a named file.
      *
      * This function opens an input file with the specified filename, 
      * reads an array of fields in symmetry-adapted form from that file, 
      * and then closes the file. 
      *
      * \param filename  name of input file
      * \param fields  array of fields (symmetry adapted basis components)
      * \param unitCell  associated crystallographic unit cell
      */
      void readFieldsBasis(std::string filename,
                           DArray< DArray<double> >& fields,
                           UnitCell<D> & unitCell) const;

      /**
      * Read a single field in basis format from an input stream.
      *
      * This function reads a single field in symmetry adapted basis 
      * format from an input stream.
      *
      * \param in  input stream (i.e., input file stream)
      * \param field  array to store the field (basis format)
      * \param unitCell  associated crystallographic unit cell
      */
      void readFieldBasis(std::istream& in, 
		          DArray<double>& field,
                          UnitCell<D> & unitCell) const;

      /**
      * Read a single field in basis format from a named file.
      *
      * This function opens an input file with the specified filename,
      * reads a single field in symmetry adapted basis format from that 
      * file, and and then closes the file. 
      *
      * \param filename  name of input file
      * \param field  array to store the field (basis format)
      * \param unitCell  associated crystallographic unit cell
      */
      void readFieldBasis(std::string filename, 
		          DArray<double>& field,
                          UnitCell<D> & unitCell) const;

      /**
      * Write an array of fields in basis format to an output stream.
      *
      * This function writes field components in a symmetry adapted basis
      * to an output stream.
      *
      * \param out  output stream (i.e., output file)
      * \param fields  array of fields (symmetry adapted basis components)
      * \param unitCell  associated crystallographic unit cell
      */
      void writeFieldsBasis(std::ostream& out,
                            DArray< DArray<double> > const & fields,
                            UnitCell<D> const & unitCell) const;

      /**
      * Write an array of fields in basis format to a named file.
      *
      * This function opens an output file with the specified filename,
      * writes components in symmetry-adapted form to that file, and then
      * closes the file.
      *
      * \param filename  name of input file
      * \param fields  array of fields (symmetry adapted basis components)
      * \param unitCell  associated crystallographic unit cell
      */
      void writeFieldsBasis(std::string filename,
                            DArray< DArray<double> > const & fields,
                            UnitCell<D> const & unitCell) const;

      /**
      * Write a single field in basis format to an output stream.
      *
      * This function writes the field in symmetry adapted basis format.
      *
      * \param out  output stream (i.e., output file)
      * \param field  field to be written (symmetry adapted basis format)
      * \param unitCell  associated crystallographic unit cell
      */
      void writeFieldBasis(std::ostream& out,
                           DArray<double> const & field,
                           UnitCell<D> const & unitCell) const;

      /**
      * Write single field in basis format to a named file. 
      *
      * This function opens an output file with the specified filename,
      * writes the field in symmetry adapted basis format to that file,
      * and closes the file.
      *
      * \param filename  name of output file
      * \param field  field to be written (symmetry adapted basis format)
      * \param unitCell  associated crystallographic unit cell
      */
      void writeFieldBasis(std::string filename,
                           DArray<double> const & field,
                           UnitCell<D> const & unitCell) const;

      ///@}
      /// \name Field File IO - Real Space Grid Format
      ///@{

      /**
      * Read array of r-grid fields from an input stream.
      *
      * The capacity of array fields is equal to nMonomer, and element
      * fields[i] is the RFT associated with monomer type i.
      *
      * \param in  input stream (i.e., input file)
      * \param fields  array of RField fields (r-space grid)
      * \param unitCell  associated crystallographic unit cell
      * \return  true iff header has a space group (isSymmetric flag)
      */
      virtual
      bool readFieldsRGrid(std::istream& in,
                           DArray<RFT>& fields,
                           UnitCell<D> & unitCell) const = 0;

      /**
      * Read an array of r-grid fields from a named file.
      *
      * This function opens an input file with the specified filename,
      * reads fields in real-space grid format from that file, and 
      * then closes the file. 
      *
      * The capacity of array fields is equal to nMonomer, and element
      * fields[i] is the RFT associated with monomer type i.
      *
      * \param filename  name of input file
      * \param fields  array of RField fields (r-space grid)
      * \param unitCell  associated crystallographic unit cell
      * \return  true iff header has a space group (isSymmetric flag)
      */
      bool readFieldsRGrid(std::string filename,
                           DArray<RFT>& fields,
                           UnitCell<D> & unitCell) const;

      /**
      * Read data for array of r-grid fields, with no header section.
      *
      * This function reads the data section of the rgrid-field format, 
      * without first reading a header.
      *
      * \param in  input file stream
      * \param fields  array of RField fields (r-space grid)
      * \param nMonomer  number of monomer types
      */
      virtual
      void readFieldsRGridData(std::istream& in,
                               DArray<RFT>& fields,
                               int nMonomer) const = 0;

      /**
      * Read single r-grid field from an input stream.
      *
      * \param in  input stream (i.e., input file)
      * \param field  single field defined on r-space grid
      * \param unitCell  associated crystallographic unit cell
      * \return  true iff header has a space group (isSymmetric flag)
      */
      virtual
      bool readFieldRGrid(std::istream &in,
                          RFT & field,
                          UnitCell<D>& unitCell) const = 0;

      /**
      * Read single r-grid field from a named file.
      *
      * \param filename  name of input file
      * \param field  single field defined on r-space grid
      * \param unitCell  associated crystallographic unit cell
      * \return  true iff header has a space group (isSymmetric flag)
      */
      bool readFieldRGrid(std::string filename,
                          RFT & field,
                          UnitCell<D>& unitCell) const;

      /**
      * Write an arra of r-grid field to an output stream.
      *
      * \param out  output stream (i.e., output file)
      * \param fields  array of RField objects (fields on r-space grid)
      * \param unitCell  associated crystallographic unit cell
      * \param writeHeader  flag to write header of file if true
      * \param isSymmetric  Do fields have a space group symmetry ?
      * \param writeMeshSize Should mesh size be written at end of header?
      */
      virtual
      void writeFieldsRGrid(std::ostream& out,
                            DArray<RFT> const & fields,
                            UnitCell<D> const & unitCell,
                            bool writeHeader = true,
                            bool isSymmetric = true,
                            bool writeMeshSize = true) const = 0;

      /**
      * Write an array of r-grid fields to a named file.
      *
      * \param filename  name of output file
      * \param fields  array of RFT objects (fields on r-space grid)
      * \param unitCell  associated crystallographic unit cell
      * \param isSymmetric Do fields have a space group symmetry ?
      */
      void writeFieldsRGrid(std::string filename,
                            DArray<RFT> const & fields,
                            UnitCell<D> const & unitCell,
                            bool isSymmetric = true) const;

      /**
      * Write a single r-grid field to an an output stream.
      * 
      * \param out  output stream
      * \param field  field defined on r-space grid
      * \param unitCell  associated crystallographic unit cell
      * \param writeHeader  should a file header be written?
      * \param isSymmetric  Does the field have a space group symmetry?
      */
      virtual
      void writeFieldRGrid(std::ostream &out,
                           RFT const & field,
                           UnitCell<D> const & unitCell,
                           bool writeHeader = true,
                           bool isSymmetric = true) const = 0;

      /**
      * Write a single r-grid field to a named file. 
      *
      * \param filename  name of output file
      * \param field  field defined on r-space grid
      * \param unitCell  associated crystallographic unit cell
      * \param isSymmetric  Does the field have a space group symmetry?
      */
      void writeFieldRGrid(std::string filename,
                           RFT const & field,
                           UnitCell<D> const & unitCell,
                           bool isSymmetric = true) const;

      ///@}
      /// \name Field File IO - Fourier Space (K-Space) Grid Format
      ///@{

      /**
      * Read an array k-grid from an input stream.
      *
      * The capacity of the array is equal to nMonomer, and element
      * fields[i] is the discrete Fourier transform of the field for
      * monomer type i.
      *
      * \param in  input stream (i.e., input file)
      * \param fields  array of KFT fields (k-space grid)
      * \param unitCell  associated crystallographic unit cell
      */
      virtual
      void readFieldsKGrid(std::istream& in,
                           DArray<KFT>& fields,
                           UnitCell<D> & unitCell) const = 0;

      /**
      * Read an array k-grid from a named file. 
      *
      * This function opens a file with name filename, reads discrete
      * Fourier components (Dft) of fields from that file, and closes
      * the file.
      *
      * The capacity of the array is equal to nMonomer, and element
      * fields[i] is the discrete Fourier transform of the field for
      * monomer type i.
      *
      * \param filename  name of input file
      * \param fields  array of KFT fields (k-space grid)
      * \param unitCell  associated crystallographic unit cell
      */
      void readFieldsKGrid(std::string filename,
                           DArray<KFT>& fields,
                           UnitCell<D> & unitCell) const;

      /**
      * Write an array k-grid fields to a output stream.
      *
      * The capacity of the array fields is equal to nMonomer. Element
      * fields[i] is the discrete Fourier transform of the field for
      * monomer type i.
      *
      * \param out  output stream (i.e., output file)
      * \param fields  array of KFT fields
      * \param unitCell  associated crystallographic unit cell
      * \param isSymmetric  Does this field have space group symmetry?
      */
      virtual
      void writeFieldsKGrid(std::ostream& out,
                            DArray<KFT> const & fields,
                            UnitCell<D> const & unitCell,
                            bool isSymmetric = true) const = 0;

      /**
      * Write an array k-grid fields to a named file.
      *
      * This function opens a file with name filename, writes discrete
      * Fourier transform components (DFT) components of fields to that
      * file, and closes the file.
      *
      * \param filename  name of output file.
      * \param fields  array of KFT fields (k-space grid)
      * \param unitCell  associated crystallographic unit cell
      * \param isSymmetric  Does this field have space group symmetry?
      */
      void writeFieldsKGrid(std::string filename,
                           DArray<KFT> const & fields,
                           UnitCell<D> const & unitCell,
                           bool isSymmetric = true) const;

      ///@}
      /// \name Field Format Conversion
      ///@{

      /**
      * Convert a field from symmetrized basis to Fourier (k-grid) format.
      *
      * \param components coefficients of symmetry-adapted basis functions
      * \param dft discrete Fourier transform of a real field
      */
      virtual
      void convertBasisToKGrid(DArray<double> const & components,
                               KFT& dft) const = 0;

      /**
      * Convert an array of fields from basis to Fourier (k-grid) format.
      *
      * The in and out parameters are arrays of fields, in which element
      * number i is the field associated with monomer type i.
      *
      * \param in  fields expanded in symmetry-adapted Fourier basis
      * \param out  fields defined as discrete Fourier transforms (k-grid)
      */
      void convertBasisToKGrid(DArray< DArray<double> > const & in,
                               DArray<KFT>& out) const;

      /**
      * Convert a field from Fourier (k-grid) to symmetrized basis form.
      *
      * If the checkSymmetry parameter is true, this function checks if
      * the input field satisfies the space group symmetry to within a
      * tolerance given by the epsilon parameter, and prints a warning to
      * Log::file() if it does not.
      *
      * \param in  discrete Fourier transform (k-grid) of a field
      * \param out  components of field in asymmetry-adapted Fourier basis
      * \param checkSymmetry  flag indicating whether to check symmetry
      * \param epsilon  error tolerance for symmetry test (if any)
      */
      virtual
      void convertKGridToBasis(KFT const & in,
                               DArray<double> & out,
                               bool checkSymmetry = true,
                               double epsilon = 1.0e-8) const = 0;

      /**
      * Convert an array of fields from Fourier (k-grid) to basis format.
      *
      * If the checkSymmetry parameter is true, this function checks if
      * the input fields all satisfies the space group symmetry to within
      * a tolerance given by the parameter epsilon, and prints a warning
      * to Log::file() for each field that does not.
      *
      * \param in  fields defined as discrete Fourier transforms (k-grid)
      * \param out  components of fields in symmetry adapted basis
      * \param checkSymmetry  flag indicate whether to check symmetry
      * \param epsilon  error tolerance for symmetry test (if any)
      */
      void convertKGridToBasis(DArray<KFT> const & in,
                               DArray< DArray<double> > & out,
                               bool checkSymmetry = true,
                               double epsilon = 1.0e-8) const;

      /**
      * Convert a single field from basis to r-grid format.
      *
      * \param in  field in symmetry adapted basis form
      * \param out field defined on real-space grid
      */
      void convertBasisToRGrid(DArray<double> const & in,
                               RFT & out) const;

      /**
      * Convert an array of fields from basis to r-grid format.
      *
      * \param in  fields in symmetry adapted basis form
      * \param out fields defined on real-space grid
      */
      void convertBasisToRGrid(DArray< DArray<double> > const & in,
                               DArray<RFT> & out) const ;

      /**
      * Convert a field from spatial grid (r-grid) to symmetrized basis.
      *
      * If the boolean checkSymmetry parameter is true, this function
      * checks whether the the input field has the correct space group
      * symmetry to within an error threshhold given by the epsilon
      * parameter. If this error threshhold is exceeded, a warning is
      * written to Log::file().
      *
      * \param in  field defined on real-space grid
      * \param out  field in symmetry adapted basis form
      * \param checkSymmetry  if true, check space group symmetry
      * \param epsilon error threshhold for symmetry test
      */
      void convertRGridToBasis(RFT const & in,
                               DArray<double> & out,
                               bool checkSymmetry = true,
                               double epsilon = 1.0e-8) const;

      /**
      * Convert an array of fields from r-grid to basis format.
      *
      * If the boolean checkSymmetry parameter is true, this function
      * checks whether the input fields all have the correct space group
      * symmetry to within an error threshhold given by the epsilon
      * parameter. If this error threshhold is exceeded by an field, a
      * warning is written to Log::file().
      *
      * \param in  fields defined on real-space grid
      * \param out  fields in symmetry adapted basis form
      * \param checkSymmetry  if true, check space group symmetry
      * \param epsilon error threshhold for symmetry test
      */
      void convertRGridToBasis(DArray<RFT> const & in,
                               DArray< DArray<double> > & out,
                               bool checkSymmetry = true,
                               double epsilon = 1.0e-8) const;

      /**
      * Convert an array of field from k-grid to r-grid format.
      *
      * This function simply calls the inverse FFT for an array of fields.
      *
      * \param in  fields in discrete Fourier format (k-grid)
      * \param out  fields defined on real-space grid (r-grid)
      */
      void convertKGridToRGrid(DArray<KFT> const & in,
                               DArray<RFT> & out) const;

      /**
      * Convert a single field from k-grid to r-grid format.
      *
      * This function simply calls the inverse FFT for a single field.
      *
      * \param in  field in discrete Fourier format (k-grid)
      * \param out  field defined on real-space grid (r-grid)
      */
      void convertKGridToRGrid(KFT const & in,
                               RFT & out) const;

      /**
      * Convert fields from spatial grid (r-grid) to k-grid format.
      *
      * This function simply calls the forward FFT repeatedly for an 
      * array of fields.
      *
      * \param in  fields defined on real-space grid (r-grid)
      * \param out  fields in discrete Fourier format (k-grid)
      */
      void convertRGridToKGrid(DArray<RFT> const & in,
                               DArray<KFT> & out) const;

      /**
      * Convert a field from spatial grid (r-grid) to k-grid format.
      *
      * This function simply calls the forward FFT for a single field.
      *
      * \param in  field defined on real-space grid (r-grid)
      * \param out  field in discrete Fourier format (k-grid)
      */
      void convertRGridToKGrid(RFT const & in,
                               KFT & out) const;

      ///@}
      /// \name Field File Format Conversion
      ///@{

      /**
      * Convert a field file from basis format to r-grid format.
      *
      * This function reads a field file in basis format, converts the
      * fields to r-grid format, and writes the fields in r-grid format
      * to a different file.
      *
      * \param inFileName name of input file (basis format)
      * \param outFileName name of output file (r-grid format)
      */
      void convertBasisToRGrid(std::string const & inFileName,
                               std::string const & outFileName) const;

      /**
      * Convert a field file from r-grid to basis format.
      *
      * This function checks if the input fields have the declared space
      * group symmetry, and prints a warning if it detects deviations
      * that exceed some small threshhold, but proceeds to attempt the
      * conversion even if such an error is detected. Converting a field
      * that does not have the declared space group symmetry to basis
      * format is a destructive operation that modifies the field in
      * unpredictable ways.
      *
      * \param inFileName name of input file (r-grid format)
      * \param outFileName name of output file (basis format)
      */
      void convertRGridToBasis(std::string const & inFileName,
                               std::string const & outFileName) const;

      /**
      * Convert a field file from Fourier (k-grid) to r-grid format.
      *
      * \param inFileName name of input file (k-grid format)
      * \param outFileName name of output file (r-grid format)
      */
      void convertKGridToRGrid(std::string const & inFileName,
                               std::string const & outFileName) const;

      /**
      * Convert a field file from r-grid to Fourier (k-grid) format.
      *
      * \param inFileName name of input file (r-grid format)
      * \param outFileName name of output file (k-grid format)
      */
      void convertRGridToKGrid(std::string const & inFileName,
                               std::string const & outFileName) const;

      /**
      * Convert a field file from Fourier (k-grid) to basis format.
      *
      * This function checks if the input fields have the declared space
      * group symmetry, and prints a warning if it detects deviations
      * that exceed some small threshhold, but proceeds to attempt the
      * conversion even if such an error is detected. Converting a field
      * that does not have the declared space group symmetry to basis
      * format is a destructive operation that modifies the field in
      * unpredictable ways.
      *
      * \param inFileName name of input file (k-grid format)
      * \param outFileName name of output file (basis format)
      */
      void convertKGridToBasis(std::string const & inFileName,
                               std::string const & outFileName) const;

      /**
      * Convert a field file from basis to Fourier (k-grid) format.
      *
      * \param inFileName name of input file (basis format)
      * \param outFileName name of output file (k-grid format)
      */
      void convertBasisToKGrid(std::string const & inFileName,
                               std::string const & outFileName) const;

      ///@}
      /// \name Test Space Group Symmetry
      ///@{

      /**
      * Check if a k-grid field has the declared space group symmetry.
      *
      * This function checks whether the discrete Fourier transform of
      * a real field is invariant under all the symmetry operations of 
      * a declared space group to within an specfiied error threshhold,
      * given by parameter epsilon. If the parameter verbose is true 
      * and the deviation from symmetry exceeds the error threshhold, 
      * errors are written to Log::file().
      *
      * \param in field in real space grid (r-grid) format
      * \param epsilon error threshold used to test for symmetry
      * \param verbose  if true, write error to Log::file()
      *
      * \return true if the field is symmetric, false otherwise
      */
      virtual
      bool hasSymmetry(KFT const & in, 
                       double epsilon = 1.0e-8,
                       bool verbose = true) const = 0;

      /**
      * Check if an r-grid field has the declared space group symmetry.
      *
      * This function checks whether a field defined on the nodes of a
      * regular real-space grid is invariant under all the symmetry 
      * operations of a declared space group to within an error 
      * threshhold given by function parameter epsilon. If parameter
      * verbose is true and the deviation from symmetry exceeds the 
      * error threshhold, errors are written to Log::file().
      *
      * \param in field in real space grid (r-grid) format
      * \param epsilon error threshold used to test for symmetry
      * \param verbose  if true, write error to Log::file()
      *
      * \return true if the field is symmetric, false otherwise
      */
      bool hasSymmetry(RFT const & in, 
                       double epsilon = 1.0e-8,
                       bool verbose = true) const;

      /**
      * Check if an r-grid field file has declared space group symmetry.
      *
      * \param inFileName name of input r-grid field file
      * \param epsilon error threshold used when testing for symmetry
      * \return true if fields all have symmetry, false otherwise
      */
      bool hasSymmetry(std::string const & inFileName,
                       double epsilon = 1.0E-8) const;

      ///@}
      /// \name Test Field Equality
      ///@{

      /**
      * Compare arrays of fields in basis format, write a report to Log file.
      *
      * Outputs maximum and root-mean-squared differences to the
      * standard Log file.
      *
      * \param field1  first array of fields (basis format)
      * \param field2  second array of fields (basis format)
      */
      void compareFieldsBasis(DArray< DArray<double> > const & field1,
                              DArray< DArray<double> > const & field2) 
      const;

      /**
      * Compare two r-grid field files, write a report to Log file.
      *
      * \param filename1  name of 1st field file
      * \param filename2  name of 2nd field file
      */
      void compareFieldsBasis(std::string const & filename1,
                              std::string const & filename2) const;

      /**
      * Compare two fields in r-grid format, write a report to Log file.
      *
      * Outputs maximum and root-mean-squared differences to the
      * standard Log file.
      *
      * \param field1  first array of fields (r-grid format)
      * \param field2  second array of fields (r-grid format)
      */
      virtual
      void compareFieldsRGrid(DArray< RFT > const & field1,
                              DArray< RFT > const & field2) const = 0;
  
      /**
      * Compare two r-grid field files, write a report to Log file.
      *
      * \param filename1  name of 1st field file
      * \param filename2  name of 2nd field file
      */
      void compareFieldsRGrid(std::string const & filename1,
                              std::string const & filename2) const; 

      ///@}
      /// \name Field Scaling (Multiplication by a Scalar)
      ///@{

      /**
      * Multiply a single field in basis format by a real scalar.
      *
      * This function takes a single real periodic field and multiplies all
      * components in place by a common real factor, thereby modifying the 
      * input.
      * 
      * \param field  field in basis form to be rescaled
      * \param factor  factor by which to multiply every field element
      */
      virtual
      void scaleFieldBasis(DArray<double>& field, 
                           double factor) const;

      /**
      * Multiply an array of fields in basis format by a real scalar.
      *
      * This function takes an array of real periodic fields and multiplies
      * all components in place by a common real scalar, thereby modifying 
      * the input.
      * 
      * \param fields  array of fields in basis form to be rescaled
      * \param factor  factor by which to multiply every field element
      */
      void scaleFieldsBasis(DArray< DArray<double> >& fields, 
                            double factor) const;

      /**
      * Multiply all fields in a basis field file by a scalar.
      *
      * \param inFileName  name of input field file
      * \param outFileName  name of file for rescaled output fields
      * \param factor  factor by which to multiply all field elements
      */
      void scaleFieldsBasis(std::string const & inFileName,
                            std::string const & outFileName,
                            double factor) const;

      /**
      * Multiply a single field in r-grid format by a scalar.
      *
      * This function takes a single real periodic field and multiplies all
      * elements in place by a common real factor, thereby modifying the 
      * input.
      * 
      * \param field  RField (r-space) field to be rescaled
      * \param factor  factor by which to multiply every field element
      */
      virtual
      void scaleFieldRGrid(RFT& field, double factor) const = 0;

      /**
      * Scale an array of r-grid fields by a scalar.
      *
      * This function takes an array of real periodic fields and multiplies
      * all elements in place by a common real scalar, thereby modifying 
      * the input.
      * 
      * \param fields  array of RField (r-space) fields to be rescaled
      * \param factor  factor by which to multiply every field element
      */
      void scaleFieldsRGrid(DArray<RFT> & fields, double factor) const;

      /**
      * Multiply all fields in an r-grid field file by a scalar.
      *
      * Read a set of fields from a file, multiply fields by a constant, 
      * and write rescaled fields to a different file. 
      *
      * \param inFileName  name of input field file
      * \param outFileName  name of file for rescaled output fields
      * \param factor  factor by which to multiply all field elements
      */
      void scaleFieldsRGrid(std::string const & inFileName,
                            const std::string & outFileName,
                            double factor) const;

      ///@}
      /// \name Replicate Unit Cell
      ///@{

      /**
      * Write r-grid fields in a replicated unit cell to std::ostream.
      *
      * This function takes an input array of periodic fields and outputs
      * them within an expanded unit cell in which the original input 
      * unit cell has been replicated a specified number of times in each 
      * direction. Results are written to an std::ostream output stream.
      *
      * Element i of the replicas IntVec<D> parameter contains the number
      * of unit cell replicas along direction i. 
      *
      * \param out  output stream (i.e., output file)
      * \param fields  array of RField (r-space) fields to be replicated
      * \param unitCell  original crystallographic unit cell
      * \param replicas  number of unit cell replicas in each direction
      */
      virtual
      void replicateUnitCell(std::ostream& out,
                             DArray<RFT> const & fields,
                             UnitCell<D> const & unitCell,
                             IntVec<D> const & replicas) const = 0;

      /**
      * Write r-grid fields in a replicated unit cell to a named file.
      *
      * This function opens output file filename, writes fields within 
      * a replicated unit cell to the file, and closes the file. See
      * documentation of the overloaded function of the same name with 
      * a std::ostream parameter, which is called internally. 
      *
      * \param filename  output file name
      * \param fields  array of RField fields (r-space grid) needs
      * \param unitCell  original crystallographic unit cell
      * \param replicas  number of unit cell replicas in each direction
      */
      void replicateUnitCell(std::string filename,
                             DArray<RFT> const & fields,
                             UnitCell<D> const & unitCell,
                             IntVec<D> const & replicas) const;

      /**
      * Write replicated fields read from one file to another.
      *
      * This function reads a field file in r-grid format, and writes
      * replicated fields to another file.
      *
      * \param inFileName  name of input field file
      * \param outFileName  name of output field file
      * \param replicas  the number of replicas in each D direction
      */
      void replicateUnitCell(std::string const & inFileName,
                             std::string const & outFileName,
                             IntVec<D> const & replicas) const;

      ///@}
      /// \name Expand Spatial Dimension 
      ///@{

      /**
      * Increase D for an array of r-grid fields, write to ostream.
      *
      * This function is used for template dimension D < 3, and allows a
      * 1D or 2D field to be transformed into a field defined in a higher
      * dimension (2D or 3D) space in which field values are indepedendent
      * of the values of coordinates associated with the added dimensions. 
      * For example, it can transform a lamellar field computed with D=1 
      * to a lamellar field field on a 3D grid (d=3), written in a format 
      * that can be read by pscf_pc when invoked with D=3.
      *
      * \param out  output stream (i.e., output file)
      * \param fields  input array of RField fields (r-space grid)
      * \param unitCell  original crystallographic unit cell
      * \param d  expanded dimension (greater than D)
      * \param newGridDimensions number of grid points in added dimensions
      */
      virtual
      void expandRGridDimension(std::ostream &out,
                                DArray<RFT > const & fields,
                                UnitCell<D> const & unitCell,
                                int d,
                                DArray<int> const& newGridDimensions) 
      const = 0;

      /**
      * Increase D for an array of r-grid fields, write to a named file.
      *
      * This function opens an output file with the specified filename,
      * writes expanded fields in RField<d> real-space grid format to 
      * that file, and then closes the file. The overloaded function of 
      * the same name with an std::ostream parameter is called internally.
      *
      * \param filename  name of output file
      * \param fields  input array of RFT objects (r-space grid)
      * \param unitCell  original crystallographic unit cell
      * \param d  expanded dimension (greater than D)
      * \param newGridDimensions  number of grid points in added dimensions
      */
      void expandRGridDimension(std::string filename,
                                DArray<RFT > const & fields,
                                UnitCell<D> const & unitCell,
                                int d,
                                DArray<int> newGridDimensions) const;

      /**
      * Increase D for an r-grid field file.
      *
      * This function reads an array of fields in D-dimensional space 
      * from an r-grid field file and writes expanded fields in 
      * d-dimensional space to another file.
      *
      * \param inFileName filename name of input field file
      * \param outFileName filename name of output field file
      * \param d  intended dimensions (d > D)
      * \param newGridDimensions number of grid points in added dimensions
      */
      void expandRGridDimension(std::string const & inFileName,
                                std::string const & outFileName,
                                int d,
                                DArray<int> newGridDimensions) const;

      ///@}
      /// \name Field File IO Utilities
      ///@{

      /**
      * Reader header of field file (fortran PSCF format)
      *
      * This reads the common part of the header for all PSCF field file
      * formats. This contains the dimension of space, the lattice 
      * system, a list of unit cell parameters, the space group name as 
      * an optional parameter, and the number of monomer types. The unit
      * cell data is used to update a UnitCell<D> that is passed as a 
      * parameter.
      *
      * Consistency checks:
      *
      * The value of "dim" in the header file must match the template
      * parameter D.  If the UnitCell<D> object passed to this function 
      * already contains a non-null lattice type, it must match the 
      * lattice system in the header file. 
      *
      * The isSymmetric parameter is set to true on return if a group 
      * name is found in the header. Presence of a group name parameter in
      * the header is optional. If the header does contain a group name
      * and a group name was declared in the parameter file (i.e., if 
      * hasGroup() is true), then these group names must match or an
      * Exception is thrown. 
      *
      * This function does NOT check if the number of monomers declared 
      * in the field file header matches the value of nMonomer set by 
      * the setNMonomer member function. On return, parameter nMonomer
      * contains the number of monomer types declared in the field file
      * header.
      *
      * If a space group was defined in the parameter file but the 
      * associated Basis object is not been initialized, this function will
      * initialize the basis by calling Basis<D>::makeBasis via a private
      * pointer, using the unit cell parameters found in the file header.
      * This function thus can modify the associated Basis object as a 
      * side effect (even though this function is marked const).  Because 
      * all member functions that read entire field files call this 
      * function to read the file header, the same statement about the
      * Basis also applies to all member functions that read fields from
      * a file.
      *
      * \param in  input stream (i.e., input file)
      * \param nMonomer  number of fields in the field file (output)
      * \param unitCell  associated crystallographic unit cell (output)
      * \param isSymmetric Is there a group name in the header? (output)
      */
      void readFieldHeader(std::istream& in, 
                           int& nMonomer,
                           UnitCell<D> & unitCell,
                           bool & isSymmetric) const;

      /**
      * Write header for field file (fortran pscf format)
      *
      * \param out  output stream (i.e., output file)
      * \param nMonomer  number of monomer types or fields
      * \param unitCell  associated crystallographic unit cell
      * \param isSymmetric Do the fields have a space group symmetry?
      */
      void writeFieldHeader(std::ostream& out, 
                            int nMonomer,
                            UnitCell<D> const & unitCell,
                            bool isSymmetric = true) const;

      ///@}
      /// \name Accessor functions (const refererence)
      ///@{

      /// Get spatial discretization mesh by const reference
      Mesh<D> const & mesh() const
      {
         UTIL_ASSERT(meshPtr_);
         return *meshPtr_;
      }

      /// Get the associated Basis by const reference
      Basis<D> const & basis() const
      {
         UTIL_ASSERT(basisPtr_);
         return *basisPtr_;
      }

      /// Get associated FileMaster by const reference
      FileMaster const & fileMaster() const
      {
         UTIL_ASSERT(fileMasterPtr_);
         return *fileMasterPtr_;
      }

      ///@}

   protected:

      /// Get the lattice type enum value by const reference
      typename UnitCell<D>::LatticeSystem const & lattice() const
      {
         UTIL_ASSERT(latticePtr_);
         return *latticePtr_;
      }

      /// Has a space group been declared externally ?
      bool hasGroup() const
      {
         UTIL_ASSERT(hasGroupPtr_);
         return *hasGroupPtr_;
      }

      /// Get associated group name string by const reference
      std::string const & groupName() const
      {
         UTIL_ASSERT(groupNamePtr_);
         return *groupNamePtr_;
      }

      /// Get associated SpaceGroup<D> by const reference
      SpaceGroup<D> const & group() const
      {
         UTIL_ASSERT(groupPtr_);
         return *groupPtr_;
      }

      /// Get FFT object by const reference
      FFT const & fft() const
      {
         UTIL_ASSERT(fftPtr_);
         return *fftPtr_;
      }

      /**
      * Check if r-grid workspace is allocated, allocate if necessary.
      */
      void checkAllocateRGrid() const;

      /**
      * Check if k-grid workspace is allocated, allocate if necessary.
      */
      void checkAllocateKGrid() const;

      /**
      * Check if basis workspace is allocated, allocate if necessary.
      *
      * If the basis is not initialized, this function peeks at the 
      * header of the field file to initialize it, and thus obtain the
      * number of basis functions.
      *
      * \param inFileName name of field file 
      */
      void checkAllocateBasis(std::string const & inFileName) const;

   private:

      /// Pointer to spatial discretization mesh
      Mesh<D> const * meshPtr_;

      /// Pointer to FFT object
      FFT const * fftPtr_;

      /// Pointer to lattice system
      typename UnitCell<D>::LatticeSystem const * latticePtr_;

      /// Pointer to boolean hasGroup (true if a space group is known)
      bool const * hasGroupPtr_;

      /// Pointer to group name string
      std::string const * groupNamePtr_;

      /// Pointer to a SpaceGroup object
      SpaceGroup<D> const * groupPtr_;

      /// Pointer to a Basis object
      Basis<D> * basisPtr_;

      /// Pointer to Filemaster (holds paths to associated I/O files)
      FileMaster const * fileMasterPtr_;

      /// Number of monomer types
      int nMonomer_;

      // Mutable work space used by functions that read field files

      /// Work array of field coefficients for all monomer types.
      mutable DArray< DArray<double> > tmpFieldsBasis_;

      /// Work array of fields on real space grid.
      mutable DArray<RFT> tmpFieldsRGrid_;

      /// Work array of fields on Fourier grid (k-grid).
      mutable DArray<KFT> tmpFieldsKGrid_;

      /// K-grid work space (single field)
      mutable KFT workDft_;

      /// Is tmpFieldsBasis_ allocated?
      mutable bool isAllocatedBasis_;

      /// Is tmpFieldsRGrid_ allocated?
      mutable bool isAllocatedRGrid_;

      /// Is tmpFieldsKGrid_ allocated?
      mutable bool isAllocatedKGrid_;

   };

} // namespace Prdc
} // namespace Pscf
#endif

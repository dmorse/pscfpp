#ifndef RPC_FIELD_IO_H
#define RPC_FIELD_IO_H

/*
* PSCF - Polymer Self-Consistent Field Theory
*
* Copyright 2016 - 2022, The Regents of the University of Minnesota
* Distributed under the terms of the GNU General Public License.
*/

#include <prdc/cpu/RFieldDft.h>            // data member
#include <prdc/crystal/UnitCell.h>         // nested LatticeSystem enum
#include <util/containers/DArray.h>        // nested function parameter

// Forward declarations for classes used only via references or pointers
namespace Util {
   class FileMaster;
}
namespace Pscf {
   template <int D> class Mesh;
   namespace Prdc {
      template <int D> class Basis;
      template <int D> class SpaceGroup;
      namespace Cpu {
         template <int D> class FFT;
         template <int D> class RField;
      }
   }
}

namespace Pscf {
namespace Rpc {

   using namespace Util;
   using namespace Pscf;
   using namespace Pscf::Prdc;
   using namespace Pscf::Prdc::Cpu;

   /**
   * File input/output operations and format conversions for fields.
   *
   * This class provides functions to read and write arrays that contain
   * fields in any of three representations (symmetry-adapted basis,
   * r-space grid, or Fourier k-space grid), and to convert among these
   * representations. The member functions that implement IO operations
   * define the file formats for these field representations.
   *
   * \ingroup Rpc_Field_Module
   */
   template <int D>
   class FieldIo
   {

   public:

      /// \name Construction, Initialization and Destruction
      ///@{

      /**
      * Constructor.
      */
      FieldIo();

      /**
      * Destructor.
      */
      ~FieldIo();

      /**
      * Create association with other objects in parent Domain.
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
                     FFT<D> const & fft,
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

      ///@}
      /// \name Field File IO - Symmetry Adapted Basis Format
      ///@{

      /**
      * Read concentration or chemical potential fields from file.
      *
      * This function reads fields in a symmetry adapted basis from
      * input stream in.
      *
      * The capacity of DArray fields is equal to nMonomer, and element
      * fields[i] is a DArray containing components of the field
      * associated with monomer type i.
      *
      * \param in  input stream (i.e., input file)
      * \param fields  array of fields (symmetry adapted basis components)
      * \param unitCell  associated crystallographic unit cell
      */
      void
      readFieldsBasis(std::istream& in, DArray< DArray<double> > & fields,
                      UnitCell<D> & unitCell) const;

      /**
      * Read concentration or chemical potential components from file.
      *
      * This function opens an input file with the specified filename,
      * reads components in symmetry-adapted form from that file, and
      * closes the file.
      *
      * \param filename  name of input file
      * \param fields  array of fields (symmetry adapted basis components)
      * \param unitCell  associated crystallographic unit cell
      */
      void readFieldsBasis(std::string filename,
                           DArray< DArray<double> > & fields,
                           UnitCell<D> & unitCell) const;

      /**
      * Read single concentration or chemical potential field from file.
      *
      * This function reads the field in symmetry adapted basis format
      * from input stream in.
      *
      * \param in  input stream (i.e., input file)
      * \param field  array to store the field (basis format)
      * \param unitCell  associated crystallographic unit cell
      */
      void readFieldBasis(std::istream& in, DArray<double>& field,
                          UnitCell<D> & unitCell) const;

      /**
      * Read single concentration or chemical potential field from file.
      *
      * This function opens an input file with the specified filename,
      * reads field in symmetry adapted basis format from that file, and
      * closes the file.
      *
      * \param filename  name of input file
      * \param field  array to store the field (basis format)
      * \param unitCell  associated crystallographic unit cell
      */
      void readFieldBasis(std::string filename, DArray<double>& field,
                          UnitCell<D> & unitCell) const;

      /**
      * Write single concentration or chemical potential field to file.
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

      /**
      * Write concentration or chemical potential field components to file.
      *
      * This function writes components in a symmetry adapted basis.
      *
      * \param out  output stream (i.e., output file)
      * \param fields  array of fields (symmetry adapted basis components)
      * \param unitCell  associated crystallographic unit cell
      */
      void writeFieldsBasis(std::ostream& out,
                            DArray< DArray<double> > const & fields,
                            UnitCell<D> const & unitCell) const;

      /**
      * Write concentration or chemical potential field components to file.
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
      * Write single concentration or chemical potential field to output
      * stream out.
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

      ///@}
      /// \name Field File IO - Real Space Grid Format
      ///@{

      /**
      * Read array of RField objects (fields on r-space grid) from istream.
      *
      * The capacity of array fields is equal to nMonomer, and element
      * fields[i] is the RField<D> associated with monomer type i.
      *
      * \param in  input stream (i.e., input file)
      * \param fields  array of RField fields (r-space grid)
      * \param unitCell  associated crystallographic unit cell
      */
      void readFieldsRGrid(std::istream& in,
                           DArray< RField<D> >& fields,
                           UnitCell<D> & unitCell) const;

      /**
      * Read array of RField objects (fields on r-space grid) from file.
      *
      * The capacity of array fields is equal to nMonomer, and element
      * fields[i] is the RField<D> associated with monomer type i.
      *
      * This function opens an input file with the specified filename,
      * reads fields in RField<D> real-space grid format from that file,
      * and then closes the file.
      *
      * \param filename  name of input file
      * \param fields  array of RField fields (r-space grid)
      * \param unitCell  associated crystallographic unit cell
      */
      void readFieldsRGrid(std::string filename,
                           DArray< RField<D> >& fields,
                           UnitCell<D> & unitCell) const;

      /**
      * Read data for array of r-grid fields, with no header section.
      *
      * This function reads the data section of the rgrid-field format, with
      * with no header.
      *
      * \param in  input file stream
      * \param fields  array of RField fields (r-space grid)
      * \param nMonomer  number of monomer types
      */
      void readFieldsRGridData(std::istream& in,
                               DArray< RField<D> >& fields,
                               int nMonomer) const;

      /**
      * Read single RField (field on an r-space grid) from an istream.
      *
      * \param in  input stream (i.e., input file)
      * \param field  fields defined on r-space grid
      * \param unitCell  associated crystallographic unit cell
      */
      void readFieldRGrid(std::istream &in,
                           RField<D> & field,
                           UnitCell<D>& unitCell) const;

      /**
      * Read single RField (field on an r-space grid) from named file.
      *
      * This function opens an input file with the specified filename,
      * reads a field in RField<D> real-space grid format, and closes
      * the file.
      *
      * \param filename  name of input file
      * \param field  fields defined on r-space grid
      * \param unitCell  associated crystallographic unit cell
      */
      void readFieldRGrid(std::string filename,
                           RField<D> & field,
                           UnitCell<D>& unitCell) const;

      /**
      * Write array of RField objects (fields on r-space grid) to ostream.
      *
      * \param out  output stream (i.e., output file)
      * \param fields  array of RField objects (fields on r-space grid)
      * \param unitCell  associated crystallographic unit cell
      * \param writeHeader  flag to write header of file if true
      * \param isSymmetric  Do fields have a space group symmetry ?
      * \param writeMeshSize Should mesh size be written at end of header?
      */
      void writeFieldsRGrid(std::ostream& out,
                            DArray< RField<D> > const & fields,
                            UnitCell<D> const & unitCell,
                            bool writeHeader = true,
                            bool isSymmetric = true,
                            bool writeMeshSize = true) const;

      /**
      * Write array of RField objects (fields on an r-space grid) to file.
      *
      * This function opens an output file with the specified filename,
      * writes fields in RField<D> real-space grid format to that file,
      * and then closes the file.
      *
      * \param filename  name of output file
      * \param fields  array of RField<D> objects (fields on r-space grid)
      * \param unitCell  associated crystallographic unit cell
      * \param isSymmetric Do fields have a space group symmetry ?
      */
      void writeFieldsRGrid(std::string filename,
                            DArray< RField<D> > const & fields,
                            UnitCell<D> const & unitCell,
                            bool isSymmetric = true) const;

      /**
      * Write a single RField (field on an r-space grid) to ostream.
      *
      * \param out  output stream
      * \param field  field defined on r-space grid
      * \param unitCell  associated crystallographic unit cell
      * \param writeHeader  should a file header be written?
      * \param isSymmetric  Does the field have a space group symmetry?
      */
      void writeFieldRGrid(std::ostream &out,
                           RField<D> const & field,
                           UnitCell<D> const & unitCell,
                           bool writeHeader = true,
                           bool isSymmetric = true) const;

      /**
      * Write a single RField (fields on an r-space grid) to a file.
      *
      * This function opens an output file with the specified filename,
      * write a field in RField<D> real-space grid format to that file,
      * and then closes the file.
      *
      * \param filename  name of output file
      * \param field  field defined on r-space grid
      * \param unitCell  associated crystallographic unit cell
      * \param isSymmetric  Does the field have a space group symmetry?
      */
      void writeFieldRGrid(std::string filename,
                           RField<D> const & field,
                           UnitCell<D> const & unitCell,
                           bool isSymmetric = true) const;

      ///@}
      /// \name Field File IO - Fourier Space (K-Space) Grid Format
      ///@{

      /**
      * Read array of RFieldDft objects (k-space fields) from istream.
      *
      * The capacity of the array is equal to nMonomer, and element
      * fields[i] is the discrete Fourier transform of the field for
      * monomer type i.
      *
      * \param in  input stream (i.e., input file)
      * \param fields  array of RFieldDft fields (k-space grid)
      * \param unitCell  associated crystallographic unit cell
      */
      void readFieldsKGrid(std::istream& in,
                           DArray< RFieldDft<D> >& fields,
                           UnitCell<D> & unitCell) const;

      /**
      * Read array of RFieldDft objects (k-space fields) from file.
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
      * \param fields  array of RFieldDft fields (k-space grid)
      * \param unitCell  associated crystallographic unit cell
      */
      void readFieldsKGrid(std::string filename,
                           DArray< RFieldDft<D> >& fields,
                           UnitCell<D> & unitCell) const;

      /**
      * Write array of RFieldDft objects (k-space fields) to file.
      *
      * The capacity of the array fields is equal to nMonomer. Element
      * fields[i] is the discrete Fourier transform of the field for
      * monomer type i.
      *
      * \param out  output stream (i.e., output file)
      * \param fields  array of RFieldDft fields
      * \param unitCell  associated crystallographic unit cell
      * \param isSymmetric  Does this field have space group symmetry?
      */
      void writeFieldsKGrid(std::ostream& out,
                            DArray< RFieldDft<D> > const & fields,
                            UnitCell<D> const & unitCell,
                            bool isSymmetric = true) const;

      /**
      * Write array of RFieldDft objects (k-space fields) to a file.
      *
      * This function opens a file with name filename, writes discrete
      * Fourier transform components (DFT) components of fields to that
      * file, and closes the file.
      *
      * \param filename  name of output file.
      * \param fields  array of RFieldDft fields (k-space grid)
      * \param unitCell  associated crystallographic unit cell
      * \param isSymmetric  Does this field have space group symmetry?
      */
      void writeFieldsKGrid(std::string filename,
                           DArray< RFieldDft<D> > const & fields,
                           UnitCell<D> const & unitCell,
                           bool isSymmetric = true) const;

      ///@}
      /// \name Field File IO Utilities
      ///@{

      /**
      * Reader header of field file (fortran pscf format)
      *
      * This reads the common part of the header for all field file
      * formats. This contains the dimension of space, the lattice
      * system, a list of unit cell parameters, the space group name
      * as an optional parameter, and the number of monomer types.
      * The unit cell data is read into the associated UnitCell<D>,
      * which is thus updated.
      *
      * The value of "dim" in the header file must match the template
      * parameter D, or an Exception is thrown.  If the UnitCell<D>
      * object passed to this function already contains a non-null
      * lattice type, it must match the lattice system in the header
      * file, or an Exception is thrown.
      *
      * Presence of a group name parameter in the header is optional
      * if a group name was declared in the parameter file (i.e., if
      * FieldIo::hasGroup() is true) and forbidden otherwise. If the
      * header contains a group name, it must match the value given
      * in the parameter file (accessed as FieldIo::groupName()), or
      * an Exception is thrown.
      *
      * If the header contains the appropriate group name but the
      * associated symmetry-adapted Fourier basis has not been
      * initialized, this function will attempt to initialize it using
      * the unit cell parameters read from file.
      *
      * \param in  input stream (i.e., input file)
      * \param nMonomer  number of fields in the field file (output)
      * \param unitCell  associated crystallographic unit cell (output)
      * \param isSymmetric Is there a group name in the header? (output)
      */
      void readFieldHeader(std::istream& in, int& nMonomer,
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
      void writeFieldHeader(std::ostream& out, int nMonomer,
                            UnitCell<D> const & unitCell,
                            bool isSymmetric = true) const;

      ///@}
      /// \name Field Format Conversion
      ///@{

      /**
      * Convert a field from symmetrized basis to Fourier grid (k-grid).
      *
      * \param components coefficients of symmetry-adapted basis functions
      * \param dft discrete Fourier transform of a real field
      */
      void convertBasisToKGrid(DArray<double> const & components,
                               RFieldDft<D>& dft) const;

      /**
      * Convert fields from symmetrized basis to Fourier grid (k-grid).
      *
      * The in and out parameters are arrays of fields, in which element
      * number i is the field associated with monomer type i.
      *
      * \param in  fields expanded in symmetry-adapted Fourier basis
      * \param out  fields defined as discrete Fourier transforms (k-grid)
      */
      void convertBasisToKGrid(DArray< DArray<double> > const & in,
                               DArray< RFieldDft<D> >& out) const;

      /**
      * Convert a field from Fourier grid (k-grid) to symmetrized basis.
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
      void convertKGridToBasis(RFieldDft<D> const & in,
                               DArray<double> & out,
                               bool checkSymmetry = true,
                               double epsilon = 1.0e-8) const;

      /**
      * Convert fields from Fourier grid (k-grid) to symmetrized basis.
      *
      * The in and out parameters are each an array of fields, in which
      * element i is the field associated with monomer type i.
      *
      * If the checkSymmetry parameter is true, this function checks if
      * the input fields all satisfies the space group symmetry to within
      * a tolerance given by the parameter epsilon, and prints a warning
      * to Log::file() if one or more fields do not.
      *
      * \param in  fields defined as discrete Fourier transforms (k-grid)
      * \param out  components of fields in symmetry adapted basis
      * \param checkSymmetry  flag indicate whether to check symmetry
      * \param epsilon  error tolerance for symmetry test (if any)
      */
      void convertKGridToBasis(DArray< RFieldDft<D> > const & in,
                               DArray< DArray<double> > & out,
                               bool checkSymmetry = true,
                               double epsilon = 1.0e-8) const;

      /**
      * Convert a field from symmetrized basis to spatial grid (r-grid).
      *
      * \param in  field in symmetry adapted basis form
      * \param out field defined on real-space grid
      */
      void convertBasisToRGrid(DArray<double> const & in,
                               RField<D> & out) const;

      /**
      * Convert fields from symmetrized basis to spatial grid (r-grid).
      *
      * \param in  fields in symmetry adapted basis form
      * \param out fields defined on real-space grid
      */
      void convertBasisToRGrid(DArray< DArray<double> > const & in,
                               DArray< RField<D> > & out) const ;

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
      void convertRGridToBasis(RField<D> const & in,
                               DArray<double> & out,
                               bool checkSymmetry = true,
                               double epsilon = 1.0e-8) const;

      /**
      * Convert fields from spatial grid (r-grid) to symmetrized basis.
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
      void convertRGridToBasis(DArray< RField<D> > const & in,
                               DArray< DArray<double> > & out,
                               bool checkSymmetry = true,
                               double epsilon = 1.0e-8) const;

      /**
      * Convert fields from k-grid (DFT) to real space (r-grid) format.
      *
      * This function simply calls the inverse FFT for an array of fields.
      * The inverse FFT provided by the underlying FFTW library overwrites
      * its input, which is why field parameter "in" is passed as a
      * non-const reference.
      *
      * \param in  fields in discrete Fourier format (k-grid)
      * \param out  fields defined on real-space grid (r-grid)
      */
      void convertKGridToRGrid(DArray< RFieldDft<D> > & in,
                               DArray< RField<D> > & out) const;

      /**
      * Convert a field from k-grid (DFT) to real space (r-grid) format.
      *
      * This function simply calls the inverse FFT for a single field.
      * The inverse FFT provided by the underlying FFTW library overwrites
      * its input, which is why argument "in" a non-const reference.
      *
      * \param in  field in discrete Fourier format (k-grid)
      * \param out  field defined on real-space grid (r-grid)
      */
      void convertKGridToRGrid(RFieldDft<D> & in,
                               RField<D> & out) const;

      /**
      * Convert fields from spatial grid (r-grid) to k-grid format.
      *
      * \param in  fields defined on real-space grid (r-grid)
      * \param out  fields in discrete Fourier format (k-grid)
      */
      void convertRGridToKGrid(DArray< RField<D> > const & in,
                               DArray< RFieldDft<D> > & out) const;

      /**
      * Convert a field from spatial grid (r-grid) to k-grid format.
      *
      * \param in  field defined on real-space grid (r-grid)
      * \param out  field in discrete Fourier format (k-grid)
      */
      void convertRGridToKGrid(RField<D> const & in,
                               RFieldDft<D> & out) const;

      ///@}
      /// \name Test Space Group Symmetry
      ///@{

      /**
      * Check if an r-grid field has the declared space group symmetry.
      *
      * This function checks whether a field defined on the nodes of a
      * regular real-space grid satisfies all the symmetries of a space
      * group to within an error threshhold given by parameter epsilon.
      * If parameter verbose is true and the deviation from symmetry
      * exceeds the error threshhold, errors are written to Log::file().
      *
      * \param in field in real space grid (r-grid) format
      * \param epsilon error threshold used to test for symmetry
      * \param verbose  if true, write error to Log::file()
      * \return true if the field is symmetric, false otherwise
      */
      bool hasSymmetry(RField<D> const & in, double epsilon = 1.0e-8,
                       bool verbose = true) const;

      /**
      * Check if a k-grid field has the declared space group symmetry.
      *
      * This function checks whether the discrete Fourier transform of
      * a real field satisfies all the symmetries of a space group to
      * within an error threshhold given by parameter epsilon. If the
      * parameter verbose is true and the deviation from symmetry
      * exceeds the error threshhold, errors are written to Log::file().
      *
      * \param in field in real space grid (r-grid) format
      * \param epsilon error threshold used to test for symmetry
      * \param verbose  if true, write error to Log::file()
      * \return true if the field is symmetric, false otherwise
      */
      bool hasSymmetry(RFieldDft<D> const & in, double epsilon = 1.0e-8,
                       bool verbose = true) const;

      ///@}
      /// \name Grid Manipulation Utilities
      ///@{

      /**
      * Expand dimension of an array of r-grid fields, write to ostream.
      *
      * This function is used for template dimension D < 3, and allows a
      * 1D or 2D field to be expanded into a higher dimensional (2D or 3D)
      * field in which field values are independent of the values of
      * coordinates associated with the added dimensions. For example, 
      * it can output a lamellar computed with D=1 on a 3D grid (d=3)
      * in a format that can be read by pscf_pc when invoked with D=3.
      *
      * \param out  output stream (i.e., output file)
      * \param fields  input array of RField fields (r-space grid)
      * \param unitCell  original crystallographic unit cell
      * \param d  expanded dimension (greater than D)
      * \param newGridDimensions number of grid points in added dimensions
      */
      void expandFieldsDimension(std::ostream &out,
                                 DArray<RField<D> > const & fields,
                                 UnitCell<D> const & unitCell,
                                 int d,
                                 DArray<int> newGridDimensions) const;

      /**
      * Expand dimensions of array of r-grid fields, write to file.
      *
      * This function opens an output file with the specified filename,
      * writes expanded fields in RField<d> real-space grid format to 
      * that file, and then closes the file. The overloaded function of 
      * the same name with a std::ostream parameter is called internally.
      *
      * \param filename  name of output file
      * \param fields  input array of RField<D> objects (r-space grid)
      * \param unitCell  original crystallographic unit cell
      * \param d  expanded dimension (greater than D)
      * \param newGridDimensions  number of grid points in added dimensions
      */
      void expandFieldsDimension(std::string filename,
                                 DArray<RField<D> > const & fields,
                                 UnitCell<D> const & unitCell,
                                 int d,
                                 DArray<int> newGridDimensions) const;

      /**
      * Write r-grid fields in a replicated unit cell to std::ostream.
      *
      * This function takes an input array of periodic fields and outputs
      * them within an expanded unit cell in which the original input unit 
      * cell has been replicated a specified number of times in each 
      * direction. Results are written to an std::ostream output stream.
      *
      * Element i of the replicas IntVec<D> parameter contains the 
      * number of unit cell replicas along direction i. 
      * 
      * \param out  output stream (i.e., output file)
      * \param fields  array of RField (r-space) fields to be replicated
      * \param unitCell  original crystallographic unit cell
      * \param replicas  number of unit cell replicas in each direction
      */
      void replicateUnitCell(std::ostream& out,
                             DArray<RField<D> > const & fields,
                             UnitCell<D> const & unitCell,
                             IntVec<D> const & replicas) const;

      /**
      * Write r-grid fields in a replicated unit cell to named file.
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
                             DArray<RField<D> > const & fields,
                             UnitCell<D> const & unitCell,
                             IntVec<D> const & replicas) const;

      ///@}

   private:

      // DFT work array for conversions basis <-> kgrid <-> r-grid
      mutable RFieldDft<D> workDft_;

      // Pointers to associated external data

      /// Pointer to spatial discretization mesh
      Mesh<D> const * meshPtr_;

      /// Pointer to FFT object
      FFT<D> const * fftPtr_;

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

      // Private accessor functions for associated external data

      /// Get spatial discretization mesh by const reference
      Mesh<D> const & mesh() const
      {
         UTIL_ASSERT(meshPtr_);
         return *meshPtr_;
      }

      /// Get FFT object by const reference
      FFT<D> const & fft() const
      {
         UTIL_ASSERT(fftPtr_);
         return *fftPtr_;
      }

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

      /// Get the associated Basis by const reference
      Basis<D> & basis() const
      {
         UTIL_ASSERT(basisPtr_);
         return *basisPtr_;
      }

      /// Get associated FileMaster by reference
      FileMaster const & fileMaster() const
      {
         UTIL_ASSERT(fileMasterPtr_);
         return *fileMasterPtr_;
      }

      /**
      * Check state of work array, allocate if necessary.
      */
      void checkWorkDft() const;

   };

   #ifndef RPC_FIELD_IO_TPP
   extern template class FieldIo<1>;
   extern template class FieldIo<2>;
   extern template class FieldIo<3>;
   #endif

} // namespace Rpc
} // namespace Pscf
#endif

#ifndef PSPC_FIELD_IO_H
#define PSPC_FIELD_IO_H

/*
* PSCF - Polymer Self-Consistent Field Theory
*
* Copyright 2016 - 2022, The Regents of the University of Minnesota
* Distributed under the terms of the GNU General Public License.
*/

#include <prdc/cpu/FFT.h>                  // member
#include <prdc/cpu/RField.h>               // function parameter
#include <prdc/cpu/RFieldDft.h>            // function parameter
#include <prdc/crystal/Basis.h>            // member
#include <prdc/crystal/SpaceGroup.h>       // member
#include <prdc/crystal/UnitCell.h>         // member

#include <pscf/mesh/Mesh.h>                // member

#include <util/misc/FileMaster.h>          // member
#include <util/containers/DArray.h>        // function parameter
#include <util/containers/Array.h>         // function parameter

namespace Pscf {
namespace Pspc
{

   using namespace Util;
   using namespace Pscf;
   using namespace Pscf::Prdc;

   /**
   * File input/output operations and format conversions for fields.
   *
   * This class provides functions to read and write arrays that contain
   * fields in any of three representations (symmetry-adapted basis, 
   * r-space grid, or Fourier k-space grid), and to convert among these 
   * representations. The functions that implement IO operations define
   * file formats for these field representations.
   *
   * \ingroup Pspc_Field_Module
   */
   template <int D>
   class FieldIo 
   {

   public:

      /**
      * Constructor.
      */
      FieldIo();

      /**
      * Destructor.
      */
      ~FieldIo();

      /**
      * Get and store addresses of associated objects.
      *
      * \param mesh  associated spatial discretization Mesh<D>
      * \param fft   associated FFT object for fast transforms
      * \param lattice  lattice system type (enumeration value)
      * \param groupName space group name string
      * \param group  associated SpaceGroup object
      * \param basis  associated Basis object
      * \param fileMaster  associated FileMaster (for file paths)
      */
      void associate(Mesh<D> const & mesh,
                     FFT<D> const & fft,
                     typename UnitCell<D>::LatticeSystem & lattice,
                     std::string & groupName,
                     SpaceGroup<D> & group,
                     Basis<D> & basis,
                     FileMaster const & fileMaster);

      /// \name Field File IO - Symmetry Adapted Basis Format
      ///@{

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
      * Read concentration or chemical potential field components from file.
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
      * Read concentration or chemical potential field components from file.
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

      ///@}
      /// \name Field File IO - Real Space Grid Format
      ///@{

      /**
      * Read single RField (field on an r-space grid) from istream.
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
      * Read array of RField objects (fields on r-space grid) from istream.
      *
      * The capacity of array fields is equal to nMonomer, and element
      * fields[i] is the RField<D> associated with monomer type i.
      * 
      * \param in  input stream (i.e., input file)
      * \param fields  array of RField fields (r-space grid)
      * \param unitCell  associated crystallographic unit cell
      */
      void readFieldsRGrid(std::istream& in, DArray< RField<D> >& fields, 
                           UnitCell<D> & unitCell) const;

      /**
      * Read array of RField objects (fields on an r-space grid) from file.
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
      * Read single RField frame from an ftmc trajectory file.
      *
      * This function opens an input trajectory file with the specified 
      * filename and reads a field frame in RField<D> real-space grid 
      * format.
      *
      * \param in  input file stream
      * \param fields  array of RField fields (r-space grid)
      * \param nMonomer  number of monomer types
      */
      void readFieldRGridData(std::istream& in, 
                              DArray< RField<D> >& fields,
                              int nMonomer) const;

      /**
      * Write a single RField (field on an r-space grid) to ostream.
      *
      * \param out  output stream 
      * \param field  field defined on r-space grid
      * \param unitCell  associated crystallographic unit cell
      * \param writeHeader  should a file header be written?
      */
      void writeFieldRGrid(std::ostream &out, 
                           RField<D> const & field, 
                           UnitCell<D> const & unitCell,
                           bool writeHeader = true) const;

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
      */
      void writeFieldRGrid(std::string filename, 
                           RField<D> const & field, 
                           UnitCell<D> const & unitCell) const;

      /**
      * Write array of RField objects (fields on r-space grid) to ostream.
      *
      * \param out  output stream (i.e., output file)
      * \param fields  array of RField fields (r-space grid)
      * \param unitCell  associated crystallographic unit cell
      * \param writeHeader  flag to write header of file if true
      */
      void writeFieldsRGrid(std::ostream& out, 
                            DArray< RField<D> > const & fields, 
                            UnitCell<D> const & unitCell,
                            bool writeHeader = true) const;

      /**
      * Write array of RField objects (fields on an r-space grid) to file.
      *
      * This function opens an output file with the specified filename, 
      * writes fields in RField<D> real-space grid format to that file, 
      * and then closes the file.
      *
      * \param filename  name of output file
      * \param fields  array of RField fields (r-space grid)
      * \param unitCell  associated crystallographic unit cell
      */
      void writeFieldsRGrid(std::string filename,
                            DArray< RField<D> > const & fields, 
                            UnitCell<D> const & unitCell) const;

      ///@}
      /// \name Field File IO - Fourier Space (K-Space) Grid Format
      ///@{

      /**
      * Read array of RFieldDft objects (k-space fields) from file.
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
      */
      void writeFieldsKGrid(std::ostream& out, 
                            DArray< RFieldDft<D> > const & fields, 
                            UnitCell<D> const & unitCell) const;
   
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
      */
      void writeFieldsKGrid(std::string filename, 
                           DArray< RFieldDft<D> > const & fields, 
                           UnitCell<D> const & unitCell) const;

      ///@}
      /// \name File IO Utilities
      ///@{

      /**
      * Reader header of field file (fortran pscf format)
      *
      * This reads the common part of the header for all field file 
      * formats. This contains the dimension of space, the unit cell, the 
      * group name and the the number of monomers. The unit cell data is
      * read into the associated UnitCell<D>, which is thus updated.
      *
      * If the associated basis is not initialized, this function will
      * attempt to initialize it using the unit cell read from file and
      * the associated group (if available) or group name. 
      * 
      * This function throws an exception if the values of "dim" read 
      * from file do not match the FieldIo template parameter D. 
      * 
      * The function does not impose any requirements on the value
      * of the input parameter nMonomer; it is overwritten and will
      * contain the value of "N_monomer" from the field file header 
      * at termination of the function.
      * 
      * If the UnitCell object passed to this function already
      * contains unit cell data, the function will check to ensure
      * that the crystal system and space group in the field file 
      * header match the data already stored. A warning will also 
      * be printed if the lattice parameters do not match. If, 
      * instead, the function is passed an empty UnitCell object, 
      * we assume that the UnitCell data is unknown or that the 
      * field file header does not need to be cross-checked with 
      * existing data. In this case, the field file header data is 
      * stored directly in the UnitCell object that was passed in.
      * 
      * \param in  input stream (i.e., input file)
      * \param nMonomer  number of fields contained in the field file
      * \param unitCell  associated crystallographic unit cell
      */
      void readFieldHeader(std::istream& in, int& nMonomer, 
                           UnitCell<D> & unitCell) const;

      /**
      * Write header for field file (fortran pscf format)
      *
      * \param out  output stream (i.e., output file)
      * \param nMonomer  number of monomer types or fields
      * \param unitCell  associated crystallographic unit cell
      */
      void writeFieldHeader(std::ostream& out, int nMonomer,
                            UnitCell<D> const & unitCell) const;

      ///@}
      /// \name Field Format Conversion
      ///@{

      /**
      * Convert a field from symmetrized basis to Fourier transform (k-grid).
      *
      * \param components coefficients of symmetry-adapted basis functions
      * \param dft discrete Fourier transform of a real field
      */
      void convertBasisToKGrid(DArray<double> const & components, 
                               RFieldDft<D>& dft) const;
   
      /**
      * Convert fields from symmetrized basis to Fourier transform (k-grid).
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
      * Convert a field from Fourier transform (kgrid) to symmetrized basis.
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
      * Convert fields from Fourier transform (k-grid) to symmetrized basis.
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
      * \param in  field defined on real-space grid
      * \param out  field in symmetry adapted basis form
      * \param checkSymmetry  boolean indicating whether to check that the 
      * symmetry of the input field matches the space group symmetry. If
      * input does not have correct symmetry, prints warning to Log::file()
      * \param epsilon  if checkSymmetry = true, epsilon is the error 
      * threshold used when comparing the k-grid and symmetry-adapted formats 
      * to determine whether field has the declared space group symmetry
      */
      void convertRGridToBasis(RField<D> const & in,
                               DArray<double> & out,
                               bool checkSymmetry = true, 
                               double epsilon = 1.0e-8) const;

      /**
      * Convert fields from spatial grid (r-grid) to symmetrized basis.
      * 
      * \param in  fields defined on real-space grid
      * \param out  fields in symmetry adapted basis form
      * \param checkSymmetry  boolean indicating whether to check that the 
      * symmetry of the input field matches the space group symmetry. If
      * input does not have correct symmetry, prints warning to Log::file()
      * \param epsilon  if checkSymmetry = true, epsilon is the error 
      * threshold used when comparing the k-grid and symmetry-adapted formats 
      * to determine whether field has the declared space group symmetry
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
      * its input, which is why argument "in" not a const reference.
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
      * its input, which is why argument "in" not a const reference.
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
      * \param in field in real space grid (r-grid) format
      * \param epsilon error threshold used when comparing the k-grid and
      * symmetry-adapted formats to determine whether field has the declared
      * space group symmetry.
      * \param verbose if field does not have symmetry and verbose = true,
      * function will write error values to Log::file().
      * \return true if the field is symmetric, false otherwise
      */
      bool hasSymmetry(RField<D> const & in, double epsilon = 1.0e-8,
                       bool verbose = true) const;

      /**
      * Check if a k-grid field has declared space group symmetry.
      *
      * \param in field in real space grid (r-grid) format
      * \param epsilon error threshold used when comparing the k-grid and
      * symmetry-adapted formats to determine whether field has the declared
      * space group symmetry.
      * \param verbose if field does not have symmetry and verbose = true,
      * function will write error values to Log::file().
      * \return true if the field is symmetric, false otherwise
      */
      bool hasSymmetry(RFieldDft<D> const & in, double epsilon = 1.0e-8,
                       bool verbose = true) const;

      ///@}

   private:

      // DFT work array for two-step conversion basis <-> kgrid <-> r-grid.
      mutable RFieldDft<D> workDft_;

      // Pointers to associated objects.

      /// Pointer to spatial discretization mesh.
      Mesh<D> const * meshPtr_;

      /// Pointer to FFT object.
      FFT<D> const * fftPtr_;

      /// Pointer to lattice system.
      typename UnitCell<D>::LatticeSystem * latticePtr_;

      /// Pointer to group name string
      std::string * groupNamePtr_;

      /// Pointer to a SpaceGroup object
      SpaceGroup<D> * groupPtr_;

      /// Pointer to a Basis object
      Basis<D> * basisPtr_;

      /// Pointer to Filemaster (holds paths to associated I/O files).
      FileMaster const * fileMasterPtr_;

      // Private accessor functions:

      /// Get spatial discretization mesh by const reference.
      Mesh<D> const & mesh() const
      {  
         UTIL_ASSERT(meshPtr_);  
         return *meshPtr_; 
      }

      /// Get FFT object by const reference.
      FFT<D> const & fft() const
      {
         UTIL_ASSERT(fftPtr_);  
         return *fftPtr_; 
      }

      /// Get group name string by const reference.
      typename UnitCell<D>::LatticeSystem & lattice() const
      {  
         UTIL_ASSERT(latticePtr_);  
         return *latticePtr_; 
      }

      /// Get group name string by const reference.
      std::string & groupName() const
      {  
         UTIL_ASSERT(groupNamePtr_);  
         return *groupNamePtr_; 
      }

      /// Get SpaceGroup by const reference.
      SpaceGroup<D> & group() const
      {
         UTIL_ASSERT(groupPtr_);  
         return *groupPtr_; 
      }

      /// Get Basis by const reference.
      Basis<D> & basis() const
      {
         UTIL_ASSERT(basisPtr_);  
         return *basisPtr_; 
      }

      /// Get FileMaster by reference.
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

   #ifndef PSPC_FIELD_IO_TPP
   extern template class FieldIo<1>;
   extern template class FieldIo<2>;
   extern template class FieldIo<3>;
   #endif

} // namespace Pspc
} // namespace Pscf
#endif

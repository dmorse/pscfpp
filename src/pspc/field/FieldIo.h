#ifndef PSPC_FIELD_IO_H
#define PSPC_FIELD_IO_H

/*
* PSCF - Polymer Self-Consistent Field Theory
*
* Copyright 2016 - 2019, The Regents of the University of Minnesota
* Distributed under the terms of the GNU General Public License.
*/

#include <pspc/field/FFT.h>                // member
#include <pspc/field/RField.h>             // function parameter
#include <pspc/field/RFieldDft.h>          // function parameter

#include <pscf/crystal/Basis.h>            // member
#include <pscf/crystal/UnitCell.h>         // member
#include <pscf/mesh/Mesh.h>                // member

#include <util/misc/FileMaster.h>          // member
#include <util/containers/DArray.h>        // function parameter
#include <util/containers/Array.h>         // function parameter

namespace Pscf {
namespace Pspc
{
   using namespace Util;
   using namespace Pscf;

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
      * \param groupName space group name string
      * \param basis  associated Basis object
      * \param fileMaster  associated FileMaster (for file paths)
      */
      void associate(Mesh<D> const & mesh,
                     FFT<D> const & fft,
                     std::string const & groupName,
                     Basis<D> const & basis,
                     FileMaster const & fileMaster);

      /// \name Field File IO - Symmetry Adapted Basis Format
      ///@{

      /**
      * Read concentration or chemical potential field components from file.
      *
      * This function reads components in a symmetry adapted basis from 
      * file in.
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
                           bool writeHeader = true) 
      const;

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
      */
      void writeFieldsRGrid(std::ostream& out, 
                            DArray< RField<D> > const & fields, 
                            UnitCell<D> const & unitCell) const;

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
      * This function throws an exception if the values of "dim" and 
      * "N_monomer" read from file do not match values of D and the 
      * input parameter nMonomer, respectively. The group name is not 
      * checked.
      * 
      * \param in  input stream (i.e., input file)
      * \param nMonomer  expected value of nMonomer
      * \param unitCell  associated crystallographic unit cell
      */
      void readFieldHeader(std::istream& in, int& nMonomer, 
                           UnitCell<D> & unitCell) const;

      /**
      * Write header for field file (fortran pscf format)
      *
      * \param out  output stream (i.e., output file)
      * \param nMonomer  number of monomer types
      * \param unitCell  associated crystallographic unit cell
      */
      void writeFieldHeader(std::ostream& out, int nMonomer,
                            UnitCell<D> const & unitCell) const;

      ///@}
      /// \name Field Format Conversion
      ///@{

      /**
      * Convert field from symmetrized basis to Fourier transform (k-grid).
      *
      * \param components coefficients of symmetry-adapted basis functions
      * \param dft discrete Fourier transform of a real field
      */
      void convertBasisToKGrid(DArray<double> const & components, 
                               RFieldDft<D>& dft) const;
   
      /**
      * Convert fields from symmetrized basis to Fourier transform (kgrid).
      * 
      * The in and out parameters are arrays of fields, in which element
      * number i is the field associated with monomer type i. 
      *
      * \param in  components of fields in symmetry adapted basis 
      * \param out fields defined as discrete Fourier transforms (k-grid)
      */
      void convertBasisToKGrid(DArray< DArray<double> > const & in,
                               DArray< RFieldDft<D> >& out) const;

      /**
      * Convert field from Fourier transform (k-grid) to symmetrized basis.
      *
      * \param in  complex DFT (k-grid) representation of a field.
      * \param out  coefficients of symmetry-adapted basis functions.
      */
      void convertKGridToBasis(RFieldDft<D> const & in, 
                               DArray<double> & out) const;

      /**
      * Convert fields from Fourier transform (kgrid) to symmetrized basis.
      * 
      * The in and out parameters are each an array of fields, in which
      * element i is the field associated with monomer type i. 
      *
      * \param in  fields defined as discrete Fourier transforms (k-grid)
      * \param out  components of fields in symmetry adapted basis 
      */
      void convertKGridToBasis(DArray< RFieldDft<D> > const & in,
                               DArray< DArray<double> > & out) const;

      /**
      * Convert fields from symmetrized basis to spatial grid (rgrid).
      * 
      * \param in  fields in symmetry adapted basis form
      * \param out fields defined on real-space grid
      */
      void convertBasisToRGrid(DArray< DArray<double> > const & in,
                               DArray< RField<D> > & out) const ;

      /**
      * Convert fields from spatial grid (rgrid) to symmetrized basis.
      * 
      * \param in  fields defined on real-space grid
      * \param out  fields in symmetry adapted basis form
      */
      void convertRGridToBasis(DArray< RField<D> > const & in,
                               DArray< DArray<double> > & out) const;

      /**
      * Convert fields from k-grid (DFT) to real space (rgrid) format.
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
      * Convert fields from spatial grid (rgrid) to k-grid format.
      * 
      * \param in  fields defined on real-space grid (r-grid)
      * \param out  fields in discrete Fourier format (k-grid)
      */
      void convertRGridToKGrid(DArray< RField<D> > const & in,
                               DArray< RFieldDft<D> > & out) const;

      ///@}
      /// \name Test Space Group Symmetry
      ///@{

      /**
      * Check if an r-grid field has the declared space group symmetry.
      *
      * \param in field in real space grid (r-grid) format
      * \return true if the field is symmetric, false otherwise
      */
      bool hasSymmetry(RField<D> const & in) const;

      /**
      * Check if a k-grid field has declared space group symmetry.
      *
      * \param in field in real space grid (r-grid) format
      * \return true if the field is symmetric, false otherwise
      */
      bool hasSymmetry(RFieldDft<D> const & in) const;

      ///@}

   private:

      // DFT work array for two-step conversion basis <-> kgrid <-> rgrid.
      mutable RFieldDft<D> workDft_;

      // Pointers to associated objects.

      /// Pointer to spatial discretization mesh.
      Mesh<D> const * meshPtr_;

      /// Pointer to FFT object.
      FFT<D> const * fftPtr_;

      /// Pointer to group name string
      std::string const * groupNamePtr_;

      /// Pointer to a Basis object
      Basis<D> const * basisPtr_;

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
      std::string const & groupName() const
      {  
         UTIL_ASSERT(groupNamePtr_);  
         return *groupNamePtr_; 
      }

      /// Get Basis by const reference.
      Basis<D> const & basis() const
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

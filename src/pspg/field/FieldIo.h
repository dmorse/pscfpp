#ifndef PSPG_FIELD_IO_H
#define PSPG_FIELD_IO_H

/*
* PSCF - Polymer Self-Consistent Field Theory
*
* Copyright 2016 - 2019, The Regents of the University of Minnesota
* Distributed under the terms of the GNU General Public License.
*/

#include <pspg/field/FFT.h>                // member
#include <pspg/field/RDField.h>             // function parameter
#include <pspg/field/RDFieldDft.h>          // function parameter

#include <pscf/crystal/Basis.h>            // member
#include <pscf/crystal/UnitCell.h>         // member
#include <pscf/mesh/Mesh.h>                // member

#include <util/misc/FileMaster.h>          // member
#include <util/containers/DArray.h>        // function parameter
#include <util/containers/Array.h>         // function parameter

namespace Pscf {
namespace Pspg
{
   using namespace Util;
   using namespace Pscf;

   /**
   * File input/output operations for fields in several file formats.
   *
   * \ingroup Pspg_Field_Module
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
      void associate(Mesh<D>& mesh,
                     FFT<D>& fft,
                     std::string& groupName,
                     Basis<D>& basis,
                     FileMaster& fileMaster);

      /// \name Field File IO
      //@{

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
      * \param in input stream (i.e., input file)
      * \param fields array of fields (symmetry adapted basis components)
      * \param unitCell  crystallographic unit cell (output)
      */
      void 
      readFieldsBasis(std::istream& in, 
                      DArray< DArray<cudaReal> >& fields,
                      UnitCell<D>& unitCell) const;

      /**
      * Read concentration or chemical potential field components from file.
      *
      * This function opens an input file with the specified filename, 
      * reads components in symmetry-adapted form from that file, and 
      * closes the file.
      *
      * \param filename name of input file
      * \param fields array of fields (symmetry adapted basis components)
      * \param unitCell  crystallographic unit cell (output)
      */
      void readFieldsBasis(std::string filename, 
                           DArray< DArray<cudaReal> >& fields,
                           UnitCell<D>& unitCell) const;

      /**
      * Write concentration or chemical potential field components to file.
      *
      * This function writes components in a symmetry adapted basis.
      *
      * \param out output stream (i.e., output file)
      * \param fields array of fields (symmetry adapted basis components)
      * \param unitCell  crystallographic unit cell 
      */
      void writeFieldsBasis(std::ostream& out, 
                            DArray< DArray<cudaReal> > const & fields,
                            UnitCell<D> const & unitCell) const;

      /**
      * Write concentration or chemical potential field components to file.
      *
      * This function opens an output file with the specified filename, 
      * writes components in symmetry-adapted form to that file, and then
      * closes the file. 
      *
      * \param filename name of input file
      * \param fields array of fields (symmetry adapted basis components)
      * \param unitCell  crystallographic unit cell 
      */
      void writeFieldsBasis(std::string filename, 
                            DArray< DArray<cudaReal> > const & fields,
                            UnitCell<D> const & unitCell) const;

      /**
      * Read array of RField objects (fields on an r-space grid) from file.
      *
      * The capacity of array fields is equal to nMonomer, and element
      * fields[i] is the RField<D> associated with monomer type i.
      * 
      * \param in input stream (i.e., input file)
      * \param fields array of RField fields (r-space grid)
      * \param unitCell  crystallographic unit cell (output)
      */
      void readFieldsRGrid(std::istream& in, 
                           DArray< RDField<D> >& fields,
                           UnitCell<D>& unitCell) const;

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
      * \param filename name of input file
      * \param fields array of RField fields (r-space grid)
      * \param unitCell  crystallographic unit cell (output)
      */
      void readFieldsRGrid(std::string filename, 
                           DArray< RDField<D> >& fields,
                           UnitCell<D>& unitCell) const;

      /**
      * Write array of RField objects (fields on an r-space grid) to file.
      *
      * \param out output stream (i.e., output file)
      * \param fields array of RField fields (r-space grid)
      * \param unitCell  crystallographic unit cell 
      */
      void writeFieldsRGrid(std::ostream& out, 
                            DArray< RDField<D> > const& fields,
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
      * \param unitCell  crystallographic unit cell 
      */
      void writeFieldsRGrid(std::string filename,
                            DArray< RDField<D> > const& fields,
                            UnitCell<D> const & unitCell) const;


      /**
      * Read a single RField objects (field on an r-space grid) from file.
      *
      * \param in  input stream
      * \param field   RField field (r-space grid)
      * \param unitCell  crystallographic unit cell 
      */
      void readFieldRGrid(std::istream& in, RDField<D> &field,
                          UnitCell<D>& unitCell) const;

      /**
      * Read a single RField objects (field on an r-space grid).
      *
      * \param filename  name of input file
      * \param field   RField field (r-space grid)
      * \param unitCell  crystallographic unit cell 
      */
      void readFieldRGrid(std::string filename, RDField<D> &field,
                          UnitCell<D>& unitCell) const;

      /**
      * Write a single RField objects (field on an r-space grid) to file.
      *
      * \param out  output stream
      * \param field   RField field (r-space grid)
      * \param unitCell  crystallographic unit cell 
      */
      void writeFieldRGrid(std::ostream& out, RDField<D> const & field,
                            UnitCell<D> const & unitCell) const;

      /**
      * Write a single RField objects (field on an r-space grid) to file.
      *
      * \param filename  name of input file
      * \param field   RField field (r-space grid)
      * \param unitCell  crystallographic unit cell 
      */
      void writeFieldRGrid(std::string filename, RDField<D> const & field,
                           UnitCell<D> const & unitCell) const;

      /**
      * Read array of RFieldDft objects (k-space fields) from file.
      *
      * The capacity of the array is equal to nMonomer, and element
      * fields[i] is the discrete Fourier transform of the field for 
      * monomer type i.
      * 
      * \param in  input stream (i.e., input file)
      * \param fields  array of RFieldDft fields (k-space grid)
      * \param unitCell  crystallographic unit cell (output)
      */
      void readFieldsKGrid(std::istream& in, 
                           DArray< RDFieldDft<D> >& fields,
                           UnitCell<D>& unitCell) const;

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
      * \param unitCell  crystallographic unit cell (output)
      */
      void readFieldsKGrid(std::string filename, 
                           DArray< RDFieldDft<D> >& fields,
                           UnitCell<D>& unitCell) const;

      /**
      * Write array of RFieldDft objects (k-space fields) to file.
      *
      * The capacity of the array fields is equal to nMonomer. Element
      * fields[i] is the discrete Fourier transform of the field for 
      * monomer type i.
      * 
      * \param out output stream (i.e., output file)
      * \param fields array of RFieldDft fields 
      * \param unitCell  crystallographic unit cell 
      */
      void writeFieldsKGrid(std::ostream& out, 
                            DArray< RDFieldDft<D> > const& fields,
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
      * \param unitCell  crystallographic unit cell 
      */
      void writeFieldsKGrid(std::string filename, 
                           DArray< RDFieldDft<D> > const& fields,
                           UnitCell<D> const & unitCell) const;

      /**
      * Write header for field file (fortran pscf format)
      *
      * \param out output stream (i.e., output file)
      * \param nMonomer number of monomer types
      * \param unitCell  crystallographic unit cell 
      */
      void writeFieldHeader(std::ostream& out, int nMonomer,
                           UnitCell<D> const & unitCell) const;

      //@}
      /// \name Field Format Conversion
      //@{

      /**
      * Convert field from symmetrized basis to Fourier transform (k-grid).
      *
      * \param components coefficients of symmetry-adapted basis functions
      * \param dft discrete Fourier transform of a real field
      */
      void convertBasisToKGrid(DArray<cudaReal> const& components, 
                               RDFieldDft<D>& dft) const;
   
      /**
      * Convert fields from symmetrized basis to Fourier transform (kgrid).
      * 
      * The in and out parameters are arrays of fields, in which element
      * number i is the field associated with monomer type i. 
      *
      * \param in  components of fields in symmetry adapted basis 
      * \param out fields defined as discrete Fourier transforms (k-grid)
      */
      void convertBasisToKGrid(DArray< DArray<cudaReal> > & in,
                               DArray< RDFieldDft<D> >& out) const;

      /**
      * Convert field from Fourier transform (k-grid) to symmetrized basis.
      *
      * \param dft complex DFT (k-grid) representation of a field.
      * \param components coefficients of symmetry-adapted basis functions.
      */
      void convertKGridToBasis(RDFieldDft<D> const& dft, 
                               DArray<cudaReal>& components) const;

      /**
      * Convert fields from Fourier transform (kgrid) to symmetrized basis.
      * 
      * The in and out parameters are each an array of fields, in which
      * element i is the field associated with monomer type i. 
      *
      * \param in  fields defined as discrete Fourier transforms (k-grid)
      * \param out  components of fields in symmetry adapted basis 
      */
      void convertKGridToBasis(DArray< RDFieldDft<D> > & in,
                               DArray< DArray<cudaReal> > & out) const;

      /**
      * Convert fields from symmetrized basis to spatial grid (rgrid).
      * 
      * \param in  fields in symmetry adapted basis form
      * \param out fields defined on real-space grid
      */
      void convertBasisToRGrid(DArray< DArray<cudaReal> >& in,
                               DArray< RDField<D> >& out) const;

      /**
      * Convert fields from spatial grid (rgrid) to symmetrized basis.
      * 
      * \param in  fields defined on real-space grid
      * \param out  fields in symmetry adapted basis form
      */
      void convertRGridToBasis(DArray< RDField<D> > & in,
                               DArray< DArray<cudaReal> > & out) const;

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
      void convertKGridToRGrid(DArray< RDFieldDft<D> > & in,
                               DArray< RDField<D> > & out) const;

      /**
      * Convert fields from spatial grid (rgrid) to k-grid format.
      * 
      * \param in  fields defined on real-space grid (r-grid)
      * \param out  fields in discrete Fourier format (k-grid)
      */
      void convertRGridToKGrid(DArray< RDField<D> > & in,
                              DArray< RDFieldDft<D> > & out) const;

      //@}

   private:

      // DFT work array for two-step conversion basis <-> kgrid <-> rgrid.
      mutable RDFieldDft<D> workDft_;

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
         // UTIL_ASSERT(meshPtr_);  
         return *meshPtr_; 
      }

      /// Get FFT object by const reference.
      FFT<D> const & fft() const
      {
         // UTIL_ASSERT(fftPtr_);  
         return *fftPtr_; 
      }

      /// Get group name string by const reference.
      std::string const & groupName() const
      {  return *groupNamePtr_; }

      /// Get Basis by const reference.
      Basis<D> const & basis() const
      {
         // UTIL_ASSERT(basisPtr_);  
         return *basisPtr_; 
      }

      /// Get FileMaster by const reference.
      FileMaster const & fileMaster() const
      {  
         UTIL_ASSERT(fileMasterPtr_);  
         return *fileMasterPtr_; 
      }

      /**
      * Reader header of field file (fortran pscf format)
      *
      * \param in input  stream (i.e., input file)
      * \param nMonomer number of monomers read from file (output)
      * \param unitCell  crystallographic unit cell (output)
      */
      void readFieldHeader(std::istream& in,
                           int& nMonomer,
                           UnitCell<D>& unitCell) const;

      /**
      * Check state of work array, allocate if necessary.
      */
      void checkWorkDft() const;

   };

   #ifndef PSPG_FIELD_IO_TPP
   extern template class FieldIo<1>;
   extern template class FieldIo<2>;
   extern template class FieldIo<3>;
   #endif

} // namespace Pspc
} // namespace Pscf
#endif

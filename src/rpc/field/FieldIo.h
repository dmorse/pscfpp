#ifndef RPC_FIELD_IO_H
#define RPC_FIELD_IO_H

/*
* PSCF - Polymer Self-Consistent Field
*
* Copyright 2015 - 2025, The Regents of the University of Minnesota
* Distributed under the terms of the GNU General Public License.
*/

#include <rp/field/FieldIoBase.h>      // base class template
#include <pscf/cpu/CppTp.h>          // base class template
#include <prdc/field/cpu/RFieldDft.h>  // base class member

namespace Pscf {
namespace Rp {

   using namespace Util;
   using namespace Prdc;
   using namespace Prdc::Cpu;

   /**
   * File input/output operations and format conversions for fields.
   *
   * The public interface of this subclass is identical to that of the 
   * base class, Rp::FieldIoBase. Member functions defined in this class 
   * all implement pure virtual functions declared by the base class for 
   * which CPU and GPU partial specializations of derived class template
   * Rp::FieldIo require different implementations. Specialization is 
   * often required because the GPU implementation requires data transfer 
   * bewteen the between CPU and GPU.
   *
   * \ingroup Rp_Field_Module
   */
   template <int D>
   class FieldIo<D, CppTp<D> > 
    : public  FieldIoBase< D, CppTp<D> >
   {

   public:

      FieldIo() = default;
      virtual ~FieldIo() = default;

      /**
      * Read array of RField objects (r-grid fields) from a stream.
      *
      * See documentation of analogous function in Rp::FieldIoBase.
      *
      * \param in  input file stream 
      * \param fields  array of RField fields (r-space grid)
      * \param unitCell  associated crystallographic unit cell
      * \return true iff the header contains a space group (isSymmetric)
      */
      bool readFieldsRGrid(std::istream& in,
                           DArray< RField<D, CppTp<D> > >& fields,
                           UnitCell<D> & unitCell) const override;

      /**
      * Read data for an array of r-grid fields, with no header section.
      *
      * See documentation of analogous function in Rp::FieldIoBase.
      *
      * \param in  input file stream
      * \param fields  array of RField fields (r-space grid)
      * \param nMonomer  number of monomer types
      */
      void readFieldsRGridData(std::istream& in,
                               DArray< RField<D, CppTp<D> > >& fields,
                               int nMonomer) const override;

      /**
      * Read a single RField (field on an r-space grid) from a stream.
      *
      * See documentation of analogous function in Rp::FieldIoBase.
      *
      * \param in  input file stream 
      * \param field  fields defined on r-space grid
      * \param unitCell  associated crystallographic unit cell
      * \return true iff the header contains a space group (isSymmetric)
      */
      bool readFieldRGrid(std::istream &in,
                          RField<D, CppTp<D> > & field,
                          UnitCell<D>& unitCell) const override;

      /**
      * Write array of RField objects (fields on r-space grid) to a stream.
      *
      * See documentation of analogous function in Rp::FieldIoBase.
      *
      * \param out  output stream (i.e., output file)
      * \param fields  array of RField objects (fields on r-space grid)
      * \param unitCell  associated crystallographic unit cell
      * \param writeHeader  flag to write file header if true
      * \param isSymmetric  Do fields have a space group symmetry ?
      * \param writeMeshSize  Should mesh size be written in header?
      */
      void writeFieldsRGrid(std::ostream& out,
                            DArray< RField<D, CppTp<D> > > const & fields,
                            UnitCell<D> const & unitCell,
                            bool writeHeader = true,
                            bool isSymmetric = true,
                            bool writeMeshSize = true) const override;

      /**
      * Write a single RField (field on an r-space grid) to a stream.
      *
      * See documentation of analogous function in Rp::FieldIoBase.
      *
      * \param out  output stream
      * \param field  field defined on r-space grid
      * \param unitCell  associated crystallographic unit cell
      * \param writeHeader  Should a file header be written?
      * \param isSymmetric  Does the field have a space group symmetry?
      */
      void writeFieldRGrid(std::ostream &out,
                           RField<D, CppTp<D> > const & field,
                           UnitCell<D> const & unitCell,
                           bool writeHeader = true,
                           bool isSymmetric = true) 
      const override;

      /**
      * Read array of RFieldDft objects (k-space fields) from a stream.
      *
      * See documentation of analogous function in Rp::FieldIoBase.
      *
      * \param in  input stream (i.e., input file)
      * \param fields  array of RFieldDft fields (k-space grid)
      * \param unitCell  associated crystallographic unit cell
      */
      void readFieldsKGrid(std::istream& in,
                           DArray< RFieldDft<D> >& fields,
                           UnitCell<D> & unitCell) 
      const override;

      /**
      * Write array of RFieldDft objects (k-space fields) to file.
      *
      * See documentation of analogous function in Rp::FieldIoBase.
      *
      * \param out  output stream (i.e., output file)
      * \param fields  array of RFieldDft fields
      * \param unitCell  associated crystallographic unit cell
      * \param isSymmetric  Does this field have space group symmetry?
      */
      void writeFieldsKGrid(std::ostream& out,
                            DArray< RFieldDft<D> > const & fields,
                            UnitCell<D> const & unitCell,
                            bool isSymmetric = true) 
      const override;

      /**
      * Convert a field from symmetrized basis to Fourier grid (k-grid).
      *
      * See documentation of analogous function in Rp::FieldIoBase.
      *
      * \param components  coefficients of in symmetry-adapted basis
      * \param dft  discrete Fourier transform of a real field
      */
      void convertBasisToKGrid(DArray<double> const & components,
                               RFieldDft<D>& dft) 
      const override;

      /**
      * Convert a field from Fourier (k-grid) to symmetrized basis form.
      *
      * See documentation of analogous function in Rp::FieldIoBase.
      *
      * \param in  discrete Fourier transform (k-grid) of a field
      * \param out  components of field in asymmetry-adapted Fourier basis
      * \param checkSymmetry  flag indicating whether to check symmetry
      * \param epsilon  error tolerance for symmetry test (if any)
      */
      void convertKGridToBasis(RFieldDft<D> const & in,
                               DArray<double> & out,
                               bool checkSymmetry = true,
                               double epsilon = 1.0e-8) 
      const override;

      /**
      * Check if a k-grid field has the declared space group symmetry.
      *
      * See documentation of analogous function in Rp::FieldIoBase.
      *
      * \param in field in real space grid (r-grid) format
      * \param epsilon error threshold used to test for symmetry
      * \param verbose  if true, write error to Log::file()
      * \return true iff the field is symmetric to within tolerance
      */
      bool hasSymmetry(RFieldDft<D> const & in, 
                       double epsilon = 1.0e-8,
                       bool verbose = true) 
      const override;

      /**
      * Compare two fields in r-grid format, output a report.
      *
      * Outputs maximum and root-mean-squared differences to the standard
      * Log file.
      *
      * \param field1  first array of fields (r-grid format)
      * \param field2  second array of fields (r-grid format)
      */
      void compareFieldsRGrid(DArray< RField<D, CppTp<D> > > const & field1,
                              DArray< RField<D, CppTp<D> > > const & field2) 
      const override;

      /**
      * Rescale a single r-grid field by a scalar factor.
      *
      * See documentation of analogous function in Rp::FieldIoBase.
      * Multiplication is done in-place, and so modifies the input.
      *
      * \param field  real space (r-grid) field (in-out)
      * \param factor  real scalar by which to multiply all elements
      */
      void scaleFieldRGrid(RField<D, CppTp<D> >& field, double factor) 
      const override;

      /**
      * Expand spatial dimension of an array of r-grid fields.
      *
      * See documentation of analogous function in Rp::FieldIoBase.
      *
      * \param out  output file stream 
      * \param fields  input array of D-dimensional r-grid fields
      * \param unitCell  original D-dimensional unit cell
      * \param d  expanded spatial dimension (d > D)
      * \param newGridDimensions  number of grid points in added dimensions
      */
      void expandRGridDimension(
                          std::ostream &out,
                          DArray< RField<D, CppTp<D> > > const & fields,
                          UnitCell<D> const & unitCell,
                          int d,
                          DArray<int> const& newGridDimensions) 
      const override;

      /**
      * Write r-grid fields in a replicated unit cell to std::ostream.  
      *
      * See documentation of analogous function in Rp::FieldIoBase.
      *
      * \param out  output file stream 
      * \param fields  array of RField (r-space) fields to be replicated
      * \param unitCell  original crystallographic unit cell
      * \param replicas  number of unit cell replicas in each direction
      */ 
      void replicateUnitCell(
                          std::ostream& out,
                          DArray< RField<D, CppTp<D> > > const & fields,
                          UnitCell<D> const & unitCell,
                          IntVec<D> const & replicas) 
      const override;

      /// Alias for base class
      using FieldIoBase =  Rp::FieldIoBase<D, CppTp<D> >;

      // Inherited public member functions (to avoid hiding)
      using FieldIoBase::writeFieldBasis;
      using FieldIoBase::writeFieldsBasis;
      using FieldIoBase::readFieldsRGrid;
      using FieldIoBase::readFieldsRGridData;
      using FieldIoBase::readFieldRGrid;
      using FieldIoBase::writeFieldsRGrid;
      using FieldIoBase::writeFieldRGrid;
      using FieldIoBase::readFieldsKGrid;
      using FieldIoBase::writeFieldsKGrid;
      using FieldIoBase::replicateUnitCell;
      using FieldIoBase::expandRGridDimension;
      using FieldIoBase::convertBasisToKGrid;
      using FieldIoBase::convertKGridToBasis;
      using FieldIoBase::hasSymmetry;
      using FieldIoBase::compareFieldsRGrid;
      using FieldIoBase::scaleFieldsRGrid;

   };

} // namespace Rp
} // namespace Pscf 

// Explicit instantiation declarations
namespace Pscf {
   namespace Rp {
      extern template class FieldIoBase<1, CppTp<1> >;
      extern template class FieldIoBase<2, CppTp<2> >;
      extern template class FieldIoBase<3, CppTp<3> >;
      extern template class FieldIo<1, CppTp<1> >;
      extern template class FieldIo<2, CppTp<2> >;
      extern template class FieldIo<3, CppTp<3> >;
   }
} 
#endif

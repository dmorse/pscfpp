#ifndef RPG_FIELD_IO_H
#define RPG_FIELD_IO_H

/*
* PSCF - Polymer Self-Consistent Field Theory
*
* Copyright 2015 - 2025, The Regents of the University of Minnesota
* Distributed under the terms of the GNU General Public License.
*/

#include <prdc/field/FieldIoReal.h>     // base class template
#include <prdc/cuda/RField.h>           // template parameter
#include <prdc/cuda/RFieldDft.h>        // template parameter
#include <prdc/cuda/FFT.h>              // template parameter

// Forward declarations for classes used only via references or pointers
namespace Util {
   class FileMaster;
   template <typename T> class DArray;
}
namespace Pscf {
   namespace Prdc {
      template <int D> class UnitCell;
   }
}

namespace Pscf {
namespace Rpg {

   using namespace Util;
   using namespace Pscf::Prdc;
   using namespace Pscf::Prdc::Cuda;

   /**
   * File input/output operations and format conversions for fields.
   *
   * Please refer to the documentation of the base class Prdc::FieldIoReal 
   * for more complete API documentation for this class template, for 
   * reasons discussed below.
   *
   * Class template Rpg::FieldIo<int D> is derived from a partial 
   * specialization of template Prdc::FieldIoReal<D, RFRT, RFKT, FFFT> 
   * that is implemented using classes RFRT=RField<D>, RFKT=RFieldDft<D>, 
   * and FFFT=FFT<D> that are all defined in the Prdc::Cuda subspace, and 
   * that can use GPU hardware. Rpg::FieldIo is thus a specialization of
   * the FieldIoReal template with GPU acceleration. An analogous class
   * named Rpc::FieldIo that is designed for standard CPU hardware is
   * defined in the Pscf::Rpc namespace 
   *
   * The public interface of Rpg::FieldIo is identical to that of the
   * base class template Prdc::FieldIoReal. All member functions defined 
   * in this Rpg::FieldIo are reimplemented virtual functions that are 
   * declared and documented in Prdc::FieldIoReal, but that have trivial
   * do-nothing implementations in the base class. These are all functions 
   * for which different implementations are required for the CPU and GPU 
   * variants, and which are thus also reimplemented in Rpg::FieldIo.
   *
   * \ingroup Rpg_Field_Module
   */
   template <int D>
   class FieldIo 
     : public  FieldIoReal< D, RField<D>, RFieldDft<D>, FFT<D> >
   {

   public:

      typedef FieldIoReal<D, RField<D>, RFieldDft<D>, FFT<D> > Base;

      // Inherited public member functions
      using Base::associate;
      using Base::setFileMaster;
      using Base::readFieldsBasis;
      using Base::readFieldBasis;
      using Base::writeFieldBasis;
      using Base::writeFieldsBasis;
      using Base::readFieldsRGrid;
      using Base::readFieldsRGridData;
      using Base::readFieldRGrid;
      using Base::writeFieldsRGrid;
      using Base::writeFieldRGrid;
      using Base::readFieldsKGrid;
      using Base::writeFieldsKGrid;
      using Base::convertBasisToKGrid;
      using Base::convertKGridToBasis;
      using Base::convertBasisToRGrid;
      using Base::convertRGridToBasis;
      using Base::convertKGridToRGrid;
      using Base::convertRGridToKGrid;
      using Base::hasSymmetry;
      using Base::scaleFieldsBasis;
      using Base::scaleFieldsRGrid;
      using Base::replicateUnitCell;
      using Base::expandRGridDimension;
      using Base::readFieldHeader;
      using Base::writeFieldHeader;

      /**
      * Default constructor.
      */
      FieldIo();

      /**
      * Destructor.
      */
      ~FieldIo();

      /**
      * Read array of RField objects (r-grid fields) from a stream.
      *
      * See documentation of analogous function in Prdc::FieldIoReal.
      *
      * \param in  input file stream 
      * \param fields  array of RField fields (r-space grid)
      * \param unitCell  associated crystallographic unit cell
      * \return true iff the header contains a space group (isSymmetric)
      */
      bool readFieldsRGrid(std::istream& in,
                           DArray< RField<D> >& fields,
                           UnitCell<D> & unitCell) const;

      /**
      * Read data for an array of r-grid fields, with no header section.
      *
      * See documentation of analogous function in Prdc::FieldIoReal.
      *
      * \param in  input file stream
      * \param fields  array of RField fields (r-space grid)
      * \param nMonomer  number of monomer types
      */
      void readFieldsRGridData(std::istream& in,
                               DArray< RField<D> >& fields,
                               int nMonomer) const;

      /**
      * Read a single RField (field on an r-space grid) from a stream.
      *
      * See documentation of analogous function in Prdc::FieldIoReal.
      *
      * \param in  input file stream 
      * \param field  fields defined on r-space grid
      * \param unitCell  associated crystallographic unit cell
      * \return true iff the header contains a space group (isSymmetric)
      */
      bool readFieldRGrid(std::istream &in,
                          RField<D> & field,
                          UnitCell<D>& unitCell) const;

      /**
      * Write array of RField objects (fields on r-space grid) to a stream.
      *
      * See documentation of analogous function in Prdc::FieldIoReal.
      *
      * \param out  output stream (i.e., output file)
      * \param fields  array of RField objects (fields on r-space grid)
      * \param unitCell  associated crystallographic unit cell
      * \param writeHeader  flag to write file header if true
      * \param isSymmetric  Do fields have a space group symmetry ?
      * \param writeMeshSize  Should mesh size be written in header?
      */
      void writeFieldsRGrid(std::ostream& out,
                            DArray< RField<D> > const & fields,
                            UnitCell<D> const & unitCell,
                            bool writeHeader = true,
                            bool isSymmetric = true,
                            bool writeMeshSize = true) const;

      /**
      * Write a single RField (field on an r-space grid) to a stream.
      *
      * See documentation of analogous function in Prdc::FieldIoReal.
      *
      * \param out  output stream
      * \param field  field defined on r-space grid
      * \param unitCell  associated crystallographic unit cell
      * \param writeHeader  Should a file header be written?
      * \param isSymmetric  Does the field have a space group symmetry?
      */
      void writeFieldRGrid(std::ostream &out,
                           RField<D> const & field,
                           UnitCell<D> const & unitCell,
                           bool writeHeader = true,
                           bool isSymmetric = true) const;

      /**
      * Read array of RFieldDft objects (k-space fields) from a stream.
      *
      * See documentation of analogous function in Prdc::FieldIoReal.
      *
      * \param in  input stream (i.e., input file)
      * \param fields  array of RFieldDft fields (k-space grid)
      * \param unitCell  associated crystallographic unit cell
      */
      void readFieldsKGrid(std::istream& in,
                           DArray< RFieldDft<D> >& fields,
                           UnitCell<D> & unitCell) const;

      /**
      * Write array of RFieldDft objects (k-space fields) to file.
      *
      * See documentation of analogous function in Prdc::FieldIoReal.
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
      * Convert a field from symmetrized basis to Fourier grid (k-grid).
      *
      * See documentation of analogous function in Prdc::FieldIoReal.
      *
      * \param components  coefficients of in symmetry-adapted basis
      * \param dft  discrete Fourier transform of a real field
      */
      void convertBasisToKGrid(DArray<double> const & components,
                               RFieldDft<D>& dft) const;

      /**
      * Convert a field from Fourier (k-grid) to symmetrized basis form.
      *
      * See documentation of analogous function in Prdc::FieldIoReal.
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
      * Check if a k-grid field has the declared space group symmetry.
      *
      * See documentation of analogous function in Prdc::FieldIoReal.
      *
      * \param in field in real space grid (r-grid) format
      * \param epsilon error threshold used to test for symmetry
      * \param verbose  if true, write error to Log::file()
      * \return true iff the field is symmetric to within tolerance
      */
      bool hasSymmetry(RFieldDft<D> const & in, 
                       double epsilon = 1.0e-8,
                       bool verbose = true) const;

      /**
      * Rescale a single field in basis format by a scalar factor.
      *
      * See documentation of analogous function in Prdc::FieldIoReal.
      * Multiplication is done in-place, and so modifies the input.
      *
      * \param field  field in basis format (in-out)
      * \param factor  real scalar by which to multiply all components
      */
      void scaleFieldBasis(DArray<double>& field, double factor) const;

      /**
      * Rescale a single r-grid field by a scalar factor.
      *
      * See documentation of analogous function in Prdc::FieldIoReal.
      * Multiplication is done in-place, and so modifies the input.
      *
      * \param field  real space (r-grid) field (in-out)
      * \param factor  real scalar by which to multiply all elements
      */
      void scaleFieldRGrid(RField<D>& field, double factor) const;
      
      /**
      * Expand spatial dimension of an array of r-grid fields.
      *
      * See documentation of analogous function in Prdc::FieldIoReal.
      *
      * \param out  output file stream 
      * \param fields  input array of D-dimensional r-grid fields
      * \param unitCell  original D-dimensional unit cell
      * \param d  expanded spatial dimension (d > D)
      * \param newGridDimensions  number of grid points in added dimensions
      */
      void expandRGridDimension(
                          std::ostream &out,
                          DArray<RField<D> > const & fields,
                          UnitCell<D> const & unitCell,
                          int d,
                          DArray<int> const& newGridDimensions) const;

      /**
      * Write r-grid fields in a replicated unit cell to std::ostream.  
      *
      * See documentation of analogous function in Prdc::FieldIoReal.
      *
      * \param out  output file stream 
      * \param fields  array of RField (r-space) fields to be replicated
      * \param unitCell  original crystallographic unit cell
      * \param replicas  number of unit cell replicas in each direction
      */ 
      void replicateUnitCell(
                          std::ostream& out,
                          DArray< RField<D> > const & fields,
                          UnitCell<D> const & unitCell,
                          IntVec<D> const & replicas) const;

   protected:

      using Base::mesh;
      using Base::fft;
      using Base::lattice;
      using Base::hasGroup;
      using Base::groupName;
      using Base::group;
      using Base::basis;
      using Base::fileMaster;

   };

   #ifndef RPG_FIELD_IO_TPP
   extern template class FieldIo<1>;
   extern template class FieldIo<2>;
   extern template class FieldIo<3>;
   #endif

} // namespace Rpg

#ifndef RPG_FIELD_IO_TPP
namespace Prdc {
   using namespace Pscf::Prdc::Cuda;
   extern template class FieldIoReal<1, RField<1>, RFieldDft<1>, FFT<1>>;
   extern template class FieldIoReal<2, RField<2>, RFieldDft<2>, FFT<2>>;
   extern template class FieldIoReal<3, RField<3>, RFieldDft<3>, FFT<3>>;
} 
#endif

} // namespace Pscf
#endif

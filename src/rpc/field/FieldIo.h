#ifndef RPC_FIELD_IO_H
#define RPC_FIELD_IO_H

/*
* PSCF - Polymer Self-Consistent Field Theory
*
* Copyright 2016 - 2022, The Regents of the University of Minnesota
* Distributed under the terms of the GNU General Public License.
*/

#include <prdc/field/FieldIoReal.h>     // base class template

#include <prdc/cpu/RField.h>            // template parameter
#include <prdc/cpu/RFieldDft.h>         // template parameter
#include <prdc/cpu/FFT.h>               // template parameter

// Forward declarations for classes used only via references or pointers
namespace Util {
   class FileMaster;
   template <typename T> class DArray;
}
namespace Pscf {
   template <int D> class Mesh;
   namespace Prdc {
      template <int D> class Basis;
      template <int D> class SpaceGroup;
      template <int D> class UnitCell;
      namespace Cpu {
         template <int D> class FFT;
      }
   }
}

namespace Pscf {
namespace Rpc {

   using namespace Util;
   using namespace Pscf::Prdc;
   using namespace Pscf::Prdc::Cpu;

   /**
   * File input/output operations and format conversions for fields.
   *
   * Please refer to the documentation of the base class Prdc::FieldIoReal 
   * for more complete API documentation for this class template, for
   * reasons discussed below.
   *
   * Class template Rpc::FieldIo<int D> is derived from a partial 
   * specialization of template Prdc::FieldIoReal<D, RFRT, RFKT, FFFT> 
   * that is implemented using classes RField<D>, RFieldDft<D>, and 
   * FFT<D> that are all defined in the Prdc::Cpu subspace, and that 
   * use only conventional use CPU hardware. Rpc::FieldIo is thus a 
   * specialization of FieldIoReal for CPU hardware. An analogous
   * class template named named Rpg::FieldIo that is designed to use 
   * a GPU is defined in the Pscf::Rpg namespace 
   *
   * The pubiic interface of Rpc::FieldIo is identical to that of the
   * base class template Prdc::FieldIoReal. All member functions defined 
   * in this class template are virtual functions that are defined and
   * documented in Prdc::FieldIoReal, but provided nontrivial 
   * implementations here. These are all functions for which different 
   * implementations are required for the CPU and GPU variants, and
   * which are thus also reimplemented in Rpg::FieldIo.
   *
   * \ingroup Rpc_Field_Module
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
      using Base::replicateUnitCell;
      using Base::expandRGridDimension;
      using Base::readFieldHeader;
      using Base::writeFieldHeader;

      /**
      * Constructor.
      */
      FieldIo();

      /**
      * Destructor.
      */
      ~FieldIo();

      /**
      * Read array of RField objects (r-grid fields) from a stream.
      *
      * See Pscf::Prdc::FieldIoReal::readFieldsRGrid .
      *
      * \param in  input stream (i.e., input file)
      * \param fields  array of RField fields (r-space grid)
      * \param unitCell  associated crystallographic unit cell
      */
      void readFieldsRGrid(std::istream& in,
                           DArray< RField<D> >& fields,
                           UnitCell<D> & unitCell) const;

      /**
      * Read data for an array of r-grid fields, with no header section.
      *
      * See Pscf::Prdc::FieldIoReal::readFieldsRGridData .
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
      * \param in  input stream (i.e., input file)
      * \param field  fields defined on r-space grid
      * \param unitCell  associated crystallographic unit cell
      */
      void readFieldRGrid(std::istream &in,
                           RField<D> & field,
                           UnitCell<D>& unitCell) const;

      /**
      * Write array of RField objects (fields on r-space grid) to a stream.
      *
      * \param out  output stream (i.e., output file)
      * \param fields  array of RField objects (fields on r-space grid)
      * \param unitCell  associated crystallographic unit cell
      * \param writeHeader  flag to write header of file if true
      * \param isSymmetric  Do fields have a space group symmetry ?
      * \param writeMeshSize Should mesh size be written in header?
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
      * Read array of RFieldDft objects (k-space fields) from a stream.
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
      * \param components  coefficients of in symmetry-adapted basis
      * \param dft  discrete Fourier transform of a real field
      */
      void convertBasisToKGrid(DArray<double> const & components,
                               RFieldDft<D>& dft) const;

      /**
      * Convert a field from Fourier (k-grid) to symmetrized basis form.
      *
      * If the checkSymmetry parameter is true, this function checks if
      * the input field satisfies the space group symmetry to within a
      * tolerance given by the epsilon parameter, and prints a warning 
      * to Log::file() if it does not.
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
      bool hasSymmetry(RFieldDft<D> const & in, 
                       double epsilon = 1.0e-8,
                       bool verbose = true) const;

      /**
      * Expand dimension of an array of r-grid fields, write to ostream.
      *
      * This function is used for template dimension D < 3, and allows a
      * 1D or 2D field to be expanded into a higher dimensional 2D or 3D
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
      void expandRGridDimension(
                          std::ostream &out,
                          DArray<RField<D> > const & fields,
                          UnitCell<D> const & unitCell,
                          int d,
                          DArray<int> const& newGridDimensions) const;

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

   #ifndef RPC_FIELD_IO_TPP
   extern template class FieldIo<1>;
   extern template class FieldIo<2>;
   extern template class FieldIo<3>;
   #endif

} // namespace Rpc

#ifndef RPC_FIELD_IO_TPP
namespace Prdc {
   using namespace Cpu;
   extern template class FieldIoReal<1, RField<1>, RFieldDft<1>, FFT<1>>;
   extern template class FieldIoReal<2, RField<2>, RFieldDft<2>, FFT<2>>;
   extern template class FieldIoReal<3, RField<3>, RFieldDft<3>, FFT<3>>;
} 
#endif

} // namespace Pscf
#endif

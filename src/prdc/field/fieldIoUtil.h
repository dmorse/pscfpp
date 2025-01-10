#ifndef PRDC_FIELD_IO_UTIL_H
#define PRDC_FIELD_IO_UTIL_H

/*
* PSCF - Polymer Self-Consistent Field Theory
*
* Copyright 2016 - 2022, The Regents of the University of Minnesota
* Distributed under the terms of the GNU General Public License.
*/

#include <util/containers/DArray.h>   
#include <pscf/math/IntVec.h>   

// Forward declarations for classes used only via references or pointers
namespace Pscf {
   template <int D> class Mesh;
   namespace Prdc {
      template <int D> class UnitCell;
      template <int D> class Basis;
   }
}

namespace Pscf {
namespace Prdc {

   using namespace Util;
   using namespace Pscf;

   // Utilities for checking field and array dimensions

   /**
   * Check allocation of a single field, allocate if necessary.
   *
   * Template parameter FT is a field type, such as RField<D> or 
   * RFieldDft<D>, that has an allocate statement that takes an
   * IntVec<D> of mesh dimensions.
   *
   * On successful exit, the mesh dimensions for the field are
   * equal to those given in parameter meshDimensions.
   *
   * If the field is allocated on entry, the above condition is
   * checked, and an Exception is thrown if is violated.  If the
   * fields is not allocated on entry, it is allocated with the
   * correct dimensions.
   *
   * \param field  field object of type FT
   * \param dimensions  required mesh dimensions
   */
   template <int D, class FT>
   void checkAllocateField(FT& field, IntVec<D> const& dimensions); 

   /**
   * Check allocation of an array of fields, allocate if necessary.
   *
   * Template parameter FT is a field type, such as RField<D> or 
   * RFieldDft<D>, that has an allocate function that takes an 
   * IntVec<D> of mesh dimensions.
   *
   * On successful exit, the capacity of the DArray fields is equal
   * to nMonomer, and the mesh dimensions for each field are equal to
   * those given in parameter meshDimensions.
   *
   * If the fields array is allocated on entry, the above conditions
   * are checked, and an Exception is thrown if any are not met. 
   * If the arrray is not allocated on entry, it is allocated with 
   * the required dimensions.
   *
   * \param fields  DArray of fields of type FT
   * \param dimensions  required mesh dimensions
   * \param nMomoner  number of monomer types (in)
   */
   template <int D, class FT>
   void checkAllocateFields(DArray< FT >& fields, 
                            IntVec<D> const& dimensions,
                            int nMonomer);

   /**
   * Check allocation of a DArray of 1D arrays, allocate if necessary.
   *
   * Template parameter AT is an allocatable array type, such as 
   * Util::DArray<double> or Pscf::HostDArray<double>, that has an allocate 
   * function that takes an integer capacity as its only parameter.
   *
   * On successful exit, the capacity of the container arrays is equal to 
   * the parameter nMonomer, and the capacity for each element of type AT 
   * is equal to the parameter capacity.
   *
   * If the container arrays is allocated on entry, the above conditions are
   * checked, and an Exception is thrown if any are violated.  If arrays is
   * not allocated on entry, it is allocated with the required dimensions.
   *
   * \param arrays  DArray of array of type AT (in/out)
   * \param dimensions  required mesh dimensions (in)
   * \param nMomoner  number of monomer types (in)
   */
   template <int D, class AT>
   void checkAllocateArrays(DArray< AT >& arrays, 
                            IntVec<D> const& dimensions,
                            int nMonomer);


   // Templates for RGrid data IO

   /**
   * Read data for array of r-grid fields, with no header section.
   *
   * This function reads the data section of the rgrid-field format
   * for multiple monomer types, with no header. 
   *
   * The template parameter AT must be an array type that provides
   * an overloaded [] subscript operator.
   *
   * \param in  input file stream
   * \param fields  array of r-grid fields (r-space grid) (out)
   * \param nMonomer  number of monomer types (in)
   */
   template <int D, class AT>
   void readRGridData(std::istream& in,
                      DArray<AT>& fields,
                      IntVec<D> const& dimensions,
                      int nMonomer);
 
   /**
   * Read data for a single r-grid field, with no header section.
   *
   * This function reads the data section of an rgrid-field format
   * for a single monomer type, with no header. 
   *
   * The template parameter AT must be an array type that provides
   * an overloaded [] subscript operator.
   *
   * \param in  input file stream
   * \param field  array containing a single r-grid field (out)
   */
   template <int D, class AT>
   void readRGridData(std::istream& in, 
                      AT& field,
                      IntVec<D> const& dimensions);
 
   /**
   * Write data for array of r-grid fields, with no header section.
   *
   * This function writes the data section of the rgrid-field format
   * for a multiple monomer types, with no header section.
   *
   * The template parameter AT must be an array type that provides 
   * an overloaded [] subscript operator.
   *
   * \param out  output file stream
   * \param fields  array of r-grid fields (r-space grid) (in)
   */
   template <int D, class AT>
   void writeRGridData(std::ostream& out, 
                       DArray<AT> const& fields,
                       IntVec<D> const& dimensions,
                       int nMonomer);
 
   /**
   * Write data for a single r-grid field, with no header section.
   *
   * This function writes the data section of an rgrid-field format
   * for a single monomer type, with no header. 
   *
   * The template parameter AT must be an array type that provides
   * an overloaded [] subscript operator.
   *
   * \param out  output file stream
   * \param field  array containing a single r-grid field (out)
   */
   template <int D, class AT>
   void writeRGridData(std::ostream& in, 
                       AT const& field,
                       IntVec<D> const& dimensions);

   // Templates for KGrid data IO
 
   /**
   * Read data for array of k-grid fields, with no header section.
   *
   * This function reads the data section of the kgrid-field format
   * for multiple monomer types, with no header. 
   *
   * The template parameter AT must be an array type that provides
   * an overloaded [] subscript operator.
   *
   * \param in  input file stream
   * \param fields  array of k-grid fields (r-space grid) (out)
   * \param dftDimensions dimensions of real DFT mesh (in)
   * \param nMonomer  number of monomer types (in)
   */
   template <int D, class AT>
   void readKGridData(std::istream& in,
                      DArray<AT>& fields,
                      IntVec<D> const& dftDimensions,
                      int nMonomer);
 
   /**
   * Read data for a single k-grid field, with no header section.
   *
   * This function reads the data section of an kgrid-field format
   * for a single monomer type, with no header. 
   *
   * The template parameter AT must be an array type that provides
   * an overloaded [] subscript operator.
   *
   * \param in  input file stream
   * \param field  array containing a single k-grid field (out)
   * \param dftDimensions dimensions of real DFT mesh (in)
   */
   template <int D, class AT>
   void readKGridData(std::istream& in, 
                      AT& field,
                      IntVec<D> const& dftDimensions);
 
   /**
   * Write data for array of k-grid fields, with no header section.
   *
   * This function writes the data section of the kgrid-field format
   * for a multiple monomer types, with no header section.
   *
   * The template parameter AT must be an array type that provides 
   * an overloaded [] subscript operator.
   *
   * \param out  output file stream
   * \param fields  array of k-grid fields (r-space grid) (in)
   * \param dftDimensions dimensions of real DFT mesh (in)
   * \param nMonomer  number of monomer types (in)
   */
   template <int D, class AT>
   void writeKGridData(std::ostream& out, 
                       DArray<AT> const& fields,
                       IntVec<D> const& dftDimensions,
                       int nMonomer);
 
   /**
   * Write data for a single k-grid field, with no header section.
   *
   * This function writes the data section of an kgrid-field format
   * for a single monomer type, with no header. 
   *
   * The template parameter AT must be an array type that provides
   * an overloaded [] subscript operator.
   *
   * \param out  output file stream
   * \param field  array containing a single k-grid field (in)
   * \param dftDimensions dimensions of real DFT mesh (in)
   */
   template <int D, class AT>
   void writeKGridData(std::ostream& in, 
                       AT const& field,
                       IntVec<D> const& dftDimensions);

   // Templates for Io of symmetrized basis format

   /**
   * Read a set of fields in basis format.
   *
   * The field header should be read before calling this function to 
   * obtain the number of basis functions in the input file, which is
   * passed to this function as parameter nStarIn.
   *
   * On input, the capacity of the fields array must equal the number
   * of monomer types represented by data columns in the input file.
   * The DArray<double> arrays associated with different monomer types
   * must all of the same nonzero capacity, denoted by fieldCapacity.
   * The number of basis functions that are processes is the smallar 
   * of nStarIn (the number available in the file) and fieldCapacity
   * (the number for which there is room in the arrays).
   * 
   * \param in  input file stream
   * \param fields  array of field components
   * \param unitCell associated crystallographic unit cell
   * \param mesh associated computational mesh
   * \param basis associated symmetry adapted basis
   * \param nStarIn number of stars declared in headers
   */
   template <int D>
   void readBasisData(std::istream& in,
                      DArray< DArray<double> >& fields,
                      UnitCell<D> const& unitCell,
                      Mesh<D> const& mesh,
                      Basis<D> const& basis,
                      int nStarIn);

   /**
   * Write array of fields in basis format, without a header.
   *
   * The number of monomer types is given by the capacity of the fields
   * array. On entry, the DArray<double> arrays associated with different 
   * monomer types must all of the same nonzero capacity, denoted by 
   * fieldCapacity. The number of basis functions written is the 
   * lesser fieldCapacity and the number of uncancelled basis functions
   * in the basis. 
   *
   * \param in  input file stream
   * \param fields  array of field components
   * \param basis associated symmetry adapted basis
   */
   template <int D>
   void writeBasisData(std::ostream &out,
                       DArray<DArray<double> > const & fields,
                       Basis<D> const & basis);

   /**
   * Convert a real field from symmetrized basis to Fourier grid.
   *
   * \param components coefficients of symmetry-adapted basis functions
   * \param dft discrete Fourier transform of a real field
   */
   template <int D, class AT>
   void convertBasisToKGrid(DArray<double> const & components,
                            AT& dft,
                            Basis<D> const& basis,
                            IntVec<D> const& dftDimensions);

   /**
   * Convert a real field from Fourier grid to symmetrized basis.
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
   template <int D, class AT>
   void convertKGridToBasis(AT const & in,
                            DArray<double> & out,
                            Basis<D> const& basis,
                            IntVec<D> const& dftDimensions,
                            bool checkSymmetry = true,
                            double epsilon = 1.0e-8);

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
   template <int D, class AT>
   bool hasSymmetry(AT const& in, 
                    Basis<D> const& basis,
                    IntVec<D> const& dftDimensions,
                    double epsilon = 1.0e-8,
                    bool verbose = true);


} // namespace Prdc
} // namespace Pscf
#include "fieldIoUtil.tpp"
#endif

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

#include <string>

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
   * RFieldDft<D>, that has an allocate function that takes an 
   * IntVec<D> of mesh dimensions, and a dimensions function that
   * returns an IntVec<D> of mesh dimensions.
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
   * IntVec<D> of mesh dimensions, and a dimensions function that
   * returns an IntVec<D> of mesh dimensions.
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
   * \param fields  DArray of fields of type FT (in/out)
   * \param dimensions  required mesh dimensions (in)
   * \param nMomoner  required number of monomer types (in)
   */
   template <int D, class FT>
   void checkAllocateFields(DArray< FT >& fields, 
                            IntVec<D> const& dimensions,
                            int nMonomer);

   /**
   * Inspect dimensions of a DArray of fields, each of type FT.
   *
   * Template parameter AT is an allocatable array type, such as 
   * Util::DArray<double> or Pscf::HostDArray<double>, that has an allocate 
   * function that takes an integer capacity as its only parameter.
   *
   * An Exception is thrown if fields is not allocated, or if the fields
   * do not all have the same positive list of dimensions.
   *
   * \param fields  DArray of fields (in)
   * \param dimensions  dimensions mesh for each field (out)
   * \param nMonomer  number of fields, or monomer types (out)
   */
   template <int D, class FT>
   void inspectFields(DArray<FT> const & fields,
                      IntVec<D> & dimensions,
                      int & nMonomer);

   /**
   * Check allocation of a DArray of 1D arrays, allocate if necessary.
   *
   * Template parameter AT is an allocatable array type, such as 
   * Util::DArray<double> or Pscf::HostDArray<double>, that has an allocate 
   * function that takes an integer capacity as its only parameter.
   *
   * On successful exit, the capacity of the arrays container is equal to 
   * the parameter nMonomer, and the capacity for each element of type AT 
   * is equal to the parameter capacity.
   *
   * If the container arrays is allocated on entry, the above conditions are
   * checked, and an Exception is thrown if any are violated.  If arrays is
   * not allocated on entry, it is allocated with the required dimensions.
   *
   * \param arrays  DArray of array of type AT (in/out)
   * \param capacity  required capacity of each array (in)
   * \param nMonomer  required number of arrays, or monomer types (in)
   */
   template <int D, class AT>
   void checkAllocateArrays(DArray< AT >& arrays, 
                            int capacity,
                            int nMonomer);

   /**
   * Inspect dimensions of a DArray of 1D arrays, each of type AT.
   *
   * An Exception is thrown if the arrays container is not allocated, or 
   * if the 1D arrays do not all have the same positive capacity.
   *
   * \param arrays  DArray of arrays (in)
   * \param capacity  capacity of each array (out)
   * \param nMonomer  number of arrays, or monomer types (out)
   */
   template <class AT>
   void inspectArrays(DArray<AT> const & arrays,
                      int & capacity,
                      int & nMonomer);


   // Templates for RGrid data IO

   /**
   * Read mesh dimensions from a field file header.
   * 
   * An Exception is thrown if the mesh dimensions in the file do not
   * match those given by the input vector meshDimensions.
   *
   * \param in  input stream
   * \param meshDimensions  expected values for mesh dimensions
   */
   template <int D>
   void readMeshDimensions(std::istream& in,
                           IntVec<D> const& meshDimensions);

   /**
   * Write mesh dimensions to a field file header.
   * 
   * \param in  input stream
   * \param meshDimensions  vector of integer mesh dimensions
   */
   template <int D>
   void writeMeshDimensions(std::ostream &out,
                            IntVec<D> const& meshDimensions);

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
   * Read the number of basis functions from a basis field file.
   *
   * \param in input stream
   * \return value of nBasis obtained from the file
   */
   int readNBasis(std::istream& in);
   
   /**
   * Write the number of basis functions to a basis field file.
   *
   * \param out output stream
   * \param value  value of nBasis 
   */
   void writeNBasis(std::ostream& in, int nBasis);
   
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

   // Field manipulation utilities

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
   * \param meshDimensions dimensions of original mesh for fields
   * \param unitCell  original crystallographic unit cell
   * \param replicas  number of unit cell replicas in each direction
   */
   template <int D, class AT>
   void replicateUnitCell(std::ostream& out,
                          DArray<AT > const & fields,
                          IntVec<D> const & meshDimensions,
                          UnitCell<D> const & unitCell,
                          IntVec<D> const & replicas);

} // namespace Prdc
} // namespace Pscf
#include "fieldIoUtil.tpp"
#endif

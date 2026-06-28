#ifndef PRDC_CFIELD_IO_H
#define PRDC_CFIELD_IO_H

/*
* PSCF - Polymer Self-Consistent Field
*
* Copyright 2015 - 2025, The Regents of the University of Minnesota
* Distributed under the terms of the GNU General Public License.
*/

#include <pscf/math/IntVec.h>   // Template with a default parameter

// Forward class declarations 
namespace Util {
   template <typename T> class DArray;
}

namespace Pscf {
namespace Prdc {

   using namespace Util;

   // Templates for complex field data IO

   /**
   * Read data for an array of complex fields, with no header section.
   *
   * This function reads the data section of a complex field file for
   * multiple monomer types, with no header section.
   *
   * The template parameter AT must be an array type that provides an
   * overloaded [] subscript operator that returns a complex number of
   * some type CT, where CT supports assignment from real and imaginary
   * parts via a function assign(CT&, double const&, double const&). 
   * Valid complex array types include DArray<fftw_complex> and 
   * HostDArray<cudaReal>.
   *
   * \ingroup Prdc_Field_Module
   *
   * \param in  input file stream
   * \param fields  array of complex fields (out)
   * \param dimensions  vector of mesh dimensions (in)
   */
   template <int D, class AT>
   void readCFieldsData(std::istream& in,
                        DArray< AT >& fields,
                        IntVec<D> const & dimensions);
 
   /**
   * Read data for a single complex field, with no header section.
   *
   * This function reads the data section of a complex field file for
   * a single monomer type, with no header section.
   *
   * The template parameter AT must be an array type that provides an
   * overloaded [] subscript operator that returns a complex number
   * of type that supports assignment from double precision real and
   * imaginary parts, assign(CT& , double const& , double const&).
   *
   * \ingroup Prdc_Field_Module
   *
   * \param in  input file stream
   * \param field  array containing a single complex field (out)
   * \param dimensions  vector of mesh dimensions (in)
   */
   template <int D, class AT>
   void readCFieldData(std::istream& in, 
                      AT& field,
                      IntVec<D> const& dimensions);
 
   /**
   * Write data for an array of complex fields, with no header section.
   *
   * This function writes the data section of the complex field file
   * format for a multiple monomer types, with no header section.
   *
   * The template parameter AT must be an array type that provides an
   * overloaded [] subscript operator that returns a complex number of
   * a type CT for which there exist real(CT const&) and imag(CT const&) 
   * functions that return real and imaginary parts. 
   *
   * \ingroup Prdc_Field_Module
   *
   * \param out  output file stream
   * \param fields  array of complex fields (in)
   * \param dimensions  vector of mesh dimensions (in)
   */
   template <int D, class AT, typename CT, typename RT>
   void writeCFieldsData(std::ostream& out, 
                         DArray< AT > const& fields,
                         IntVec<D> const& dimensions);
 
   /**
   * Write data for a single complex field, with no header section.
   *
   * This function writes the data section of complex field file format
   * for a single monomer type, with no header section.
   *
   * The template parameter AT must be an array type that provides an
   * overloaded [] subscript operator that returns a complex number a
   * type CT for which there exist real(CT const&) and imag(CT const&) 
   * functions that return real and imaginary parts. 
   *
   * \ingroup Prdc_Field_Module
   *
   * \param out  output file stream
   * \param field  array containing a single complex field (out)
   * \param dimensions  vector of mesh dimensions
   */
   template <int D, class AT, typename CT, typename RT>
   void writeCFieldData(std::ostream& out, 
                        AT const& field,
                        IntVec<D> const& dimensions);

} // namespace Prdc
} // namespace Pscf
#include "cFieldIo.tpp"
#endif

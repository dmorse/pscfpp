#ifndef PRDC_CPU_FIELD_BASIS_CONVERTER_H
#define PRDC_CPU_FIELD_BASIS_CONVERTER_H

/*
* PSCF Package 
*
* Copyright 2016 - 2022, The Regents of the University of Minnesota
* Distributed under the terms of the GNU General Public License.
*/

#include <util/containers/DMatrix.h>     // member

namespace Util {
   template <typename T> class DArray;
}

namespace Pscf {
namespace Prdc {
namespace Cpu {

   template <int D> class RField;
   template <int D> class RFieldDft;

   using namespace Util;

   /**
   * Tool for conversion of fields to a basis in composition space.
   *
   * A FieldBasisConverter has an associated orthogonal basis of vectors 
   * in a real vector space of dimension nMonomer, where nMonomer is the
   * number of monomer types.  This class provides functions that can
   * convert between a list of fields associated with different monomer 
   * types to a corresponding list of field components that are each
   * associated with a vector in such a basis. 
   * 
   * This class is designed primarily for use in field theory simulations 
   * in which it useful to define field components associated with
   * eigenvectors of an nMonomer x nMonomer symmetric matrix that 
   * describes interactions among monomers of different types.
   * 
   * \ingroup Prdc_Cpu_Module
   */
   template <int D>
   class FieldBasisConverter 
   {

   public:

      /**
      * Default constructor.
      */
      FieldBasisConverter();

      /**
      * Constructor.
      *
      * The parameter basis is an nMonomer times nMonomer matrix of 
      * orthogonal basis vectors in which the first index (row index) 
      * identifies a basis vector and the second index identifies a 
      * monomer type of a particular vector component. The basis 
      * vectors must be orthogonal, and must be normalized such that 
      * the sum of the squares of elements of each basis vector is 
      * equal to the value of the parameter normSq. 
      *
      * Note that matrix basis is passed by value, and a local copy is made
      * internally during construction. 
      * 
      * \param basis  Matrix of basis vectors
      * \param normSq  Sum of squares of elements of each basis vector
      */
      FieldBasisConverter(DMatrix<double> basis, double normSq);

      /**
      * Constructor.
      *
      * Same as the primary constructor, but with parameter normSq passed 
      * as an integer. This integer is converted internally to a floating
      * point value. 
      * 
      * \param basis  Matrix of eigenvectors
      * \param normSq  Sum of squares of elements of each basis vector
      */
      FieldBasisConverter(DMatrix<double> basis, int normSq);

      /**
      * Destructor.
      */
      virtual ~FieldBasisConverter();

      /**
      * Set or reset the basis after construction.
      *
      * The parameters have same meaning as parameters of non-default 
      * constructors.
      *
      * \param basis  Matrix of vectors.
      * \param normSq  Sum of squares of elements of each basis vector
      */
      void setBasis(DMatrix<double> basis, double normSq);

      /**
      * Check validity (orthogonality and normalization) of the basis.
      *
      * This function returns the maximum absolute magnitude of the 
      * difference between the dot product of any two basis vectors and 
      * the expected value. The expected value is zero for any two distinct 
      * basis vectors and normSq for a dot product of a basis vector with 
      * itself.
      */
      double maxBasisError() const;

      /**
      * Convert a set of monomer fields to field basis components.
      *
      * The value of the component associated with a specific basis 
      * vector at each grid point (parameter out) is obtained by taking 
      * a dot product of that basis vector with the vector of values of 
      * fields associated with different monomer types (parameter in) at 
      * the same grid point, and dividing that dot product by normSq.
      *
      * \param in input fields, indexed by monomer type id
      * \param out output field components, indexed by basis vector id
      */
      void convertToBasis(DArray< RField<D> > const & in, 
                          DArray< RField<D> > & out) const;

      /**
      * Convert a set of field basis components to monomer fields.
      *
      * The vector of nMonomer values of output fields (parameter out) at 
      * each grid point is given by a linear superposition of basis vectors
      * with component or coefficient values given by the values of the 
      * input fields components (parameter in) at the the same grid point.
      *
      * \param in input field components, indexed by basis vector id
      * \param out output fields, indexed by monomer type id
      */
      void convertFromBasis(DArray< RField<D> > const & in, 
                            DArray< RField<D> > & out) const;

   private:

      /**
      * Basis in real vector space of dimension nMonomer.
      *
      * The matrix basis_ must be a square matrix of dimension
      * nMonomer_, in which the first (row) index identifies a 
      * a basis vector, and the second (column) index identifies
      * a monomer type. Vectors in this basis must be orthogonal.
      */
      DMatrix<double> basis_;

      /**
      * Square norm of each basis vector in basis_.
      *
      * The sum of the squares of the elements of every vector in
      * the basis must be equal to normSq_.
      */
      double normSq_;

      /**
      * Number of monomer types (dimension of basis).
      */
      int nMonomer_;

   };

   #ifndef PRDC_CPU_FIELD_BASIS_CONVERTER_TPP
   // Suppress implicit instantiation
   extern template class FieldBasisConverter<1>;
   extern template class FieldBasisConverter<2>;
   extern template class FieldBasisConverter<3>;
   #endif

} // namespace Pscf::Pspc::Cpu
} // namespace Pscf::Pspc
} // namespace Pscf
#endif

#ifndef PRDC_CPU_FIELD_BASIS_CONVERTER_H
#define PRDC_CPU_FIELD_BASIS_CONVERTER_H

/*
* PSCF Package 
*
* Copyright 2016 - 2022, The Regents of the University of Minnesota
* Distributed under the terms of the GNU General Public License.
*/

#include <util/containers/DMatrix.h>     // member
//#include <util/containers/DArray.h>      

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
   * Tool for conversion of fields to a basis in composition space
   * 
   * \ingroup Prdc_Cpu_Module
   */
   template <int D>
   class FieldBasisConverter 
   {

   public:

      /**
      * Constructor.
      *
      * The parameter basis is an nMonomer times nMonomer matrix of 
      * orthogonal basis vectors in which the first index (row index) 
      * identifies a basis vectors and the second indices the monomer 
      * associated with a particular vector component. The basis 
      * vectors must be normalized such that the sum of the squares of 
      * elements of each basis vector is equal to parameter normSq. 
      *
      * The matrix basis is passed by value, and a local copy is made
      * internally.
      * 
      * \param basis  Matrix of eigenvectors
      * \param normSq  Sum of squares of elements of each basis vector
      */
      FieldBasisConverter(DMatrix<double> basis, double normSq);

      /**
      * Constructor.
      *
      * Same as the primary constructor, but with normSq passed as an
      * integer that is converted to a floating point number internally.
      * 
      * \param basis  Matrix of eigenvectors
      * \param normSq  Sum of squares of elements of each basis vector
      */
      FieldBasisConverter(DMatrix<double> basis, int normSq);

      /**
      * Destructor.
      */
      virtual ~FieldBasisConverter();

      void convertToBasis(DArray< RField<D> >  const & in, 
                          DArray< RField<D> > & out) const;

      void convertFromBasis(DArray< RField<D> > const & in, 
                            DArray< RField<D> > & out) const;


   private:

      /*
      * Basis in real vector space of dimension nMonomer.
      * First index identifies a basis vector, 2nd is monomer type.
      */
      DMatrix<double> basis_;

      /*
      * Square norm of each basis vector in basis_.
      */
      double normSq_;

      /*
      * Number of monomer types.
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

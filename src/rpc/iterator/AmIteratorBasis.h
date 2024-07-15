#ifndef RPC_AM_ITERATOR_BASIS_H
#define RPC_AM_ITERATOR_BASIS_H

/*
* PSCF - Polymer Self-Consistent Field Theory
*
* Copyright 2016 - 2022, The Regents of the University of Minnesota
* Distributed under the terms of the GNU General Public License.
*/

#include <rpc/iterator/Iterator.h>               // base class argument
#include <rpc/iterator/ImposedFieldsGenerator.h> // member variable
#include <pscf/iterator/AmIteratorTmpl.h>        // base class template
#include <pscf/iterator/AmbdInteraction.h>       // member variable
#include <util/containers/DArray.h>              // base class argument


namespace Pscf {
namespace Rpc
{

   template <int D> class System;

   using namespace Util;

   /**
   * Rpc implementation of the Anderson Mixing iterator with symmetry.
   *
   * \ingroup Rpc_Iterator_Module
   */
   template <int D>
   class AmIteratorBasis
      : public AmIteratorTmpl< Iterator<D>, DArray<double> >
   {

   public:

      /**
      * Constructor.
      *
      * \param system System object associated with this iterator.
      */
      AmIteratorBasis(System<D>& system);

      /**
      * Destructor.
      */
      ~AmIteratorBasis();

      /**
      * Read all parameters and initialize.
      *
      * \param in input filestream
      */
      void readParameters(std::istream& in);

      /**
      * Return compressor times contributions.
      *
      * \param out  output stream for timer report
      */
      void outputTimers(std::ostream& out);

      /**
      * Return specialized sweep parameter types to add to the Sweep object
      */
      DArray<ParameterType> getParameterTypes();

      /**
      * Set the value of a specialized sweep parameter
      * 
      * \param name  name of the specialized parameter
      * \param ids  array of integer indices specifying the value to set
      * \param value  the value to which the parameter is set
      * \param success  boolean flag used to indicate if parameter was set
      */
      void setParameter(std::string name, DArray<int> ids, 
                                          double value, bool& success);

      /**
      * Get the value of a specialized sweep parameter
      * 
      * \param name  name of the specialized parameter
      * \param ids  array of integer indices specifying the value to get
      * \param success  boolean flag used to indicate if parameter was gotten
      */
      double getParameter(std::string name, DArray<int> ids, bool& success)
      const;

      // Inherited public member functions
      using AmIteratorTmpl<Iterator<D>, DArray<double> >::solve;
      using AmIteratorTmpl<Iterator<D>, DArray<double> >::clearTimers;
      using Iterator<D>::isFlexible;
      using Iterator<D>::flexibleParams;
      using Iterator<D>::setFlexibleParams;
      using Iterator<D>::nFlexibleParams;

   protected:

      // Inherited protected members
      using ParamComposite::readOptional;
      using ParamComposite::readParamCompositeOptional;
      using ParamComposite::readOptionalFSArray;
      using ParamComposite::setClassName;
      using AmIteratorTmpl< Iterator<D>, DArray<double> >::verbose;
      using Iterator<D>::system;
      using Iterator<D>::isFlexible_;
      using Iterator<D>::flexibleParams_;

      /**
      * Setup iterator just before entering iteration loop.
      *
      * \param isContinuation Is this a continuation within a sweep?
      */
      void setup(bool isContinuation);

   private:

      /// ImposedFieldsGenerator object
      ImposedFieldsGenerator<D> imposedFields_;

      /// Local copy of interaction, adapted for use AMBD residual definition
      AmbdInteraction interaction_;

      /// How are stress residuals scaled in error calculation?
      double scaleStress_;

      /**
      * Assign one field to another.
      *
      * \param a the field to be set (lhs of assignment)
      * \param b the field for it to be set to (rhs of assigment)
      */
      void setEqual(DArray<double>& a, DArray<double> const & b);

      /**
      * Compute the inner product of two vectors
      */
      double dotProduct(DArray<double> const & a, DArray<double> const & b);

      /**
      * Find the maximum magnitude element of a residual vector.
      */
      double maxAbs(DArray<double> const & hist);

      /**
      * Update the basis for residual or field vectors.
      *
      * \param basis RingBuffer of residual or field basis vectors
      * \param hists RingBuffer of past residual or field vectors
      */
      void updateBasis(RingBuffer<DArray<double> > & basis,
                       RingBuffer<DArray<double> > const & hists);

      /**
      * Add linear combination of basis vectors to trial field.
      *
      * \param trial trial vector (input-output)
      * \param basis RingBuffer of basis vectors
      * \param coeffs array of coefficients of basis vectors
      * \param nHist number of histories stored at this iteration
      */
      void addHistories(DArray<double>& trial,
                        RingBuffer<DArray<double> > const & basis,
                        DArray<double> coeffs,
                        int nHist);

      /**
      * Add predicted error to field trial.
      *
      * \param fieldTrial trial field (in-out)
      * \param resTrial predicted error for current trial
      * \param lambda Anderson-Mixing mixing
      */
      void addPredictedError(DArray<double>& fieldTrial,
                             DArray<double> const & resTrial,
                             double lambda);

      /**
      * Does the system has an initial guess for the field?
      */
      bool hasInitialGuess();

      /**
      * Compute and returns the number of elements in field vector.
      *
      * Called during allocation and then stored.
      */
      int nElements();

      /**
      * Gets the current field vector from the system.
      *
      * \param curr current field vector
      */
      void getCurrent(DArray<double>& curr);

      /**
      * Have the system perform a computation using new field.
      *
      * Solves the modified diffusion equations, computes concentrations,
      * and optionally computes stress components.
      */
      void evaluate();

      /**
      * Compute the residual vector.
      *
      * \param resid current residual vector value
      */
      void getResidual(DArray<double>& resid);

      /**
      * Updates the system field with the new trial field.
      *
      * \param newGuess trial field vector
      */
      void update(DArray<double>& newGuess);

      /**
      * Outputs relevant system details to the iteration log.
      */
      void outputToLog();

   };

   #ifndef RPC_AM_ITERATOR_BASIS_TPP
   // Suppress implicit instantiation
   extern template class AmIteratorBasis<1>;
   extern template class AmIteratorBasis<2>;
   extern template class AmIteratorBasis<3>;
   #endif

} // namespace Rpc
} // namespace Pscf
#endif

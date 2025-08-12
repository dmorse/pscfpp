#ifndef PRDC_UNIT_CELL_BASE_H
#define PRDC_UNIT_CELL_BASE_H

/*
* PSCF - Polymer Self-Consistent Field
*
* Copyright 2015 - 2025, The Regents of the University of Minnesota
* Distributed under the terms of the GNU General Public License.
*/

#include <pscf/math/RealVec.h>         // member
#include <pscf/math/IntVec.h>          // interface, with default parameters
#include <util/containers/FArray.h>    // member
#include <util/containers/FMatrix.h>   // member
#include <util/containers/FSArray.h>   // interface

// Forward declarations
namespace Util {
   template <typename T> class Signal;
   template <> class Signal<void>;
}

namespace Pscf {
namespace Prdc {

   using namespace Util;

   /**
   * Base class template for a crystallographic unit cell.
   *
   * \ingroup Prdc_Crystal_Module
   */
   template <int D>
   class UnitCellBase
   {
   public:

      /**
      * Constructor.
      */
      UnitCellBase();

      /**
      * Destructor.
      */
      ~UnitCellBase();

      /// \name Unit Cell Parameters
      ///@{

      /**
      * Set all the parameters of unit cell.
      *
      * The lattice system must already be set to a non-null value.
      *
      * \param parameters array of unit cell parameters
      */
      void setParameters(FSArray<double, 6> const & parameters);

      /**
      * Get the number of parameters in the unit cell.
      */
      int nParameter() const;

      /**
      * Get the parameters of this unit cell.
      */
      FSArray<double, 6> parameters() const;

      /**
      * Get a single parameter of this unit cell.
      *
      * \param i array index of the desired parameter
      */
      double parameter(int i) const;

      /**
      * Has this unit cell been initialized?
      *
      * A unit cell is initialized after both the lattice type and
      * the lattice parameters have been set.
      */
      bool isInitialized() const;

      ///@}
      /// \name Bravais Lattice Data
      ///@{
 
      /**
      * Get Bravais basis vector i, denoted by a_i.
      *
      * \param i  array index of the desired Bravais basis vector
      */
      const RealVec<D>& rBasis(int i) const;

      /**
      * Get reciprocal basis vector i, denoted by b_i.
      *
      * \param i  array index of the desired reciprocal basis vector
      */
      const RealVec<D>& kBasis(int i) const;

      /**
      * Get component j of derivative of rBasis vector a_i w/respect to k.
      *
      * \param k  index of cell parameter
      * \param i  index of the desired basis vector a_i
      * \param j  index of a Cartesian component of a_i
      */
      double drBasis(int k, int i, int j) const;

      /**
      * Get component j of derivative of kBasis vector b_i w/respect to k.
      *
      * \param k index of cell parameter
      * \param i array index of the desired reciprocal basis vector b_i
      * \param j index of a Cartesian component of b_i
      */
      double dkBasis(int k, int i, int j) const;

      /**
      * Get derivative of dot product ai.aj with respect to parameter k.
      *
      * \param k  index of cell parameter
      * \param i  array index of 1st Bravais basis vector a_i
      * \param j  array index of 2nd Bravais basis vector a_i
      */
      double drrBasis(int k, int i, int j) const;

      /**
      * Get derivative of dot product bi.bj with respect to parameter k.
      *
      * \param k  index of cell parameter
      * \param i  array index of 1st reciprocal basis vector b_i
      * \param j  array index of 2nd reciprocal basis vector b_i
      */
      double dkkBasis(int k, int i, int j) const;

      ///@}
      /// \name Wavevector Properties
      ///@{

      /**
      * Compute square magnitude of reciprocal lattice vector.
      *
      * \param k  vector of components of a reciprocal lattice vector
      */
      virtual double ksq(IntVec<D> const & k) const;

      /**
      * Compute derivative of square wavevector w/respect to cell parameter.
      *
      * This function computes and returns a derivative with respect to 
      * unit cell parameter number n of the square of a reciprocal lattice
      * vector with integer coefficients given by the elements of vec.
      *
      * \param vec  vector of components of a reciprocal lattice vector
      * \param n  index of a unit cell parameter
      */
      virtual double dksq(IntVec<D> const & vec, int n) const;

      ///@}
      /// \name Signal Interface
      ///@{

      /**
      * Associating an externally defined signal with this unit cell.
      *
      * \param signal  Signal triggered by unit cell modification
      */
      void setSignal(Signal<void>& signal);

      /**
      * Does this object have an associated Signal<void>?
      */
      bool hasSignal() const;

      /**
      * Get the associated Signal by non-const reference.
      */
      Signal<void>& signal();

      ///@}

   protected:

      /**
      * Array of Bravais lattice basis vectors.
      */
      FArray<RealVec<D>, D> rBasis_;

      /**
      * Array of reciprocal lattice basis vectors.
      */
      FArray<RealVec<D>, D> kBasis_;

      /**
      * Array of derivatives of rBasis.
      *
      * Element drBasis_[k](i,j) is the derivative with respect to
      * parameter k of component j of Bravais basis vector i.
      */
      FArray<FMatrix<double, D, D>, 6> drBasis_;

      /**
      * Array of derivatives of kBasis
      *
      * Element dkBasis_[k](i,j) is the derivative with respect to
      * parameter k of component j of reciprocal basis vector i.
      */
      FArray<FMatrix<double, D, D>, 6> dkBasis_;

      /**
      * Array of derivatives of a_i.a_j
      *
      * Element drrBasis_[k](i,j) is the derivative with respect
      * to parameter k of the dot product (a_i.a_j) of Bravais
      * lattice basis vectors a_i and a_j.
      */
      FArray<FMatrix<double, D, D>, 6> drrBasis_;

      /**
      * Array of derivatives of b_i.b_j
      *
      * Element dkkBasis_[k](i,j) is the derivative with respect
      * to parameter k of the dot product (b_i.b_j) of reciprocal
      * lattice basis vectors b_i and b_j.
      */
      FArray<FMatrix<double, D, D>, 6> dkkBasis_;

      /**
      * Parameters used to describe the unit cell.
      */
      FArray<double, 6> parameters_;

      /**
      * Number of parameters required to specify unit cell.
      */
      int nParameter_;

      /**
      * Has this unit cell been fully initialized?
      */
      bool isInitialized_;

      /**
      * Compute all protected data, given latticeSystem and parameters.
      *
      * All functions that reset unit cell parameters must call this 
      * function to reset depenent data, including stream insertion 
      * operators.
      *
      * Calls initializeToZero, setBasis, & computeDerivatives functions
      * internally. Also sets isInitialized flag true and notifies any 
      * observers.
      */
      void setLattice();

   private:

      /**
      * Pointer to a Signal<void> subobject.
      */
      Signal<void>* signalPtr_;

      /**
      * Initialize all arrays to zero.
      *
      * Sets all elements of the following arrays to zero:
      * rBasis_, kBasis_, drBasis_, dkBasis, drrBasis_ and dkkBasis_.
      */
      void initializeToZero();

      /**
      * Set values of rBasis, kBasis, and drBasis.
      *
      * Invoke initializeToZero before this, computeDerivatives after.
      *
      * \pre Lattice system, nParam and parameters must be set.
      */
      virtual void setBasis() = 0;

      /**
      * Compute values of dkBasis_, drrBasis_, and dkkBasis_.
      */
      void computeDerivatives();

      /**
      * Copy constructor - private and unimplemented to prohibit.
      */
      UnitCellBase(UnitCellBase<D> const & other);

      /**
      * Assignment operator - private and unimplemented to prohibit.
      */
      UnitCellBase<D>& operator = (UnitCellBase<D> const & other);

   };

   // Inline member functions

   /*
   * Get the number of Parameters in the unit cell.
   */
   template <int D>
   inline
   int UnitCellBase<D>::nParameter() const
   {  return nParameter_;  }

   /*
   * Get a single parameter of the unit cell.
   */
   template <int D>
   inline
   double UnitCellBase<D>::parameter(int i) const
   {  return parameters_[i]; }

   /*
   * Get Bravais basis vector i.
   */
   template <int D>
   inline
   const RealVec<D>& UnitCellBase<D>::rBasis(int i) const
   {  return rBasis_[i];  }

   /*
   * Get reciprocal basis vector i.
   */
   template <int D>
   inline
   const RealVec<D>& UnitCellBase<D>::kBasis(int i) const
   {  return kBasis_[i];  }

   /*
   * Get component j of derivative of r basis vector i w/respect to param k.
   */
   template <int D>
   inline
   double UnitCellBase<D>::drBasis(int k, int i, int j) const
   {  return drBasis_[k](i,j);  }

   /*
   * Get component j of derivative of r basis vector i w/respect to param k.
   */
   template <int D>
   inline
   double UnitCellBase<D>::dkBasis(int k, int i, int j) const
   {  return dkBasis_[k](i, j);  }

   /*
   * Get the derivative of ki*kj with respect to parameter k.
   */
   template <int D>
   inline
   double UnitCellBase<D>::dkkBasis(int k, int i, int j) const
   {  return dkkBasis_[k](i, j);  }

   /*
   * Get the derivative of ri*rj with respect to parameter k.
   */
   template <int D>
   inline
   double UnitCellBase<D>::drrBasis(int k, int i, int j) const
   {  return drrBasis_[k](i, j);  }

   /*
   * Is this unit cell initialized?
   */
   template <int D>
   inline
   bool UnitCellBase<D>::isInitialized() const
   {  return isInitialized_; }

   #ifndef PRDC_UNIT_CELL_BASE_TPP
   // Suppress implicit instantiation of base class
   extern template class UnitCellBase<1>;
   extern template class UnitCellBase<2>;
   extern template class UnitCellBase<3>;
   #endif

}
}
#endif

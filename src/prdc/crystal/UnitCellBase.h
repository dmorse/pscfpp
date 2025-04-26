#ifndef PRDC_UNIT_CELL_BASE_H
#define PRDC_UNIT_CELL_BASE_H

/*
* PSCF - Polymer Self-Consistent Field Theory
*
* Copyright 2015 - 2025, The Regents of the University of Minnesota
* Distributed under the terms of the GNU General Public License.
*/

#include <util/signal/Signal.h>
#include <util/containers/FArray.h>
#include <util/containers/FSArray.h>
#include <util/containers/FMatrix.h>
#include <pscf/math/IntVec.h>
#include <pscf/math/RealVec.h>
#include <util/math/Constants.h>

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
      /// \name Observer Interface
      ///@{

      /**
      * Add an Observer to this unit cell.
      *
      * This adds an instance of class "Observer" and a void member
      * function of that class to the list of observers associated with
      * this UnitCellBase<D> object. The methodPtr argument must be the
      * relative address of a void function that takes no arguments.
      * 
      * Whenever the parameters of this unit cell are set or modified
      * all observers are notified by calling the prescribed member
      * function of each observer.
      *
      * <b> Usage </b> : Suppose Observer is a class that has a member 
      * function named "receive" that should be called when the unit
      * cell parameters are set or updated. The "receive" function must
      * have a signature "void receive()". An instance of that class
      * may be added to the list of observers for a UnitCellBase<3>
      * using the following syntax:
      * \code
      *    UnitCellBase<3> unitCell;
      *    Observer observer;
      *    unitCell.addObserver(observer, &Observer::receive);
      * \endcode
      * Alternatively, one could explicitly create a pointer to the 
      * member function and then pass that to the unit cell, like this:
      * \code
      *    UnitCellBase<3> unitCell;
      *
      *    Observer observer;
      *    void (Observer::*functionPtr)() = nullptr;
      *    functionPtr = &Observer::receive;
      *
      *    v.addObserver(observer, functionPtr);
      * \endcode
      *
      * \param observer  observer object (invokes member function)
      * \param methodPtr  pointer to relevant member function
      */
      template <class Observer>
      void addObserver(Observer& observer, void (Observer::*methodPtr)());

      /**
      * Clear all observers and delete the Signal<> object.
      */
      void clearObservers();

      /**
      * Return the current number of observers.
      *
      * The number of observers is zero upon construction, is incremented
      * by 1 every time addObserver() is called, and is reset to 0 by
      * calling clearObservers().
      */
      int nObserver() const;

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
      * Does this object have a Signal<> subobject?
      */
      bool hasSignal() const;

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

   /*
   * Get the number of Parameters in the unit cell.
   */
   template <int D>
   inline
   int UnitCellBase<D>::nParameter() const
   {  return nParameter_;  }

   /*
   * Get the unit cell parameters.
   */
   template <int D>
   inline
   FSArray<double, 6> UnitCellBase<D>::parameters() const
   {
      FSArray<double, 6> parameters;
      for (int i = 0; i < nParameter_; ++i) {
         parameters.append(parameters_[i]);
      }
      return parameters;
   }

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

   // Non-inline member functions

   /*
   * Constructor.
   */
   template <int D>
   UnitCellBase<D>::UnitCellBase()
    : nParameter_(0),
      isInitialized_(false),
      signalPtr_(nullptr)
   {  initializeToZero(); }

   /*
   * Destructor.
   */
   template <int D>
   UnitCellBase<D>::~UnitCellBase()
   {
      if (hasSignal()) {
         delete signalPtr_;
      }
   }

   /*
   * Set all the parameters in the unit cell.
   */
   template <int D>
   void UnitCellBase<D>::setParameters(FSArray<double, 6> const& parameters)
   {
      UTIL_CHECK(parameters.size() == nParameter_);
      isInitialized_ = false;
      for (int i = 0; i < nParameter_; ++i) {
         parameters_[i] = parameters[i];
      }
      setLattice();
   }
   
   /*
   * Get square magnitude of reciprocal basis vector.
   */
   template <int D>
   double UnitCellBase<D>::ksq(IntVec<D> const & k) const
   {
      RealVec<D> g(0.0);
      RealVec<D> p;
      for (int i = 0; i < D; ++i) {
         p.multiply(kBasis_[i], k[i]);
         g += p;
      }
      double value = 0.0;
      for (int i = 0; i < D; ++i) {
         value += g[i]*g[i];
      }
      return value;
   }

   /*
   * Get magnitude of derivative of square of reciprocal basis vector.
   */
   template <int D>
   double UnitCellBase<D>::dksq(IntVec<D> const & vec, int n) const
   {
      double element = 0.0;
      double value = 0.0;

      for (int p = 0; p < D; ++p){
         for (int q = 0; q < D; ++q){
            element = dkkBasis(n, p, q);
            value += vec[p]*vec[q]*element;
         }
      }

      return value;
   }

   // Functions to manage observers

   template <int D>
   template <class Observer>
   void UnitCellBase<D>::addObserver(Observer& observer,
                                     void (Observer::*methodPtr)())
   {
      if (!hasSignal()) {
         signalPtr_ = new Signal<void>;
      }
      signalPtr_->addObserver(observer, methodPtr);
   }

   template <int D>
   void UnitCellBase<D>::clearObservers() 
   {
      if (hasSignal()) {
         signalPtr_->clear();
         delete signalPtr_;
         signalPtr_ = nullptr;
      }
   }

   template <int D>
   int UnitCellBase<D>::nObserver() const
   {
      if (hasSignal()) {
         return signalPtr_->nObserver();
      } else {
         return 0;
      }
   }

   // Protected member functions

   /*
   * Initialize internal arrays to zero.
   */
   template <int D>
   void UnitCellBase<D>::initializeToZero()
   {
      // Initialize all elements to zero
      int i, j, k;
      for (i = 0; i < D; ++i) {
         for (j = 0; j < D; ++j) {
            rBasis_[i][j] = 0.0;
            kBasis_[i][j] = 0.0;
         }
      }
      for (k = 0; k < 6; ++k){
         for (i = 0; i < D; ++i) {
            for (j = 0; j < D; ++j) {
               drBasis_[k](i,j) = 0.0;
               dkBasis_[k](i,j) = 0.0;
               drrBasis_[k](i,j) = 0.0;
               dkkBasis_[k](i,j) = 0.0;
            }
         }
      }
   }

   /*
   * Compute quantities involving derivatives.
   */
   template <int D>
   void UnitCellBase<D>::computeDerivatives()
   {
      // Compute dkBasis
      int p, q, r, s, t;
      for (p = 0; p < nParameter_; ++p) {
         for (q = 0; q < D; ++q) {
            for (r = 0; r < D; ++r) {

               // Loop over free indices s, t
               for (s = 0; s < D; ++s) {
                  for (t = 0; t < D; ++t) {
                     dkBasis_[p](q,r)
                       -= kBasis_[q][s]*drBasis_[p](t,s)*kBasis_[t][r];
                  }
               }
               dkBasis_[p](q,r) /= 2.0*Constants::Pi;

            }
         }
      }

      // Compute drrBasis and dkkBasis
      for (p = 0; p < nParameter_; ++p) {
         for (q = 0; q < D; ++q) {
            for (r = 0; r < D; ++r) {
               for (s = 0; s < D; ++s) {
                  drrBasis_[p](q,r) += rBasis_[q][s]*drBasis_[p](r,s);
                  drrBasis_[p](q,r) += rBasis_[r][s]*drBasis_[p](q,s);
                  dkkBasis_[p](q,r) += kBasis_[q][s]*dkBasis_[p](r,s);
                  dkkBasis_[p](q,r) += kBasis_[r][s]*dkBasis_[p](q,s);

               }
            }
         }
      }

   }

   /*
   * Set all lattice parameters.
   */
   template <int D>
   void UnitCellBase<D>::setLattice()
   {
      initializeToZero();
      setBasis();
      computeDerivatives();
      isInitialized_ = true;
      if (hasSignal()) {
         signalPtr_->notify();
      }
   }

   template <int D>
   bool UnitCellBase<D>::hasSignal() const
   {  return bool(signalPtr_); }

}
}
#endif

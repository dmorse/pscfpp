#ifndef RPG_AM_COMPRESSOR_TPP
#define RPG_AM_COMPRESSOR_TPP

/*
* PSCF - Polymer Self-Consistent Field
*
* Copyright 2015 - 2025, The Regents of the University of Minnesota
* Distributed under the terms of the GNU General Public License.
*/

#include "AmCompressor.h"
#include <rpg/system/System.h>
#include <prdc/cuda/RField.h>
#include <prdc/cuda/resources.h>
#include <pscf/math/IntVec.h>
#include <util/global.h>

namespace Pscf {
namespace Rpg{

   using namespace Util;

   // Constructor
   template <int D>
   AmCompressor<D>::AmCompressor(System<D>& system)
    : Compressor<D>(system),
     isAllocated_(false)
   { ParamComposite::setClassName("AmCompressor"); }

   // Destructor
   template <int D>
   AmCompressor<D>::~AmCompressor()
   {}

   // Read parameters from file
   template <int D>
   void AmCompressor<D>::readParameters(std::istream& in)
   {
      // Call parent class readParameters
      AmIteratorTmpl<Compressor<D>,
                     DeviceArray<cudaReal> >::readParameters(in);
      AmIteratorTmpl<Compressor<D>,
                     DeviceArray<cudaReal> >::readErrorType(in);
   }


   // Initialize just before entry to iterative loop.
   template <int D>
   void AmCompressor<D>::setup(bool isContinuation)
   {
      const int nMonomer = system().mixture().nMonomer();
      const int meshSize = system().domain().mesh().size();
      const IntVec<D> dimensions = system().domain().mesh().dimensions();

      // Allocate memory required by AM algorithm if not done earlier.
      AmIteratorTmpl<Compressor<D>,
                     DeviceArray<cudaReal> >::setup(isContinuation);

      // Allocate memory required by compressor if not done earlier.
      if (!isAllocated_) {
         w0_.allocate(nMonomer);
         wFieldTmp_.allocate(nMonomer);
         for (int i = 0; i < nMonomer; ++i) {
            w0_[i].allocate(dimensions);
            wFieldTmp_[i].allocate(dimensions);
         }
         isAllocated_ = true;
      }

      // Store value of initial guess chemical potential fields
      for (int i = 0; i < nMonomer; ++i) {
         VecOp::eqV(w0_[i], system().w().rgrid(i));
      }
   }

   template <int D>
   int AmCompressor<D>::compress()
   {
      int solve = Base::solve();
      //mdeCounter_ = Base::totalItr();
      return solve;
   }

   template<int D>
   double AmCompressor<D>::computeLambda(double r)
   {
      return 1.0;
   }

   // Private virtual functions that interact with the parent System

   /*
   * Compute and return the number of elements in a field vector.
   */
   template <int D>
   int AmCompressor<D>::nElements()
   {  return system().domain().mesh().size(); }

   /*
   * Does the system have an initial field guess?
   */
   template <int D>
   bool AmCompressor<D>::hasInitialGuess()
   {  return system().w().hasData(); }

   /*
   * Get the current field from the system.
   */
   template <int D>
   void AmCompressor<D>::getCurrent(DeviceArray<cudaReal>& curr)
   {
      /*
      * The field that we are adjusting is the Langrange multiplier field
      * with number of grid pts components.The current value is the difference
      * between w and w0_ for the first monomer (any monomer should give the
      * same answer)
      */
      VecOp::subVV(curr, system().w().rgrid(0), w0_[0]);
   }

   /*
   * Perform the main system computation (solve the MDE).
   */
   template <int D>
   void AmCompressor<D>::evaluate()
   {
      system().compute();
      ++mdeCounter_;
   }

   /*
   * Compute the residual for the current system state.
   */
   template <int D>
   void AmCompressor<D>::getResidual(DeviceArray<cudaReal>& resid)
   {
      const int nMonomer = system().mixture().nMonomer();

      // Initialize resid to c field of species 0 minus 1
      VecOp::subVS(resid, system().c().rgrid(0), 1.0);

      // Add other c fields to get SCF residual vector elements
      for (int i = 1; i < nMonomer; i++) {
         VecOp::addEqV(resid, system().c().rgrid(i));
      }
   }

   /*
   * Update the current system field coordinates.
   */
   template <int D>
   void AmCompressor<D>::update(DeviceArray<cudaReal>& newGuess)
   {
      // Convert back to field format
      const int nMonomer = system().mixture().nMonomer();

      // New field is the w0_ + the newGuess for the Lagrange multiplier field
      for (int i = 0; i < nMonomer; i++) {
         VecOp::addVV(wFieldTmp_[i], w0_[i], newGuess);
      }

      // set system r grid
      system().w().setRGrid(wFieldTmp_);
   }

   template<int D>
   void AmCompressor<D>::outputToLog()
   {}

   template<int D>
   void AmCompressor<D>::outputTimers(std::ostream& out) const
   {
      // Output timing results, if requested.
      out << "\n";
      out << "Compressor times contributions:\n";
      Base::outputTimers(out);
   }

   template<int D>
   void AmCompressor<D>::clearTimers()
   {
      Base::clearTimers();
      mdeCounter_ = 0;
   }

   // Private virtual functions for vector math

   /*
   * Assign one array to another.
   */
   template <int D>
   void AmCompressor<D>::setEqual(DeviceArray<cudaReal>& a,
                                  DeviceArray<cudaReal> const & b)
   {
      UTIL_CHECK(b.capacity() == a.capacity());
      VecOp::eqV(a, b);
   }

   /*
   * Compute and return inner product of two vectors.
   */
   template <int D>
   double AmCompressor<D>::dotProduct(DeviceArray<cudaReal> const & a,
                                      DeviceArray<cudaReal> const & b)
   {
      UTIL_CHECK(a.capacity() == b.capacity());
      return Reduce::innerProduct(a, b);
   }

   /*
   * Compute and return maximum element of a vector.
   */
   template <int D>
   double AmCompressor<D>::maxAbs(DeviceArray<cudaReal> const & a)
   {  return Reduce::maxAbs(a); }

   /*
   * Compute the vector difference a = b - c
   */
   template <int D>
   void AmCompressor<D>::subVV(DeviceArray<cudaReal>& a,
                               DeviceArray<cudaReal> const & b,
                               DeviceArray<cudaReal> const & c)
   {
      VecOp::subVV(a, b, c);
   }

   /*
   * Composite a += b*c for vectors a and b, scalar c
   */
   template <int D>
   void AmCompressor<D>::addEqVc(DeviceArray<cudaReal>& a,
                                 DeviceArray<cudaReal> const & b,
                                 double c)
   {
      UTIL_CHECK(a.capacity() == b.capacity());
      VecOp::addEqVc(a, b, c);
   }

}
}
#endif

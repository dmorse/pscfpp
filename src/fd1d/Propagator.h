#ifndef FD1D_PROPAGATOR_H
#define FD1D_PROPAGATOR_H

/*
* PFTS - Polymer Field Theory Simulator
*
* Copyright 2013, David Morse (morse012@.umn.edu)
* Distributed under the terms of the GNU General Public License.
*/

#include <chem/PropagatorTmpl.h>       // base class template
#include <fd1d/TridiagonalSolver.h>    // member
#include <util/containers/DArray.h>    // member template

namespace Fd1d{ 

   using namespace Util;
   using namespace Chem;

   class Propagator : public PropagatorTmpl<Propagator>
   {

   public:

      // Public typedefs

      /**
      * Chemical potential field type.
      */
      typedef DArray<double> WField;

      /**
      * Monomer concentration field type.
      */
      typedef DArray<double> CField;

      /**
      * Propagator q-field type.
      */
      typedef DArray<double> QField;

      // Member functions

      /**
      * Constructor.
      */
      Propagator();

      /**
      * Destructor.
      */
      ~Propagator();

      /**
      * Initialize grid and allocate required memory.
      *
      * \param xMin minimum value of space coordinate
      * \param xMax minimum value of space coordinate
      * \param nx number of spatial grid points
      * \param ns number of contour steps, # grid points - 1
      */
      void setGrid(double xMin, double xMax, int nx, int ns);

      /**
      * Solve the modified diffusion equation (MDE) for this block.
      *
      * This function computes an initial QField at the head of this
      * block, and then solves the modified diffusion equation for 
      * the block to propagate from the head to the tail. The initial
      * QField at the head is computed by pointwise multiplication of
      * of the tail QFields of all source propagators.
      *
      * \param w chemical potential field for relevant monomer type
      */
      void solve(const WField& w);
  
      /**
      * Solve the MDE for a specified initial condition.
      *
      * This function solves the modified diffusion equation for this
      * block with a specified initial condition, which is given by 
      * head parameter of the function. The function is intended for 
      * use in testing.
      *
      * \param w chemical potential field for relevant monomer type
      * \param head initial condition of QField at head of block
      */
      void solve(const WField& w, const QField& head);
  
      /**
      * Compute unnormalized concentration for block by integration.
      *
      * Upon return, grid point r of the integral array contains the 
      * integral int ds q(r,s)q^{*}(r,L-s), q(r,s) is the solution 
      * obtained from this propagator, q^{*} is the solution of the
      * partner propagator, which solves the MDE for the same block
      * in the opposite direction, and s is contour variable that is
      * integrated over the domain 0 < s < L, where L is the length
      * of the block. 
      *
      * \param prefactor multiplying integral
      * \param integral contour integral of propagator product.
      */ 
      void integrate(double prefactor, CField& integral);

      /**
      * Compute and return partition function for the molecule.
      *
      * This function computes the partition function Q for the 
      * molecule as a spatial average of the initial/head Qfield 
      * for this propagator and the final/tail Qfield of its
      * partner. 
      */ 
      double computeQ();

      /**
      * Return q-field at beginning of block (initial condition).
      */
      const QField& head() const;

      /**
      * Return q-field at end of block.
      */
      const QField& tail() const;

   protected:

      /**
      * Set Crank-Nicholson solver before main integration loop.
      */
      void setupSolver(const WField& w);

      /**
      * Compute initial QField at head from tail QFields of sources.
      */
      void computeHead();

      /**
      * One step of integration loop, from i to i+1.
      *
      * \param iStep time step index, in range 0 to ns.
      */
      void step(int iStep);
 
   private:
      
      DArray<QField> qFields_;

      // Elements of tridiagonal matrices used in propagation
      DArray<double> dA_;
      DArray<double> dB_;
      DArray<double> uA_;
      DArray<double> uB_;

      // Work vector
      DArray<double> v_;

      // Solver used in Crank-Nicholson algorithm
      TridiagonalSolver solver_;

      /// Monomer statistical segment length.
      double step_;

      /// Minimum value of spatial coordinate.
      double xMin_;

      /// Maximum value of spatial coordinate.
      double xMax_;

      /// Spatial grid step size.
      double dx_;

      /// Contour length step size.
      double ds_;

      /// Number of contour length steps = # grid points - 1.
      int ns_;

      /// Number of spatial grid points.
      int nx_;

      /// Has an specified initial condition (head).
      bool hasHead_;

   };

   /*
   * Return q-field at beginning of block.
   */
   inline Propagator::QField const& Propagator::head() const
   {  return qFields_[0]; }

   /*
   * Return q-field at end of block, after solution.
   */
   inline Propagator::QField const& Propagator::tail() const
   {  return qFields_[ns_-1]; }

} 
#endif

#ifndef RPG_FOURIER_MOVE_TPP
#define RPG_FOURIER_MOVE_TPP

/*
* PSCF - Polymer Self-Consistent Field Theory
*
* Copyright 2016 - 2022, The Regents of the University of Minnesota
* Distributed under the terms of the GNU General Public License.
*/

#include "FourierMove.h"
#include "McMove.h" 
#include <rpg/fts/VecOpFts.h>
#include <rpg/System.h>      
#include <prdc/crystal/shiftToMinimum.h>
#include <prdc/cuda/resources.h>
#include <pscf/mesh/MeshIterator.h>
#include <pscf/math/IntVec.h>
#include <pscf/cuda/CudaRandom.h>
#include <util/random/Random.h>
#include <util/param/ParamComposite.h>
#include <util/global.h>

#include <iostream>
#include <complex>
#include <random>
#include <cmath>

namespace Pscf {
namespace Rpg 
{

   using namespace Util;

   /*
   * Constructor.
   */
   template <int D>
   FourierMove<D>::FourierMove(McSimulator<D>& simulator) 
    : McMove<D>(simulator),
      isAllocated_(false)
   { setClassName("FourierMove"); }

   /*
   * Destructor, empty default implementation.
   */
   template <int D>
   FourierMove<D>::~FourierMove()
   {}

   /*
   * ReadParameters
   */
   template <int D>
   void FourierMove<D>::readParameters(std::istream &in)
   {
      //Read the probability
      readProbability(in);
      // attampt move range [A, -A]
      read(in, "A", stepSize_);
      read(in, "F*", fStar_);
      read(in, "tau", tau_);
      read(in, "N", n_);
      read(in, "volumeFraction", f_);
      read(in, "statisticalSegmentLength", b_);
   }
   
   template <int D>
   void FourierMove<D>::setup()
   {  
      McMove<D>::setup();
      const int nMonomer = system().mixture().nMonomer();
      IntVec<D> const & dimensions = system().mesh().dimensions();
      
      if (!isAllocated_){
         wFieldTmp_.allocate(nMonomer);
         randomFieldR_.allocate(dimensions);
         randomFieldK_.allocate(dimensions);
         wKGrid_.allocate(nMonomer);
         wRGrid_.allocate(nMonomer);
         sqrtSq_.allocate(dimensions);
         for (int i = 0; i < nMonomer; ++i) {
            wFieldTmp_[i].allocate(dimensions);
            wKGrid_[i].allocate(dimensions);
            wRGrid_[i].allocate(dimensions);
         }
         isAllocated_ = true;
      }
      computeRgSquare();
      computeStructureFactor();
   }
   
   /**
   * Compute radius of gyration Rg^2 = Nb^2/6
   */
   template <int D>
   void FourierMove<D>::computeRgSquare()
   {
      rgSquare_ = n_ * b_ * b_ / 6.0;
   }

   /*
   * Debye function g(f,x)
   */
   template <int D>
   double FourierMove<D>::computeDebye(double f, double x)
   {
      double fx = f * x;
      return 2.0 * (fx + exp(-fx) - 1.0) / (x * x);
   }
   
   /*
   * F(x,f) = g(1,x)/{g(f,x)g(1-f,x) - [g(1,x) - g(f,x) - g(1-f,x)]^2/4}
   */
   template <int D>
   double FourierMove<D>::computeF(double x)
   {
      //Only for diblock
      if (f_ > 0.5){
         f_ = 1.0 - f_;
      }
      
      double g1 = computeDebye(1.0, x);
      double gf = computeDebye(f_, x);
      double g1f = computeDebye(1.0 - f_, x);
      double denominator = gf * g1f - pow(g1 - gf- g1f, 2.0)/4.0;
      return g1 / denominator;
   }
   
   /*
   * S(q)/N Fredrickson-Helfand structure factor. S(q)/N = 1/(F(x) - F* + tau_)
   */  
   template <int D>
   double FourierMove<D>::computeS(double qSquare)
   {
      //Input is square magnitude of reciprocal basis vector
      double x = qSquare * rgSquare_;
      double denominator =  computeF(x) - fStar_ + tau_;
      // Nbar = Rho^(-2)*a^6* N && Rho = 1/vMonomer && vMonomer = 1/sqrt(Nbar)
      const double vMonomer = system().mixture().vMonomer();
      return 1.0 / (vMonomer * denominator);
   }
   
   template <int D>
   void FourierMove<D>::computeStructureFactor()
   {
      int meshSize = system().domain().mesh().size();
      HostDArray<cudaReal> sqTemp(meshSize);
      MeshIterator<D> itr;
      itr.setDimensions(system().mesh().dimensions());
      IntVec<D> G; IntVec<D> Gmin; 
      double qSq; double S_; 
      for (itr.begin(); !itr.atEnd(); ++itr){
         if (itr.rank() == 0){
            sqTemp[itr.rank()] = 0;
         } else {
            G = itr.position();
            Gmin = shiftToMinimum(G, system().mesh().dimensions(), system().unitCell());
            qSq = system().unitCell().ksq(Gmin);
            // Compute Fredrickson-Helfand structure factor
            S_ = computeS(qSq);
            sqTemp[itr.rank()] = sqrt(S_);
         }
      }
      sqrtSq_ = sqTemp; // copy sqTemp to device
   }
   
   /*
   * Attempt unconstrained move
   */
   template <int D>
   void FourierMove<D>::attemptMove()
   {
   
      const int nMonomer = system().mixture().nMonomer();
      const int meshSize = system().domain().mesh().size();
      
      for (int i = 0; i < nMonomer; ++i) {
         VecOp::eqV(wRGrid_[i], system().w().rgrid(i));
      }

      // Convert real grid to KGrid format
      for (int i = 0; i < nMonomer; ++i) {
         system().fft().forwardTransform(wRGrid_[i], wKGrid_[i]);
      }
      
      for (int i = 0; i < nMonomer; ++i){
         // Generate random number for real part
         cudaRandom().uniform(randomFieldR_);
         // Generate random numbers between [-stepSize_,stepSize_]
         VecOpFts::mcftsScale(randomFieldR_, stepSize_);
         // Generate random numbers between [-stepSize_*S(q)^(1/2), stepSize_*S(q)^(1/2)]
         VecOp::mulEqV(randomFieldR_, sqrtSq_);
         
         // Generate random number for imaginary part
         cudaRandom().uniform(randomFieldK_);
         // Generate random numbers between [-stepSize_,stepSize_]
         VecOpFts::mcftsScale(randomFieldK_, stepSize_);
         // Generate random numbers between [-stepSize_*S(q)^(1/2), stepSize_*S(q)^(1/2)]
         VecOp::mulEqV(randomFieldK_, sqrtSq_);
   
         // Attempt Move in ourier (k-grid) format
         VecOpFts::fourierMove(wKGrid_[i], randomFieldR_, randomFieldK_);
      }
      
      // Convert back to real space (destroys data in wKGrid_)
      for (int i = 0; i < nMonomer; ++i) {
         system().fft().inverseTransformUnsafe(wKGrid_[i], wFieldTmp_[i]);
      }
      // Update attemptMove
      system().setWRGrid(wFieldTmp_);
   }


   /*
   * Trivial default implementation - do nothing
   */
   template <int D>
   void FourierMove<D>::output()
   {}
   
   template<int D>
   void FourierMove<D>::outputTimers(std::ostream& out)
   {
      // Output timing results, if requested.
      out << "\n";
      out << "FourierMove time contributions:\n";
      McMove<D>::outputTimers(out);
   }

}
}
#endif

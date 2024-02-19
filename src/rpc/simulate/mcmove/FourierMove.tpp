#ifndef RPC_FOURIER_MOVE_TPP
#define RPC_FOURIER_MOVE_TPP

/*
* PSCF - Polymer Self-Consistent Field Theory
*
* Copyright 2016 - 2022, The Regents of the University of Minnesota
* Distributed under the terms of the GNU General Public License.
*/

#include "FourierMove.h"
#include "McMove.h" 
#include <rpc/System.h>      
#include <prdc/crystal/shiftToMinimum.h>
#include <pscf/mesh/MeshIterator.h>
#include <pscf/math/IntVec.h>
//#include <pscf/math/RealVec.h>
#include <util/random/Random.h>
#include <util/param/ParamComposite.h>
#include <util/global.h>

#include <iostream>
#include <complex>
#include <random>
#include <cmath>

//#include <util/archives/Serializable_includes.h>
//#include <util/format/Int.h>

namespace Pscf {
namespace Rpc 
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
         wKGrid_.allocate(nMonomer);
         for (int i = 0; i < nMonomer; ++i) {
            wFieldTmp_[i].allocate(dimensions);
            wKGrid_[i].allocate(dimensions);
         }
         isAllocated_ = true;
      }
      structureFactors_.allocate(system().mesh().size());
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
      MeshIterator<D> itr;
      itr.setDimensions(system().mesh().dimensions());
      IntVec<D> G; IntVec<D> Gmin; 
      double qSq; double S_; 
      for (itr.begin(); !itr.atEnd(); ++itr){
         if (itr.rank() == 0){
            structureFactors_[itr.rank()] = 0;
         }
         // Obtain square magnitude of reciprocal basis vector
         if (itr.rank() >= 1){
            G = itr.position();
            Gmin = shiftToMinimum(G, system().mesh().dimensions(), system().unitCell());
            qSq = system().unitCell().ksq(Gmin);
            // Compute Fredrickson-Helfand structure factor
            S_ = computeS(qSq);
            structureFactors_[itr.rank()] = S_;
         }
      }
   }
   
   
   
   /*
   * Attempt unconstrained move
   */
   template <int D>
   void FourierMove<D>::attemptMove()
   {
   
      const int nMonomer = system().mixture().nMonomer();
      // Convert real grid to KGrid format
      for (int i = 0; i < nMonomer; ++i) {
         system().fft().forwardTransform(system().w().rgrid()[i], wKGrid_[i]);
      }
      // Complex Random number Generator
      // Set up the Mersenne Twister engine  
      std::mt19937_64 engine(std::random_device{}());
      // Set up the uniform distribution for real and imaginary parts
      std::uniform_real_distribution<double> real_dist(-stepSize_, stepSize_);
      std::uniform_real_distribution<double> imag_dist(-stepSize_, stepSize_);
      
      /// Move step size in Fourier space Delta_W(q) is randomly selected from 
      /// uniform distribution [-stepSize_*S(q)^(1/2), stepSize_*S(q)^(1/2)]      
      double S_; 
      for (int i = 0; i < nMonomer; i++){
         MeshIterator<D> itr;
         itr.setDimensions(wKGrid_[0].dftDimensions());
         for (itr.begin(); !itr.atEnd(); ++itr){
            // Obtain square magnitude of reciprocal basis vector
            if (itr.rank() >= 1){
               // Compute Fredrickson-Helfand structure factor
               S_ = structureFactors_[itr.rank()];
               // Generate random complex number
               std::complex<double> z(real_dist(engine), imag_dist(engine));
               // Attempt Move in ourier (k-grid) format
               wKGrid_[i][itr.rank()][0] += z.real() * sqrt(S_); 
               wKGrid_[i][itr.rank()][1] += z.imag() * sqrt(S_);
            }
         }
      }
      
      // Convert Fourier (k-grid) format back to Real grid format
      for (int i = 0; i < nMonomer; ++i) {
         system().fft().inverseTransform(wKGrid_[i], wFieldTmp_[i]);
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
      out << "Fourier Move times contributions:\n";
      McMove<D>::outputTimers(out);
   }

}
}
#endif

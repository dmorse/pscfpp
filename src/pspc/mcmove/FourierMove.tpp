#ifndef PSPC_FOURIER_MOVE_TPP
#define PSPC_FOURIER_MOVE_TPP

/*
* PSCF - Polymer Self-Consistent Field Theory
*
* Copyright 2016 - 2022, The Regents of the University of Minnesota
* Distributed under the terms of the GNU General Public License.
*/

#include "FourierMove.h"
#include "McMove.h" 
#include <util/param/ParamComposite.h>
#include <pspc/System.h>      
#include <util/archives/Serializable_includes.h>
#include <util/random/Random.h>
#include <pscf/math/IntVec.h>
#include <pscf/math/RealVec.h>
#include <util/format/Int.h>
#include <pscf/mesh/MeshIterator.h>
#include <iostream>
#include <complex>
#include <random>





namespace Pscf {
namespace Pspc {

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
   * ReadParameters, empty default implementation.
   */
   template <int D>
   void FourierMove<D>::readParameters(std::istream &in)
   {
      //Read the probability
      readProbability(in);
      // attampt move range [A, -A]
      read(in, "A", A_);
      read(in, "F*", Fstar_);
      read(in, "epsilon", epsilon_);
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
      computeRg();
   }
   
   /**
   * Compute radius of gyration Rg = Nb^2/6
   */
   template <int D>
   void FourierMove<D>::computeRg()
   {
      int np = system().mixture().nPolymer();
      if (np > 0) {
         double nb; // number of monomer in each polymer
         double length; // fraction of monomer
         double kuhn; // statistical segment length of monomer
         int i; // molecule index
         int j; // block index
         for (i = 0; i < np; ++i) {
            nb = system().mixture().polymer(i).nBlock();
            Rg_ = 0;
            for (j = 0; j < nb; ++j) {
               Block<D>& block = system().mixture().polymer(i).block(j);
               kuhn = block.kuhn();
               length = block.length();
               Rg_ += length * kuhn * kuhn / 6;
            }
         }
      }      
   }

   /*
   * Debye function g(f,x)
   */
   template <int D>
   double FourierMove<D>::g(double f, double x)
   {
      double fx = f * x;
      return 2.0 * (fx + exp(-fx) - 1.0) / (x * x);
   }
   
   /*
   * F(x,f) = g(1,x)/{g(f,x)g(1-f,x) - [g(1,x) - g(f,x) - g(1-f,x)]^2/4}
   */
   template <int D>
   double FourierMove<D>::F(double x)
   {
      //Only for diblock
      double f = system().mixture().polymer(0).block(0).length();
      if (f > 0.5){
         f = 1.0 - f;
      }
      
      double g1 = g(1.0, x);
      double gf = g(f, x);
      double g1f = g(1.0 - f, x);
      double denominator = gf * g1f - pow(g1 - gf- g1f, 2.0)/4.0;
      return g1 / denominator;
   }
   
   /*
   * S(q)/N Fredrickson-Helfand structure factor. S(q)/N = 1/(F(x) - F* + epsilon_)
   */  
   template <int D>
   double FourierMove<D>::S(double qSquare)
   {
      //Input is square magnitude of reciprocal basis vector
      double x = qSquare * Rg_*Rg_;
      double denominator =  F(x) - Fstar_ + epsilon_;
      //Log::file() << "F(x) " << F(x) << "\n";
      return 1.0 / denominator;
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
      std::uniform_real_distribution<double> real_dist(-A_, A_);
      std::uniform_real_distribution<double> imag_dist(-A_, A_);
      
      MeshIterator<D> itr(wKGrid_[0].dftDimensions());
      IntVec<D> waveBz; 
      int Id; double qSq; double S_; 
      /// Move step size in Fourier space \Delta_W(q) is randomly selected from 
      /// uniform distribution [-A*S(q)^(1/2), A*S(q)^(1/2)]      

      for (int i = 0; i < nMonomer; i++){
         for (itr.begin(); !itr.atEnd(); ++itr){
            // Obtain square magnitude of reciprocal basis vector
            if (itr.rank() >= 1){
               Id = system().basis().waveId(itr.position());
               waveBz = system().basis().wave(Id).indicesBz;
               qSq = system().unitCell().ksq(waveBz);
               // Compute Fredrickson-Helfand structure factor
               S_ = S(qSq);
               // Generate random complex number
               std::complex<double> z(real_dist(engine), imag_dist(engine));
               // Attempt Move in ourier (k-grid) format
               wKGrid_[i][itr.rank()][0] += z.real() * S_; 
               wKGrid_[i][itr.rank()][1] += z.imag() * S_;
            #if 0
            if (itr.rank() == 1){
               double x = qSq * Rg_*Rg_;
               Log::file() << "qSq " << qSq << "\n";
               Log::file() << "Rg_ " << Rg_ << "\n";
               Log::file() << "g1 " << g(1,qSq * Rg_*Rg_) << "\n";
               Log::file() << "g " << g(0.5,qSq * Rg_*Rg_) << "\n";
               Log::file() << "F " << F(qSq * Rg_*Rg_) << "\n";
               Log::file() << "S_ " << S_ << "\n";
            }
            #endif
            }
         }
      }
      //Log::file() << " FourierMove" << "\n";
      
      // Convert Fourier (k-grid) format back to Real grid format
      for (int i = 0; i < nMonomer; ++i) {
         system().fft().inverseTransform(wKGrid_[i], 
                                        wFieldTmp_[i]);
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

}
}
#endif

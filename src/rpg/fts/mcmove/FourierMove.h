#ifndef RPG_FOURIER_MOVE_H
#define RPG_FOURIER_MOVE_H

/*
* PSCF - Polymer Self-Consistent Field Theory
*
* Copyright 2016 - 2022, The Regents of the University of Minnesota
* Distributed under the terms of the GNU General Public License.
*/

#include "McMove.h"                          // base class
#include <prdc/cuda/RFieldDft.h> 
#include <prdc/cuda/RField.h>
#include <util/containers/DArray.h>

namespace Pscf {
namespace Rpg
{

   using namespace Util;
   using namespace Pscf::Prdc;
   using namespace Pscf::Prdc::Cuda;

   /**
   * FourierMove is a Monte Carlo move in fourier space
   *
   * \ingroup Rpg_Fts_McMove_Module
   */
   template <int D>
   class FourierMove : public McMove<D>
   {

   public:

      /**
      * Constructor.
      *
      * \param simulator parent McSimulator
      */
      FourierMove(McSimulator<D>& simulator);

      /**
      * Destructor.
      *
      * Empty default implementation.
      */
      ~FourierMove();

      /**
      * Read required parameters from file.
      *
      * \param in input parameter file stream
      */
      void readParameters(std::istream &in);
      
      /**
      * Output statistics for this move (at the end of simulation)
      */
      void output();
      
      /**
      * Setup before the beginning of each simulation run
      */
      void setup();
      
      /**
      * Return fourier move times contributions.
      *
      * \param out  output stream for timer statistics
      */
      void outputTimers(std::ostream& out);
      
      // Inherited public member function
      using McMove<D>::move;
      using McMove<D>::readProbability;
      using McMove<D>::clearTimers;
      using ParamComposite::read;
      using ParamComposite::setClassName;

   protected:
      
      using McMove<D>::system;
      using McMove<D>::random;
      using McMove<D>::cudaRandom;
      
      /**
      * Attempt unconstrained move.
      *
      * This function should convert system w fields from r-grid format to
      * RFieldDft format, modify the system w fields in fourier space, and
      * convert back to r-grid format, as returned by system().w().rgrid(),
      * in order apply an unconstrained attempted move. The compressor will 
      * then be applied in order to restore the density constraint.
      */
      void attemptMove();


   private:
      
      /**
      * Compute square of radius of gyration Rg = (Nb^2/6)
      */
      void computeRgSquare();
      
      /**
      * Compute Leibler's F(x,f) function. 
      * 
      * F(x,f) = g(1,x)/{g(f,x)g(1-f,x) - [g(1,x) - g(f,x) - g(1-f,x)]^2/4}
      *
      * \param x nondimensionalized value of q^2 Rg^2
      */
      double computeF(double x);
      
      /**
      * Compute Debye function g(f,x).
      * 
      * Debye function: g(f,x) = 2[fx + exp(-fx)-1]/x^2
      *
      * \param f fraction of A block
      * \param x nondimensionalized value of q^2 Rg^2
      */
      double computeDebye(double f, double x);
      
      /**
      * Compute Fredrickson-Helfand structure for specific qSquare.
      *
      * For diblock copolymers,Fredrickson-Helfand structure 
      * S(q)/N = 1/(F(x,f) - F* + epsilon ) x = q^2Rg^2. 
      * q is wave vector, Rg is radius of gyration
      */
      double computeS(double qSquare);
      
      /**
      * Compute Fredrickson-Helfand structure factor
      */
      void computeStructureFactor();
      
      /**
      * Input variable, Move step size stepSize_.
      *
      * In Fourier space Delta W(q) is randomly selected from a uniform
      * distribution [-stepSize_*S(q)^(1/2), stepSize_*S(q)^(1/2)].
      */ 
      double stepSize_;
     
      /** 
      * Input variable, F*. 
      *
      * Ex: at f = 0.5, F* = 20.990; at f = 0.45, F* = 21.396.
      */
      double fStar_;
      
      /**
      * Input variable tau_. 
      *
      * User can calculate using Equation (5.2) from reference:
      * "Fluctuation effects in the theory of microphase separation
      * in block copolymers." Fredrickson, G. H., & Helfand, E. 
      * J. Chem. Phys. 87(1), 697-705 (1987).
      */
      double tau_;
      
      /// Input variable, diblock volume fraction.
      double f_;
      
      /// Input variable, degree of polymerization N. 
      double n_;

      /// Input variable, statistical segment length b.
      double b_;
      
      /// Radius of gyration 
      double rgSquare_;
      
      /// square root of the structure factor
      RField<D> sqrtSq_; 
      
      /// Random fields in real grid between [-stepSize_, stepSize_]
      RField<D> randomFieldR_;
      
      /// Random fields in fourier space between [-stepSize_, stepSize_]
      RField<D> randomFieldK_;
      
      /// wField in in Real Space
      DArray< RField<D> > wRGrid_;
      
      /// wField in in Fourier Space
      DArray< RFieldDft<D> > wKGrid_;
      
      /// Local value of wField after attempt McMove.
      DArray< RField<D> > wFieldTmp_;
      
      /// Has the variable been allocated?
      bool isAllocated_;
      
   };
      
   #ifndef RPG_FOURIER_MOVE_TPP
   // Suppress implicit instantiation
   extern template class FourierMove<1>;
   extern template class FourierMove<2>;
   extern template class FourierMove<3>;
   #endif

}
}
#endif

#ifndef PSPC_FOURIER_MOVE_H
#define PSPC_FOURIER_MOVE_H

/*
* PSCF - Polymer Self-Consistent Field Theory
*
* Copyright 2016 - 2022, The Regents of the University of Minnesota
* Distributed under the terms of the GNU General Public License.
*/

#include "McMove.h"                          //base class
#include <util/global.h>
#include <util/containers/DArray.h>
#include <pspc/field/RField.h>
#include <pspc/field/RFieldDft.h> 
#include <util/param/ParamComposite.h>
#include <pscf/math/IntVec.h>
#include <pscf/math/RealVec.h>

namespace Pscf {
namespace Pspc
{

   using namespace Util;

   /**
   * McMove is an abstract base class for Monte Carlo moves.
   *
   * \ingroup Pspc_McMove_Module
   */
   template <int D>
   class FourierMove : public McMove<D>
   {

   public:

      /**
      * Constructor.
      *
      * \param system parent System object.
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
      * Empty default implementation.
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
      
      
      // Inherited public member function
      using McMove<D>::move;
      using McMove<D>::readProbability;
      using ParamComposite::read;
      using ParamComposite::setClassName;


   protected:
      
      using McMove<D>::system;
      using McMove<D>::random;
      /**
      *  Attempt unconstrained move.
      *
      *  This function should modify the system w fields in r-grid
      *  format, as returned by system().w().rgrid(), in order apply
      *  an unconstrained attempted move. The compressor will then be
      *  applied in order to restore the density constraint.
      *
      */
      void attemptMove();


   private:
      
      /// Move step size in Fourier space \Delta_W(q) is randomly selected from 
      /// uniform distribution [-A*S(q)^(1/2), A*S(q)^(1/2)].
      double A_;
      
      /**
      * For diblock copolymers,
      * Fredrickson-Helfand structure S(q)/N = 1/(F(x,f) - F* + epsilon )
      * x = q^2Rg^2. q is wave vector, Rg is radius of gyration
      */
      
      /// Radius of gyration 
      double Rg_;
      
      /// Input variable F*. Ex: at f = 0.5, F* = 20.990; at f = 0.45, F* = 21.396.
      double Fstar_;
      
      /// Input variable epsilon_.
      double epsilon_;
      
      /**
      * Compute radius of gyration Rg = Nb^2/6
      */
      void computeRg();
      
      /**
      * F(x,f) function. F(x,f) = g(1,x)/{g(f,x)g(1-f,x) - [g(1,x) - g(f,x) - g(1-f,x)]^2/4}
      */
      double F(double x);
      
      /**
      * g(x,f) Debye function. g(f,x) = 2[fx + exp(-fx)-1]/x^2
      */
      double g(const double f, double x);
      
      /**
      * S(q)/N Fredrickson-Helfand structure
      */ 
      double S(double q);
      
      /**
      * Get magnitude of derivative of square of reciprocal basis vector.
      */    
      double ksq(IntVec<D> const & k) const;
      
      /// wField in 
      DArray< RFieldDft<D> > wKGrid_;
      
      /// wField after attempt McMove. local variable wFieldTmp_ in Fourier Space used in attemptMove() function
      DArray< RField<D> > wFieldTmp_;
      
      /**
      * Has the variable been allocated?
      */
      bool isAllocated_;
      
   
   };
      

}
}
#endif

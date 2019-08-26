/*
* PSCF - Polymer Self-Consistent Field Theory
*
* Copyright 2016, The Regents of the University of Minnesota
* Distributed under the terms of the GNU General Public License.
*/

#include "ChiInteraction.h"
#include <pscf/math/LuSolver.h>

namespace Pscf {
   
   using namespace Util;

   /*
   * Constructor.
   */
   ChiInteraction::ChiInteraction()
    : Interaction()
   {  setClassName("ChiInteraction"); }

   /*
   * Destructor.
   */
   ChiInteraction::~ChiInteraction()
   {}

   /*
   * Read chi matrix from file. 
   */
   void ChiInteraction::readParameters(std::istream& in)
   {
      UTIL_CHECK(nMonomer() > 0);
      chi_.allocate(nMonomer(), nMonomer());
      chiInverse_.allocate(nMonomer(), nMonomer());
      indemp_.allocate(nMonomer(), nMonomer());
      readDSymmMatrix(in, "chi", chi_, nMonomer());

      if (nMonomer() == 2) {
         double det = chi_(0,0)*chi_(1, 1) - chi_(0,1)*chi_(1,0);
         double norm = chi_(0,0)*chi_(0, 0) + chi_(1,1)*chi_(1,1) 
                     + 2.0*chi_(0,1)*chi_(1,0);
         if (fabs(det/norm) < 1.0E-8) {
            UTIL_THROW("Singular chi matrix");
         }
         chiInverse_(0,1) = -chi_(0,1)/det;
         chiInverse_(1,0) = -chi_(1,0)/det;
         chiInverse_(1,1) = chi_(0,0)/det;
         chiInverse_(0,0) = chi_(1,1)/det;

      }

      else {

      LuSolver solver;
      solver.allocate(nMonomer());
      solver.computeLU(chi_); 
      solver.inverse(chiInverse_); 

      }

     /* int i, j;

      for (i = 0; i < chiInverse_.capacity1(); ++i) { //row
         for (j = 0; j < chiInverse_.capacity2(); ++j) { //coloumn
           
            std::cout <<chiInverse_(i,j)<<std::endl;

         }
      }*/


        //    std::cout <<chiInverse_.capacity1()<<std::endl;
      //std::cout <<chiInverse_.capacity2()<<std::endl;

      double sum = 0;
      int i, j, k;
      std::cout<<std::endl;

      for (i = 0; i < nMonomer(); ++i) { //coloumn
         indemp_(0,i) = 0;
         for (j = 0; j < nMonomer(); ++j) { //row     
            indemp_(0,i) -= chiInverse_(j,i);        
         }

     //std::cout <<"indemp (0,"<<i<<") = "<<indemp_(0,i)<<std::endl; 
         sum -= indemp_(0,i);
         for (k = 0; k < nMonomer(); ++k) { //row                                             
            indemp_(k,i) = indemp_(0,i);
         }
      }

     // std :: cout <<"sum = "<<sum<<std::endl;

      for (i = 0; i < nMonomer(); ++i) { //row
         for (j = 0; j < nMonomer(); ++j) { //coloumn                                             
            indemp_(i,j) /= sum;
         }
         indemp_(i,i) +=1 ; 
      }


      //for (i = 0; i < chiInverse_.capacity1(); ++i) { //row
        // for (j = 0; j < chiInverse_.capacity2(); ++j) { //coloumn
            
          //  std::cout <<"chiInverse ("<<i<<","<<j<<") = "<<chiInverse_(i,j)<<std::endl;

           // std::cout <<"indemp ("<<i<<","<<j<<") = "<<indemp_(i,j)<<std::endl;
        
        // }
     // }




   }

   /*
   * Compute and return excess Helmholtz free energy per monomer.
   */
   double ChiInteraction::fHelmholtz(Array<double> const & c) const
   {
      int i, j;
      double sum = 0.0;
      for (i = 0; i < nMonomer(); ++i) {
         for (j = 0; j < nMonomer(); ++j) {
            sum += chi_(i, j)* c[i]*c[j];
         }
      }
      return 0.5*sum;
   }

   /*
   * Compute chemical potential from monomer concentrations
   */
   void 
   ChiInteraction::computeW(Array<double> const & c, 
                            Array<double>& w) const
   {
      int i, j;
      for (i = 0; i < nMonomer(); ++i) {
         w[i] = 0.0;
         for (j = 0; j < nMonomer(); ++j) {
            w[i] += chi_(i, j)* c[j];
         }
      }
   }

   /*
   * Compute concentrations and xi from chemical potentials.
   */ 
   void 
   ChiInteraction::computeC(Array<double> const & w, 
                            Array<double>& c, double& xi)
   const
   {
      double sum1 = 0.0;
      double sum2 = 0.0;
      int i, j;
      for (i = 0; i < nMonomer(); ++i) {
         for (j = 0; j < nMonomer(); ++j) {
            sum1 += chiInverse_(i, j)*w[j];
            sum2 += chiInverse_(i, j);
         }
      }
      xi = (sum1 - 1.0)/sum2;
      for (i = 0; i < nMonomer(); ++i) {
         c[i] = 0;
         for (j = 0; j < nMonomer(); ++j) {
            c[i] += chiInverse_(i, j)*( w[j] - xi );
         }
      }
   }

   /*
   * Compute Langrange multiplier from chemical potentials.
   */ 
   void 
   ChiInteraction::computeXi(Array<double> const & w, double& xi)
   const
   {
      double sum1 = 0.0;
      double sum2 = 0.0;
      int i, j;
      for (i = 0; i < nMonomer(); ++i) {
         for (j = 0; j < nMonomer(); ++j) {
            sum1 += chiInverse_(i, j)*w[j];
            sum2 += chiInverse_(i, j);
         }
      }
      xi = (sum1 - 1.0)/sum2;
   }

   /*
   * Return dWdC = chi matrix.
   */
   void 
   ChiInteraction::computeDwDc(Array<double> const & c, Matrix<double>& dWdC)
   const
   {
      int i, j;
      for (i = 0; i < nMonomer(); ++i) {
         for (j = 0; j < nMonomer(); ++j) {
            dWdC(i, j) = chi_(i, j);
         }
      }
   }

} // namespace Pscf

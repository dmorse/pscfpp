/*
* PSCF - Polymer Self-Consistent Field
*
* Copyright 2015 - 2026, The Regents of the University of Minnesota
* Distributed under the terms of the GNU General Public License.
*/

#include "FhInteraction.h"
#include <pscf/interaction/Interaction.h>
#include <pscf/math/LuSolver.h>

namespace Pscf {

   using namespace Util;

   /*
   * Constructor.
   */
   FhInteraction::FhInteraction()
    : nMonomer_(0)
   {  setClassName("FhInteraction"); }

   /*
   * Copy constructor.
   */
   FhInteraction::FhInteraction(FhInteraction const & other)
    : nMonomer_(0)
   {  
      setClassName("FhInteraction"); 
      if (other.nMonomer() > 0) {
         setNMonomer(other.nMonomer());
         setChi(other.chi());
      }
   }

   /*
   * Constructor, copy from from Pscf::Interaction
   */
   FhInteraction::FhInteraction(Interaction const & other)
    : nMonomer_(0)
   {  
      UTIL_CHECK(other.nMonomer() > 0);
      setClassName("FhInteraction"); 
      setNMonomer(other.nMonomer());
      setChi(other.chi());
   }

   /*
   * Destructor.
   */
   FhInteraction::~FhInteraction()
   {}

   /*
   * Assignment from FhInteraction
   */
   FhInteraction& 
   FhInteraction::operator = (FhInteraction const & other)
   {  
      if (this != &other) {
         if (other.nMonomer() == 0) {
            UTIL_CHECK(nMonomer_ == 0);
         } else {
            if (nMonomer_ == 0) {
               setNMonomer(other.nMonomer());
            }
            UTIL_CHECK(nMonomer_ == other.nMonomer());
            setChi(other.chi());
         }
      }
      return *this;
   }

   /*
   * Assignment from Interaction
   */
   FhInteraction& FhInteraction::operator = (Interaction const & other)
   {  
      UTIL_CHECK(other.nMonomer() > 0);
      if (nMonomer_ == 0) {
         setNMonomer(other.nMonomer());
      }
      UTIL_CHECK(nMonomer_ == other.nMonomer());
      setChi(other.chi());
      return *this;
   }

   /*
   * Set the number of monomer types.
   */
   void FhInteraction::setNMonomer(int nMonomer)
   {  
      UTIL_CHECK(nMonomer_ == 0);
      UTIL_CHECK(nMonomer > 0);
      nMonomer_ = nMonomer; 
      chi_.allocate(nMonomer, nMonomer);
      chiInverse_.allocate(nMonomer, nMonomer);
      setChiZero();
   }

   /*
   * Read chi matrix from file.
   */
   void FhInteraction::readParameters(std::istream& in)
   {
      UTIL_CHECK(nMonomer() > 0);
      readDSymmMatrix(in, "chi", chi_, nMonomer());

      // Compute relevant AM iterator quantities.
      updateMembers();
   }

   void FhInteraction::updateMembers()
   {
      UTIL_CHECK(chi_.isAllocated());
      UTIL_CHECK(chiInverse_.isAllocated());

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

      } else {
         LuSolver solver;
         solver.allocate(nMonomer());
         solver.computeLU(chi_);
         solver.inverse(chiInverse_);
      }

   }

   /*
   * Set values for all chi parameters.
   */
   void FhInteraction::setChi(Matrix<double> const & other)
   {
      UTIL_CHECK(nMonomer_ > 0);
      UTIL_CHECK(other.capacity1() == nMonomer_);
      UTIL_CHECK(other.capacity2() == nMonomer_);
      double value1, value2, diff;
      int i, j;
      for (i = 0; i < nMonomer_; ++i) {
         chi_(i, i) = other(i, i);
      }
      if (nMonomer_ > 1) {
         for (i = 0; i < nMonomer_; ++i) {
            for (j = 0; j < i; ++j) {
               value1 = other(i, j);
               value2 = other(j, i);
               diff = std::abs( value1 - value2 );
	       if (diff > 1.0E-10) {
                  Log::file() << "Diff     = " << diff << "\n";
                  Log::file() << " i, j    = " << i << "  " << j << "\n";
                  Log::file() << "chi(i,j) = " << value1 << "\n";
                  Log::file() << "chi(j,i) = " << value2 << "\n";
                  UTIL_THROW("Error: Asymmetric chi matrix");
               }
               chi_(i, j) = value1;
               chi_(j, i) = value1;
            }
         }
      }
      updateMembers();
   }

   /*
   * Set a single chi parameter.
   */
   void FhInteraction::setChi(int i, int j, double chi)
   {  
      UTIL_CHECK(nMonomer_ > 0);
      UTIL_CHECK(i >= 0);
      UTIL_CHECK(i < nMonomer_);
      UTIL_CHECK(j >= 0);
      UTIL_CHECK(j < nMonomer_);
      chi_(i,j) =  chi; 
      if (i != j) {
         chi_(j,i) = chi;
      }

      // Compute relevant AM iterator quantities. 
      updateMembers();
   }

   /*
   * Set all elements of the chi matrix to zero.
   */
   void FhInteraction::setChiZero()
   {
      UTIL_CHECK(nMonomer_ > 0);
      int i, j;
      for (i = 0; i < nMonomer_; i++) {
         for (j = 0; j < nMonomer_; j++) {
            chi_(i, j) = 0.0;
            chiInverse_(i, j) = 0.0;
         }
      }
   }

   /*
   * Compute and return excess Helmholtz free energy per monomer.
   */
   double FhInteraction::fHelmholtz(Array<double> const & c) const
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
   FhInteraction::computeW(Array<double> const & c,
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
   FhInteraction::computeC(Array<double> const & w,
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
   FhInteraction::computeXi(Array<double> const & w, double& xi)
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
   FhInteraction::computeDwDc(Array<double> const & c, 
                               Matrix<double>& dWdC) const
   {
      int i, j;
      for (i = 0; i < nMonomer(); ++i) {
         for (j = 0; j < nMonomer(); ++j) {
            dWdC(i, j) = chi_(i, j);
         }
      }
   }

} // namespace Pscf

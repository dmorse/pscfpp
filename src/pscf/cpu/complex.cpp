/*
* PSCF - Polymer Self-Consistent Field 
*
* Copyright 2015 - 2025, The Regents of the University of Minnesota
* Distributed under the terms of the GNU General Public License.
*/

#include "complex.h"

namespace Pscf {

   // Exponentiation and logarithm

   /*
   * Exponentation of a ffts_complex variable, z = exp(a).
   */
   void assignExp(fftw_complex & z, fftw_complex const & a)
   {
      std::complex<double> arg = std::complex<double>(a[0], a[1]); 
      std::complex<double> result = std::exp(arg);
      z[0] = result.real();
      z[1] = result.imag();
   }

   /*
   * Logarithm of an fftw_complex variable, z = log(a).
   */
   void assignLog(fftw_complex & z, fftw_complex const & a)
   {  
      std::complex<double> arg = std::complex<double>(a[0], a[1]); 
      std::complex<double> result = std::log(arg);
      z[0] = result.real();
      z[1] = result.imag();
   }

   /*
   * Stream extraction operator for fftw_complex.
   */
   std::istream& operator >> (std::istream& is, fftw_complex& z)
   {
      is >> z[0] >> z[1];
      return is;
   }

   /*
   * Stream insertion operator for fftw_complex.
   */
   std::ostream& operator << (std::ostream& os, fftw_complex const & z)
   {
      os <<  z[0] << "  " << z[1];
      return os;
   }

}

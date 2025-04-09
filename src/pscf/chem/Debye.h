#ifndef PSCF_DEBYE_H
#define PSCF_DEBYE_H

/*
* PSCF - Polymer Self-Consistent Field Theory
*
* Copyright 2016 - 2025, The Regents of the University of Minnesota
* Distributed under the terms of the GNU General Public License.
*/

namespace Pscf {

/**
* Functions used to compute intramolecular correlation functions.
* 
* \ingroup Pscf_Chem_Module
*/
namespace Debye {

   /**
   * Compute and return a homopolymer correlation function.
   * 
   * This function returns the intramolecular correlation function
   * for a homopolymer of specified length and statistical segment
   * length. The result for the thread model can be expressed as a
   * function 
   * \f[
   *   d = (length)^2 g(x)
   * \f] 
   * where x =  ksq * length * kuhn * kuhn / 6 and
   * \f[
   *     g(x) \equiv 2[ e^{-x} - 1 + x ]/x^2
   * \f]
   * is the Debye function.
   *
   * \ingroup Pscf_Chem_Module
   *
   * \param ksq  square of wavenumber
   * \param length  contour length of polymer or block
   * \param kuhn  statistical segement length
   */
   double d(double ksq, double length, double kuhn);

   /**
   * Compute and return one-sided correlation factor for one block.
   *
   * This function returns the function
   * \f[
   *   e = (length) h(x)
   * \f]
   * where x =  ksq * length * kuhn * kuhn / 6, and
   * \f[
   *     h(x) \equiv [ 1 - e^{-x} ]/x 
   * \f]
   * is another dimensionless function of x.
   *
   * \ingroup Pscf_Chem_Module
   *
   * \param ksq  square of wavenumber
   * \param length  contour length of block
   * \param kuhn  statistical segement length
   */
   double e(double ksq, double length, double kuhn);

} // namespace Debye
} // namespace Pscf
#endif

#ifndef PSCF_DEBYE_H
#define PSCF_DEBYE_H

/*
* PSCF - Polymer Self-Consistent Field Theory
*
* Copyright 2015 - 2025, The Regents of the University of Minnesota
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
   * Compute and return intrablock correlation function (thread model)
   * 
   * This function returns the intramolecular correlation function for a
   * homopolymer of specified length and statistical segment length. The
   * The result for the thread model can be expressed as a function 
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
   double dt(double ksq, double length, double kuhn);

   /**
   * Compute and return an intrablock correlation function (bead model)
   * 
   * This function returns the intramolecular correlation function for 
   * a homopolymer of specified length and statistical segment length. 
   * The result for the bead model can be expressed as a function 
   * \f[
   *     g(x) \equiv 2[ e^{-yN} - 1 + N(1-e^{-y}) ]/(1-e^{-y})^2
   * \f]
   * where y =  ksq * kuhn * kuhn / 6 and
   *
   * \ingroup Pscf_Chem_Module
   *
   * \param ksq  square of wavenumber
   * \param nBead  number of beads in the block (converted to double)
   * \param kuhn  statistical segement length
   */
   double db(double ksq, double nBead, double kuhn);

   /**
   * Compute and return one-sided factor for one block (thread model).
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
   double et(double ksq, double length, double kuhn);

   /**
   * Compute and return one-sided factor for one block (thread model).
   *
   * This function returns the function
   * \f[
   *     e(x) \equiv ( 1 - e^{-Ny} )/(1-e^{-y})
   * \f]
   * where y =  ksq * kuhn * kuhn / 6 and
   *
   * \ingroup Pscf_Chem_Module
   *
   * \param ksq  square of wavenumber
   * \param nBead  number of beads in the block (converted to double)
   * \param kuhn  statistical segement length
   */
   double eb(double ksq, double nBead, double kuhn);

} // namespace Debye
} // namespace Pscf
#endif

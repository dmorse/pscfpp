#ifndef PSCF_CORRELATION_DEBYE_H
#define PSCF_CORRELATION_DEBYE_H

/*
* PSCF - Polymer Self-Consistent Field
*
* Copyright 2015 - 2026, The Regents of the University of Minnesota
* Distributed under the terms of the GNU General Public License.
*/

namespace Pscf {
namespace Correlation {

   /**
   * Compute and return intrablock correlation function (thread model)
   * 
   * This function returns the intramolecular correlation function for a
   * homopolymer of specified length and statistical segment length. The
   * result for the thread model can be expressed as a function 
   * \f[
   *   d(k) = L^2 g(x)
   * \f] 
   * where \f$ x =  k^2 L b^2 / 6 \f$, for L = length, b = kuhn, and
   * \f$ k^{2} \f$ = ksq, and
   * \f[
   *     g(x) \equiv 2[ e^{-x} - 1 + x ]/x^2
   * \f]
   * is the Debye function. This function also gives the intra-block
   * correlation function for a Gaussian block of contour given by the 
   * parameter "length" and statistical segment length given by "kuhn".
   *
   * \ingroup Pscf_Correlation_Module
   *
   * \param ksq  square of wavenumber
   * \param length  contour length of polymer or block
   * \param kuhn  statistical segement length
   */
   double dt(double ksq, double length, double kuhn);

   /**
   * Compute and return an intrablock correlation function (bead model)
   * 
   * This function returns the intramolecular correlation function for a
   * homopolymer of specified length and statistical segment length. The
   * result for the bead model can be expressed as a function 
   * \f[
   *     d(k) \equiv 2[ e^{-yN} - 1 + N(1-e^{-y}) ]/(1-e^{-y})^2
   * \f]
   * where \f$ y = k^2 b^2 / 6 \f$, for N = nBead, b = kuhn,  and 
   * \f$ k^{2} = \f$ ksq.  This function also gives the intra-block 
   * correlation function for block of N beads.
   *
   * \ingroup Pscf_Correlation_Module
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
   *   e(k) = L [ 1 - e^{-x} ] / x
   * \f]
   * where \f$ x =  k^2 L b^2 / 6 \f$, for L = length, b = kuhn, and
   * \f$ k^{2} \f$ = ksq.
   *
   * The intra-block correlation function \f$ \omega_{ij}(k) \f$ for two
   * distinct blocks with block indices i and j can be expressed in the
   * thread model as a product
   * \f[
   *    \omega_{ij}(k) = e^{-k^{2}R_{ij}^{2}/6} e_{i}(k) e_{j}(k)
   * \f]
   * where \f$ R_{ij}^{2} \f$ is the mean-squared end-to-end length of a
   * sequence of other blocks that form a path connecting the two blocks 
   * of interest (if any), and where \f$ e_{i}(k) \f$ is the one sided 
   * factor returned by this function for block \f$ i \f$. In the thread
   * model \f$ R_{ij}^{2} = 0 \f$ for blocks that both terminate at a 
   * shared vertex.
   *
   * \ingroup Pscf_Correlation_Module
   *
   * \param ksq  square of wavenumber
   * \param length  contour length of block
   * \param kuhn  statistical segement length
   */
   double et(double ksq, double length, double kuhn);

   /**
   * Compute and return one-sided factor for one block (bead model).
   *
   * This function returns the function
   * \f[
   *     e(x) \equiv ( 1 - e^{-Ny} )/(1-e^{-y})
   * \f]
   * where \f$ y =  k^2 b^2 / 6 \f$, for N = nBead, b = kuhn,  and 
   * \f$ k^{2} \f$ = ksq.
   *
   * The intramolecular correlation function \f$ \omega_{ij}(k) \f$ for two
   * distinct blocks with block indices i and j can be expressed in the
   * bead model as a product
   * \f[
   *    \omega_{ij}(k) = e^{-k^{2}R_{ij}^{2}/6} e_{i}(k) e_{j}(k)
   * \f]
   * where \f$ R_{ij}^{2} \f$ is the mean-squared end-to-end length of the
   * sequence of bonds that lie along a path connecting the two blocks, 
   * and where \f$ e_{i}(k) \f$ is the one sided factor returned by this 
   * function for block \f$ i \f$. In the bead model, the value of 
   * \f$ R_{ij}^{2} \f$ is a sum of squares of statistical segment lengths 
   * for all bonds along the path that connects the end monomers of the 
   * two blocks. In the bead model, if blocks \f$ i \f$ and \f$ j \f$
   * terminate at shared vertex, and have statistical segment lengths 
   * \f$ b_{i} \f$ and \f$ b_{j} \f$, then they are taken to be connected 
   * by a connecting bond with a effective squared statistical segment 
   * length \f$ R^{2}_{ij} = ( b_{i}^{2} + b_{j}^{2} )/2 \f$.
   *
   * \ingroup Pscf_Correlation_Module
   *
   * \param ksq  square of wavenumber
   * \param nBead  number of beads in the block (converted to double)
   * \param kuhn  statistical segement length within the block
   */
   double eb(double ksq, double nBead, double kuhn);

} // namespace Correlation
} // namespace Pscf
#endif

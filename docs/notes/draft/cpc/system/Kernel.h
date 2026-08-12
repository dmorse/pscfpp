#ifndef CP_KERNEL_H
#define CP_KERNEL_H

/*
* PSCF - Polymer Self-Consistent Field
*
* Copyright 2015 - 2025, The Regents of the University of Minnesota
* Distributed under the terms of the GNU General Public License.
*/

#include <util/param/ParamComposite.h>    // base class

namespace Pscf {
namespace Cp {

   using namespace Util;

   /**
   * Kernel for interactions of nonzero range.
   *
   * \ingroup Cp_System_Module
   */
   class Kernel : public ParamComposite
   {

   public:

      /**
      * Constructor.
      */
      Kernel();

      /**
      * Destructor.
      */
      virtual ~Kernel();

      // Modifiers

      /**
      * Read model parameters.
      *
      * \pre Must be called after setNMonomer.
      */
      virtual void readParameters(std::istream& in);

      // Accessors

      /** 
      * Return the range of binary interactions.
      */
      double range() const; 

      /**
      * Compute the Fourier-space kernel for particle smearing.
      *
      * \param kSq  squared wavenumber (input)
      * \return Fourier space kernel for particle smearing
      */
      double f(double kSq) const;

      /**
      * Compute the Fourier-space kernel for pair interactions.
      * 
      * This must be defined such that f = g*g for fixed wavenumber.
      *
      * \param kSq  squared wavenumber (input)
      * \return damping Fourier-space factor for pair interactions
      */
      double g(double kSq) const;

   private:

      // Range of interaction
      double range_;

      // Constant prefactor used in calculation of g
      double cg_;

      // Constant prefactor used in calculation of f
      double cf_;

   };

   // Inline function

   inline double Kernel::range() const
   {  return range_; }

} // namespace Cp
} // namespace Pscf
#endif

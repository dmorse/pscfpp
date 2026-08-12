/*
* PSCF - Polymer Self-Consistent Field
*
* Copyright 2015 - 2025, The Regents of the University of Minnesota
* Distributed under the terms of the GNU General Public License.
*/

#include "Kernel.h"
#include <cmath>

namespace Pscf {
namespace Cp {

   using namespace Util;

   /*
   * Constructor.
   */
   Kernel::Kernel()
    : range_(1.0),
      cg_(0.25),
      cf_(0.5)
   {  setClassName("Kernel"); }

   /*
   * Destructor.
   */
   Kernel::~Kernel()
   {}

   /*
   * Read parameters from file.
   */
   void Kernel::readParameters(std::istream& in)
   {
      read(in, "range", range_);

      cg_ = -0.25/(range_ * range_);
      cf_ = 2.0*cg_;
   }

   /*
   * Compute & return Fourier-space kernel for pair interactions.
   */
   double Kernel::f(double kSq) const
   {
      return std::exp(cf_ * kSq);
   }

   /*
   * Compute & return Fourier-space kernel for particle smearing.
   */
   double Kernel::g(double kSq) const
   {
      return std::exp(cg_ * kSq);
   }

} // namespace Cp
} // namespace Pscf

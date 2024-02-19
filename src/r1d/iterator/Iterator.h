#ifndef R1D_ITERATOR_H
#define R1D_ITERATOR_H

/*
* PSCF - Polymer Self-Consistent Field Theory
*
* Copyright 2016 - 2022, The Regents of the University of Minnesota
* Distributed under the terms of the GNU General Public License.
*/

#include <util/param/ParamComposite.h>    // base class
#include <r1d/SystemAccess.h>            // base class
#include <util/global.h>                  

namespace Pscf {
namespace R1d
{

   using namespace Util;

   /**
   * Base class for iterative solvers for SCF equations.
   *
   * \ingroup R1d_Iterator_Module
   */
   class Iterator : public ParamComposite, public SystemAccess
   {

   public:

      /**
      * Default constructor.
      */
      Iterator();

      /**
      * Constructor.
      * 
      * \param system parent System object
      */
      Iterator(System& system);

      /**
      * Destructor.
      */
      ~Iterator();

      /**
      * Iterate to solution.
      *
      * \param isContinuation true iff part of sweep, and not first step.
      * \return error code: 0 for success, 1 for failure.
      */
      virtual int solve(bool isContinuation = false) = 0;

   };

} // namespace R1d
} // namespace Pscf
#endif

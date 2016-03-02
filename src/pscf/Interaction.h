#ifndef PSCF_INTERACTION_H
#define PSCF_INTERACTION_H

/*
* PSCF - Polymer Self-Consistent Field Theory
*
* Copyright 2016, The Regents of the University of Minnesota
* Distributed under the terms of the GNU General Public License.
*/

#include <util/param/ParamComposite.h>    // base class
#include <util/containers/Array.h>        // argument
#include <util/global.h>                  

namespace Pscf {

   using namespace Util;

   /**
   * Base class for excess free energy models.
   *
   * \ingroup Pscf_Fd1d_Module
   */
   class Interaction : public ParamComposite
   {

   public:

      /**
      * Constructor.
      */
      Interaction();

      /**
      * Destructor.
      */
      virtual ~Interaction();

      /**
      * Set the number of monomer types.
      *
      * \param nMonomer number of monomer types.
      */
      void setNMonomer(int nMonomer);

      /**
      * Compute wField from and pressure p.
      *
      * \return error code: 0 for success, 1 for failure.
      */
      virtual 
      void computeW(Array<double> const & cField, double p, 
                    Array<double>& wField)
      {}

      /**
      * Compute wField from and pressure p.
      *
      * \return error code: 0 for success, 1 for failure.
      */
      virtual 
      void computeC(Array<double> const & wField, 
                    Array<double>& cField)
      {}

      /**
      * Get number of monomer types.
      */
      int nMonomer() const;

   private:

      /// Number of monomers
      int nMonomer_;

   };

   /**
   * Get number of monomer types.
   */
   inline 
   int Interaction::nMonomer() const
   {  return nMonomer_; }

} // namespace Pscf
#endif

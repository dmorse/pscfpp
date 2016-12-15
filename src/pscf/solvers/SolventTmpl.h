#ifndef PSCF_SOLVENT_TMPL_H
#define PSCF_SOLVENT_TMPL_H

/*
* PSCF - Polymer Self-Consistent Field Theory
*
* Copyright 2016, The Regents of the University of Minnesota
* Distributed under the terms of the GNU General Public License.
*/

#include <pscf/chem/Species.h>
#include <util/param/ParamComposite.h>

namespace Pscf
{ 

   using namespace Util;

   /**
   * Template for a class representing a solvent species.
   *
   * Template argument TP is a propagator class. This is 
   * only used to define the data types for concentration 
   * and chemical potential fields.
   *
   * \ingroup Pscf_Solvers_Module
   */
   template <class TP>
   class SolventTmpl : public Species, public ParamComposite
   {
   public:

      /**
      * Monomer concentration field.
      */
      typedef typename TP::CField CField;

      /** 
      * Monomer chemical potential field.
      */
      typedef typename TP::WField WField;

      /**
      * Constructor.
      */
      SolventTmpl()
      {}
   
      /**
      * Constructor.
      */
      ~SolventTmpl()
      {}
   
      /**
      * Compute monomer concentration field and phi and/or mu.
      *
      * Pure virtual function: Must be implemented by subclasses.
      * Upon return, concentration field, phi and mu are all set.
      *
      * \param wField monomer chemical potential field.
      */
      virtual void compute(WField const & wField ) = 0;
   
      /**
      * Get monomer concentration field for this solvent.
      */
      const CField& concentration() const
      {  return concentration_;  }
   
   protected:

      CField concentration_;
   
   };

} 
#endif 

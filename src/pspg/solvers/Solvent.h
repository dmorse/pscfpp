#ifndef PSPG_SOLVENT_H
#define PSPG_SOLVENT_H

/*
* PSCF - Polymer Self-Consistent Field Theory
*
* Copyright 2016, The Regents of the University of Minnesota
* Distributed under the terms of the GNU General Public License.
*/

#include <pscf/chem/Species.h>
#include <util/param/ParamComposite.h>
#include <pspg/field/RDField.h>

namespace Pscf { 
namespace Pspg { 

   using namespace Util;

   /**
   * Class representing a solvent species.
   *
   * \ingroup Pspg_Solvers_Module
   */
   template <int D>
   class Solvent : public Species, public ParamComposite
   {
   public:

      /**
      * Monomer concentration field.
      */
      typedef RDField<D> CField;

      /** 
      * Monomer chemical potential field.
      */
      typedef RDField<D> WField;

      /**
      * Constructor.
      */
      Solvent()
      {}
   
      /**
      * Constructor.
      */
      ~Solvent()
      {}
   
      /**
      * Compute monomer concentration field and phi and/or mu.
      *
      * Pure virtual function: Must be implemented by subclasses.
      * Upon return, concentration field, phi and mu are all set.
      *
      * \param wField monomer chemical potential field.
      */
      virtual void compute(WField const & wField )
      {};
   
      /**
      * Get monomer concentration field for this solvent.
      */
      const CField& concentration() const
      {  return concentration_;  }
   
   private:

      CField concentration_;
   
   };

}
} 
#endif 

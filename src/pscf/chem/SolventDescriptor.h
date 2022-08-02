#ifndef PSCF_SOLVENT_DESCRIPTOR_H
#define PSCF_SOLVENT_DESCRIPTOR_H

/*
* PSCF - Polymer Self-Consistent Field Theory
*
* Copyright 2016 - 2019, The Regents of the University of Minnesota
* Distributed under the terms of the GNU General Public License.
*/

#include <pscf/chem/Species.h>             // base class
#include <util/param/ParamComposite.h>     // base class

namespace Pscf { 

   using namespace Util;

   /**
   * Descriptor for a solvent species.
   *
   * \ingroup Pscf_Chem_Module
   */
   class SolventDescriptor : public Species, public ParamComposite
   {

   public:

      /**
      * Constructor.
      */
      SolventDescriptor();
   
      /**
      * Constructor.
      */
      ~SolventDescriptor();
   
      /**
      * Read and initialize.
      *
      * \param in input parameter stream
      */
      virtual void readParameters(std::istream& in);

      /// \name Setters (set member data)
      ///@{

      /**
      * Set value of phi (volume fraction), if ensemble is closed.
      *
      * \throw Exception if ensemble is open
      * \param phi desired volume fraction for this species
      */
      void setPhi(double phi);

      /**
      * Set value of mu (chemical potential), if ensemble is closed.
      *
      * \throw Exception if ensemble is closed
      * \param mu  desired chemical potential for this species
      */
      void setMu(double mu);

      /**
      * Set the monomer id for this solvent.
      *
      * \param monomerId  integer id of monomer type, in [0,nMonomer-1]
      */ 
      void setMonomerId(int monomerId);
  
      /**
      * Set the size or volume of this solvent species.
      *
      * The ``size" is steric volume / reference volume.
      *
      * \param size volume of solvent
      */ 
      void setSize(double size);

      ///@}
      /// \name Accessors (getters)
      ///@{
 
      /**
      * Get the monomer type id.
      */ 
      int monomerId() const;
  
      /**
      * Get the size (number of monomers) in this solvent.
      */
      double size() const;

      ///@}

      // Inherited accessor functions 
      using Pscf::Species::phi;
      using Pscf::Species::mu;
      using Pscf::Species::q;
      using Pscf::Species::ensemble;

   protected:

      // Inherited data members
      using Pscf::Species::phi_;
      using Pscf::Species::mu_;
      using Pscf::Species::q_;
      using Pscf::Species::ensemble_;

      /// Identifier for the associated monomer type.
      int monomerId_;

      /// Size of this block = volume / monomer reference volume. 
      double size_;

   };
   
   // Inline member functions

   /*
   * Get the monomer type id.
   */ 
   inline int SolventDescriptor::monomerId() const
   {  return monomerId_; }

   /*
   * Get the size (number of monomers) in this block.
   */
   inline double SolventDescriptor::size() const
   {  return size_; }

}
#endif 

#ifndef PSCF_SOLVENT_TMPL_H
#define PSCF_SOLVENT_TMPL_H

/*
* PSCF - Polymer Self-Consistent Field Theory
*
* Copyright 2016 - 2019, The Regents of the University of Minnesota
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
   * \ingroup Pscf_Solver_Module
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
      * Set the monomer id.
      *
      * \param monomerId integer id of monomer type (>=0)
      */ 
      void setMonomerId(int monomerId);
  
      /**
      * Set the size of this block.
      *
      * The ``size" is steric volume / reference volume.
      *
      * \param size block size 
      */ 
      void setSize(double size);
  
      /**
      * Compute monomer concentration field and phi and/or mu.
      *
      * Pure virtual function: Must be implemented by subclasses.
      * Upon return, concentration field, phi and mu are all set.
      *
      * \param wField monomer chemical potential field.
      */
      virtual void compute(WField const & wField ) = 0;
   
      //@}
      /// \name Accessors (getters)
      //@{
    
      /**
      * Get the monomer type id.
      */ 
      int monomerId() const;
  
      /**
      * Get the size (number of monomers) in this block.
      */
      double size() const;

      /**
      * Get monomer concentration field for this solvent.
      */
      const CField& concentration() const;
      //{  return concentration_;  }
   
      //@}

   protected:

      CField concentration_;
   
   private:
  
      /// Identifier for the associated monomer type.
      int monomerId_;

      /// Size of this block = volume / monomer reference volume. 
      double size_;

      #if 0
      friend 
      std::istream& 
      operator >> (std::istream& in, SolventTmpl<TP> &block);

      friend 
      std::ostream& 
      operator << (std::ostream& out, const SolventTmpl<TP> &block);
      #endif

   };

   // Inline member functions

   /*
   * Get the monomer type id.
   */ 
   template <class TP>
   inline int SolventTmpl<TP>::monomerId() const
   {  return monomerId_; }

   /*
   * Get the size (number of monomers) in this block.
   */
   template <class TP>
   inline double SolventTmpl<TP>::size() const
   {  return size_; }
    
   /**
   * Get monomer concentration field for this solvent.
   */
   template <class TP>
   typename SolventTmpl<TP>::CField const & SolventTmpl<TP>::concentration() const
   {  return concentration_;  }
   
   // Non-inline functions

   #if 0
   /**
   * istream extractor for a SolventTmpl<TP>.
   *
   * \param in  input stream
   * \param block  SolventTmpl<TP> to be read from stream
   * \return modified input stream
   */
   template <class TP>
   std::istream& operator >> (std::istream& in, SolventTmpl<TP> &block);

   /**
   * ostream inserter for a SolventTmpl<TP>.
   *
   * \param out  output stream
   * \param block  SolventTmpl<TP> to be written to stream
   * \return modified output stream
   */
   template <class TP>
   std::ostream& 
   operator << (std::ostream& out, const SolventTmpl<TP> &block);

   /*
   * Serialize to/from an archive.
   */
   template <class Archive>
   void SolventTmpl<TP>::serialize(Archive& ar, unsigned int)
   {
      ar & monomerId_;
      ar & size_;
   }
   #endif
    
} 
#endif 

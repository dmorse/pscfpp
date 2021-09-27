#ifndef PSCF_SOLVENT_DESCRIPTOR_H
#define PSCF_SOLVENT_DESCRIPTOR_H

/*
* PSCF - Polymer Self-Consistent Field Theory
*
* Copyright 2016 - 2019, The Regents of the University of Minnesota
* Distributed under the terms of the GNU General Public License.
*/

#include <pscf/chem/Species.h>

namespace Pscf
{ 

   using namespace Util;

   /**
   * Chemical description of a solvent species.
   *
   * \ingroup Pscf_Chem_Module
   */
   class SolventDescriptor : public Species
   {
   public:

      /**
      * Constructor.
      */
      SolventDescriptor()
      {}
   
      /**
      * Constructor.
      */
      ~SolventDescriptor()
      {}
   
      /**
      * Set the monomer id for this solvent.
      *
      * \param monomerId integer id of monomer type (>=0)
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
  
      //@}
      /// \name Accessors (getters)
      //@{
 
      /**
      * Get the monomer type id.
      */ 
      int monomerId() const;
  
      /**
      * Get the size (number of monomers) in this solvent.
      */
      double size() const;

      //@}

   private:
  
      /// Identifier for the associated monomer type.
      int monomerId_;

      /// Size of this block = volume / monomer reference volume. 
      double size_;

      #if 0
      friend 
      std::istream& 
      operator >> (std::istream& in, SolventDescriptor &block);

      friend 
      std::ostream& 
      operator << (std::ostream& out, const SolventDescriptor &block);
      #endif

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
    
   // Non-inline functions

   #if 0
   /**
   * istream extractor for a SolventDescriptor.
   *
   * \param in  input stream
   * \param block  SolventDescriptor to be read from stream
   * \return modified input stream
   */
   std::istream& operator >> (std::istream& in, SolventDescriptor &block);

   /**
   * ostream inserter for a SolventDescriptor.
   *
   * \param out  output stream
   * \param block  SolventDescriptor to be written to stream
   * \return modified output stream
   */
   std::ostream& 
   operator << (std::ostream& out, const SolventDescriptor &block);

   /*
   * Serialize to/from an archive.
   */
   template <class Archive>
   void SolventDescriptor::serialize(Archive& ar, unsigned int)
   {
      ar & monomerId_;
      ar & size_;
   }
   #endif
    
} 
#endif 

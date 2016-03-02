#ifndef PSCF_BLOCK_DESCRIPTOR_H
#define PSCF_BLOCK_DESCRIPTOR_H

/*
* PSCF - Polymer Self-Consistent Field Theory
*
* Copyright 2016, The Regents of the University of Minnesota
* Distributed under the terms of the GNU General Public License.
*/

#include <util/containers/Pair.h>

#include <iostream>

namespace Pscf
{ 

   using namespace Util;

   /**
   * A linear homopolymer block within a block copolymer.
   *
   * \ingroup Pscf_Base_Module
   */
   class BlockDescriptor
   {
   public:
  
      /**
      * Constructor.
      */ 
      BlockDescriptor();
  
      /**
      * Serialize to/from archive.
      *
      * \param ar input or output Archive
      * \param versionId archive format version index
      */ 
      template <class Archive>
      void serialize(Archive& ar, unsigned int versionId);

      /// \name Setters
      //@{
    
      /**
      * Set the id for this block.
      *
      * \param id integer index for this block
      */ 
      void setId(int id);
  
      /**
      * Set indices of associated vertices.
      *
      * \param vertexAId integer id of vertex A
      * \param vertexBId integer id of vertex B
      */ 
      void setVertexIds(int vertexAId, int vertexBId);
  
      /**
      * Set the monomer id.
      *
      * \param monomerId integer id of monomer type (>=0)
      */ 
      void setMonomerId(int monomerId);
  
      /**
      * Set the length of this block.
      *
      * The ``length" is steric volume / reference volume.
      *
      * \param length block length (number of monomers).
      */ 
      void setLength(double length);
  
      //@}
      /// \name Accessors (getters)
      //@{
    
      /**
      * Get the id of this block.
      */ 
      int id() const;
  
      /**
      * Get the monomer type id.
      */ 
      int monomerId() const;
  
      /**
      * Get the pair of associated vertex ids.
      */ 
      const Pair<int>& vertexIds() const;
  
      /**
      * Get id of an associated vertex.
      *
      * \param i index of vertex (0 or 1)
      */ 
      int vertexId(int i) const;
  
      /**
      * Get the length (number of monomers) in this block.
      */
      double length() const;

      //@}

   private:
  
      /// Identifier for this block, unique within the polymer.
      int id_;

      /// Identifier for the associated monomer type.
      int monomerId_;

      /// Two propagators, one for each direction. 
      Pair<int> vertexIds_;

      /// Length of this block = volume / monomer reference volume. 
      double length_;

      friend 
      std::istream& operator >> (std::istream& in, BlockDescriptor &block);

      friend 
      std::ostream& operator << (std::ostream& out, const BlockDescriptor &block);

   };

   /**
   * istream extractor for a BlockDescriptor.
   *
   * \param in  input stream
   * \param block  BlockDescriptor to be read from stream
   * \return modified input stream
   */
   std::istream& operator >> (std::istream& in, BlockDescriptor &block);

   /**
   * ostream inserter for a BlockDescriptor.
   *
   * \param out  output stream
   * \param block  BlockDescriptor to be written to stream
   * \return modified output stream
   */
   std::ostream& operator << (std::ostream& out, const BlockDescriptor &block);

   // Inline member functions

   /*
   * Get the id of this block.
   */ 
   inline int BlockDescriptor::id() const
   {  return id_; }

   /*
   * Get the monomer type id.
   */ 
   inline int BlockDescriptor::monomerId() const
   {  return monomerId_; }

   /*
   * Get the pair of associated vertex ids.
   */ 
   inline const Pair<int>& BlockDescriptor::vertexIds() const
   {  return vertexIds_; }

   /*
   * Get id of an associated vertex.
   */ 
   inline int BlockDescriptor::vertexId(int i) const
   {  return vertexIds_[i]; }

   /*
   * Get the length (number of monomers) in this block.
   */
   inline double BlockDescriptor::length() const
   {  return length_; }
    
   /*
   * Serialize to/from an archive.
   */
   template <class Archive>
   void BlockDescriptor::serialize(Archive& ar, unsigned int)
   {
      ar & id_;
      ar & monomerId_;
      ar & vertexIds_;
      ar & length_;
   }
    
} 
#endif 
